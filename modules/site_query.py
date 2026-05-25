# -*- coding: utf-8 -*-

import csv
import re
import subprocess
from pathlib import Path
from collections.abc import Iterable

from modules.annotate_mut import FIELDNAMES, annotate_candidate_row
from modules.bam_allele_freq import BASES, _count_pileup_events, _run_mpileup, ensure_bam_index
from modules.gene_database import (
    cds_index_to_genomic_pos,
    load_gene_database,
    reverse_complement,
    translate_dna,
)
from modules.utils import ensure_dir, resolve_project_path


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def _resolve_query_positions(gene_record: dict,
                             genomic_pos: int | None,
                             cds_pos: int | None,
                             aa_pos: int | None):
    provided = [value is not None for value in [genomic_pos, cds_pos, aa_pos]].count(True)
    if provided != 1:
        raise ValueError("Provide exactly one of --genomic_pos, --cds_pos, or --aa_pos.")

    chrom = gene_record.get("resolved_chrom") or gene_record["chrom"]
    if genomic_pos is not None:
        return [{
            "chrom": chrom,
            "pos": genomic_pos,
            "query_type": "genomic_pos",
            "query_label": str(genomic_pos),
        }]

    if cds_pos is not None:
        if cds_pos < 1:
            raise ValueError("--cds_pos must be 1-based and >= 1")
        pos = cds_index_to_genomic_pos(gene_record, cds_pos - 1)
        if pos is None:
            raise ValueError(f"CDS position is outside {gene_record['symbol']}: {cds_pos}")
        return [{
            "chrom": chrom,
            "pos": pos,
            "query_type": "cds_pos",
            "query_label": str(cds_pos),
        }]

    if aa_pos < 1:
        raise ValueError("--aa_pos must be 1-based and >= 1")

    positions = []
    first_cds_index = (aa_pos - 1) * 3
    for codon_offset in range(3):
        cds_index = first_cds_index + codon_offset
        pos = cds_index_to_genomic_pos(gene_record, cds_index)
        if pos is None:
            raise ValueError(f"AA position is outside {gene_record['symbol']}: {aa_pos}")
        positions.append({
            "chrom": chrom,
            "pos": pos,
            "query_type": "aa_pos",
            "query_label": f"{aa_pos}:codon_base_{codon_offset + 1}",
        })
    return positions


def _as_position_list(value):
    if value is None:
        return []
    if isinstance(value, str):
        return [int(item.strip()) for item in value.split(",") if item.strip()]
    if isinstance(value, Iterable):
        return list(value)
    return [value]


def _resolve_many_query_positions(gene_record: dict,
                                  genomic_pos: int | list[int] | None,
                                  cds_pos: int | list[int] | None,
                                  aa_pos: int | list[int] | None):
    query_groups = [
        ("genomic_pos", _as_position_list(genomic_pos)),
        ("cds_pos", _as_position_list(cds_pos)),
        ("aa_pos", _as_position_list(aa_pos)),
    ]
    provided = [(name, values) for name, values in query_groups if values]
    if len(provided) != 1:
        raise ValueError("Provide exactly one of --genomic_pos, --cds_pos, or --aa_pos.")

    query_type, values = provided[0]
    queries = []
    seen = set()
    for value in values:
        value = int(value)
        if value < 1:
            raise ValueError(f"{query_type} values must be 1-based and >= 1")
        for query in _resolve_query_positions(
            gene_record,
            genomic_pos=value if query_type == "genomic_pos" else None,
            cds_pos=value if query_type == "cds_pos" else None,
            aa_pos=value if query_type == "aa_pos" else None,
        ):
            key = (query["chrom"], query["pos"], query["query_type"], query["query_label"])
            if key in seen:
                continue
            seen.add(key)
            queries.append(query)
    return queries


def _resolve_aa_codon_queries(gene_record: dict, aa_positions: list[int]):
    chrom = gene_record.get("resolved_chrom") or gene_record["chrom"]
    cds_sequence = gene_record["cds_sequence"]
    queries = []
    seen = set()
    for aa_pos in aa_positions:
        aa_pos = int(aa_pos)
        if aa_pos < 1:
            raise ValueError("--aa_pos must be 1-based and >= 1")

        codon_start = (aa_pos - 1) * 3
        ref_codon = cds_sequence[codon_start:codon_start + 3].upper()
        if len(ref_codon) != 3:
            raise ValueError(f"AA position is outside {gene_record['symbol']}: {aa_pos}")

        genomic_positions = []
        for codon_offset in range(3):
            cds_index = codon_start + codon_offset
            pos = cds_index_to_genomic_pos(gene_record, cds_index)
            if pos is None:
                raise ValueError(f"AA position is outside {gene_record['symbol']}: {aa_pos}")
            genomic_positions.append(pos)

        key = tuple(genomic_positions)
        if key in seen:
            continue
        seen.add(key)
        queries.append({
            "chrom": chrom,
            "aa_pos": aa_pos,
            "cds_start": codon_start + 1,
            "genomic_positions": genomic_positions,
            "ref_codon": ref_codon,
            "ref_aa": translate_dna(ref_codon),
        })
    return queries


def _run_samtools_view(bam_path: Path,
                       chrom: str,
                       start: int,
                       end: int,
                       min_mapq: int,
                       logger):
    region = f"{chrom}:{start}-{end}"
    cmd = [
        "bash", "-lc",
        f"samtools view -q {min_mapq} -F 2308 {bam_path} {region}"
    ]
    logger.info(f"CMD: {' '.join(map(str, cmd))}")
    try:
        completed = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        if exc.stdout:
            logger.error(exc.stdout.strip())
        if exc.stderr:
            logger.error(exc.stderr.strip())
        raise RuntimeError(
            "samtools view failed. Check that the BAM has an index and the queried contig exists in the BAM."
        ) from exc
    return completed.stdout.splitlines()


def _read_bases_at_positions(sam_line: str, target_positions: set[int]):
    parts = sam_line.rstrip("\n").split("\t")
    if len(parts) < 11:
        return {}

    ref_pos = int(parts[3])
    cigar = parts[5]
    sequence = parts[9].upper()
    read_pos = 0
    bases = {}

    for length_text, op in CIGAR_RE.findall(cigar):
        length = int(length_text)
        if op in ("M", "=", "X"):
            for offset in range(length):
                current_ref_pos = ref_pos + offset
                if current_ref_pos in target_positions:
                    base_index = read_pos + offset
                    if base_index < len(sequence):
                        bases[current_ref_pos] = sequence[base_index]
            ref_pos += length
            read_pos += length
        elif op in ("I", "S"):
            read_pos += length
        elif op in ("D", "N"):
            ref_pos += length
        elif op in ("H", "P"):
            continue

        if target_positions.issubset(bases.keys()):
            break

    return bases


def _read_events_at_positions(sam_line: str, target_positions: set[int]):
    parts = sam_line.rstrip("\n").split("\t")
    if len(parts) < 11:
        return {"bases": {}, "insertions": {}}

    ref_pos = int(parts[3])
    cigar = parts[5]
    sequence = parts[9].upper()
    read_pos = 0
    bases = {}
    insertions = {}

    for length_text, op in CIGAR_RE.findall(cigar):
        length = int(length_text)
        if op in ("M", "=", "X"):
            for offset in range(length):
                current_ref_pos = ref_pos + offset
                if current_ref_pos in target_positions:
                    base_index = read_pos + offset
                    if base_index < len(sequence):
                        bases[current_ref_pos] = sequence[base_index]
            ref_pos += length
            read_pos += length
        elif op == "I":
            anchor_pos = ref_pos - 1
            inserted_sequence = sequence[read_pos:read_pos + length]
            if anchor_pos in target_positions and inserted_sequence:
                insertions.setdefault(anchor_pos, []).append(inserted_sequence)
            read_pos += length
        elif op == "S":
            read_pos += length
        elif op in ("D", "N"):
            ref_pos += length
        elif op in ("H", "P"):
            continue

    return {"bases": bases, "insertions": insertions}


def _coding_base_from_genomic_base(base: str, strand: str):
    base = base.upper()
    if strand == "-":
        return reverse_complement(base)
    return base


def _codon_rows_for_aa_queries(gene: str,
                               gene_record: dict,
                               codon_queries: list[dict],
                               sam_lines: list[str],
                               requested_alt: str | None,
                               min_alt_freq: float):
    target_positions = {
        pos
        for query in codon_queries
        for pos in query["genomic_positions"]
    }
    read_events = [
        _read_events_at_positions(line, target_positions)
        for line in sam_lines
    ]
    rows = []

    for query in codon_queries:
        codon_counts = {}
        insertion_counts = {}
        for read_event in read_events:
            position_bases = read_event["bases"]
            if not all(pos in position_bases for pos in query["genomic_positions"]):
                continue
            codon = "".join(
                _coding_base_from_genomic_base(position_bases[pos], gene_record["strand"])
                for pos in query["genomic_positions"]
            )
            if set(codon).issubset(set(BASES)):
                codon_counts[codon] = codon_counts.get(codon, 0) + 1

            for anchor_pos in query["genomic_positions"]:
                for inserted_sequence in read_event["insertions"].get(anchor_pos, []):
                    inserted_coding = _coding_base_from_genomic_base(inserted_sequence, gene_record["strand"])
                    if set(inserted_coding).issubset(set(BASES)):
                        insertion_counts[inserted_coding] = insertion_counts.get(inserted_coding, 0) + 1

        depth = sum(codon_counts.values())
        ref_codon = query["ref_codon"]
        ref_depth = codon_counts.get(ref_codon, 0)

        if requested_alt:
            requested = requested_alt.upper()
            if len(requested) == 3 and set(requested).issubset(set(BASES)):
                codons = [requested]
            else:
                codons = [
                    codon for codon in codon_counts
                    if translate_dna(codon) == requested
                ]
                if not codons:
                    codons = [ref_codon] if depth == 0 else []
        else:
            codons = []
            if ref_codon:
                codons.append(ref_codon)
            for codon, count in sorted(codon_counts.items()):
                if codon == ref_codon:
                    continue
                allele_freq = count / depth if depth else 0
                if allele_freq >= min_alt_freq:
                    codons.append(codon)

        pos_label = "-".join(str(pos) for pos in query["genomic_positions"])
        for codon in codons:
            alt_depth = codon_counts.get(codon, 0)
            allele_freq = alt_depth / depth if depth else 0
            alt_aa = translate_dna(codon)
            ref_aa = query["ref_aa"]
            aa_pos = query["aa_pos"]
            aa_change = f"{ref_aa}{aa_pos}{alt_aa}"
            if ref_aa == alt_aa:
                aa_change = f"{ref_aa}{aa_pos}{ref_aa} (synonymous)"

            note = f"site_query:aa_pos={aa_pos}; codon-spanning reads only"
            if depth == 0:
                note += "; no reads mapped across requested codon"
            elif codon == ref_codon:
                note += "; reference codon"
            elif alt_depth == 0:
                note += "; requested codon/amino acid not observed"

            rows.append({
                "Gene name": gene_record.get("description") or gene,
                "Gene symbol": gene,
                "gene": gene,
                "chrom": query["chrom"],
                "pos": pos_label,
                "ref": ref_codon,
                "alt": codon,
                "qual": "",
                "filter": "PASS" if depth else "NO_MAPPING",
                "depth": depth,
                "ref_depth": ref_depth,
                "alt_depth": alt_depth,
                "allele_freq": f"{allele_freq:.6f}",
                "variant_type": "REF" if codon == ref_codon else "CODON",
                "source": "bam_codon_site_query",
                "cds_pos": query["cds_start"],
                "codon_pos": "1-3",
                "ref_codon": ref_codon,
                "alt_codon": codon,
                "aa_pos": aa_pos,
                "ref_aa": ref_aa,
                "alt_aa": alt_aa,
                "aa_change": aa_change,
                "note": note,
            })

        if requested_alt and "+" in requested_alt:
            requested_insertions = [requested_alt.split("+", 1)[1].upper()]
        elif requested_alt and len(requested_alt) != 3:
            requested_insertions = []
        else:
            requested_insertions = [
                sequence for sequence, count in sorted(insertion_counts.items())
                if depth and count / depth >= min_alt_freq
            ]

        for inserted_sequence in requested_insertions:
            if not inserted_sequence:
                continue
            alt_depth = insertion_counts.get(inserted_sequence, 0)
            allele_freq = alt_depth / depth if depth else 0
            if len(inserted_sequence) % 3 == 0:
                inserted_aa = translate_dna(inserted_sequence)
                aa_change = f"{query['ref_aa']}{query['aa_pos']}_ins{inserted_aa}"
                alt_aa = inserted_aa
                note = f"site_query:aa_pos={query['aa_pos']}; codon-spanning reads only; in-frame insertion"
            else:
                aa_change = f"{query['ref_aa']}{query['aa_pos']}fs"
                alt_aa = "frameshift"
                note = f"site_query:aa_pos={query['aa_pos']}; codon-spanning reads only; frameshift insertion"
            if depth == 0:
                note += "; no reads mapped across requested codon"
            elif alt_depth == 0:
                note += "; requested insertion not observed"

            rows.append({
                "Gene name": gene_record.get("description") or gene,
                "Gene symbol": gene,
                "gene": gene,
                "chrom": query["chrom"],
                "pos": pos_label,
                "ref": ref_codon,
                "alt": f"{ref_codon}+{inserted_sequence}",
                "qual": "",
                "filter": "PASS" if depth else "NO_MAPPING",
                "depth": depth,
                "ref_depth": ref_depth,
                "alt_depth": alt_depth,
                "allele_freq": f"{allele_freq:.6f}",
                "variant_type": "INS",
                "source": "bam_codon_site_query",
                "cds_pos": query["cds_start"],
                "codon_pos": "1-3",
                "ref_codon": ref_codon,
                "alt_codon": f"{ref_codon}+{inserted_sequence}",
                "aa_pos": query["aa_pos"],
                "ref_aa": query["ref_aa"],
                "alt_aa": alt_aa,
                "aa_change": aa_change,
                "note": note,
            })

    return rows


def _pileup_line_by_position(lines):
    by_pos = {}
    for line in lines:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 5:
            by_pos[int(parts[1])] = parts
    return by_pos


def _rows_for_site(gene: str,
                   gene_record: dict,
                   gene_db: dict,
                   query: dict,
                   pileup_parts,
                   requested_alt: str | None,
                   min_alt_freq: float):
    if pileup_parts:
        chrom, pos, ref_base, raw_depth, read_bases = pileup_parts[:5]
        events = _count_pileup_events(read_bases, ref_base)
        counts = events["bases"]
        insertions = events["insertions"]
    else:
        chrom = query["chrom"]
        pos = str(query["pos"])
        ref_base = ""
        counts = {base: 0 for base in BASES}
        insertions = {}

    ref = ref_base.upper()
    depth = sum(counts.values())
    ref_depth = counts.get(ref, 0) if ref else 0

    requested_alt = requested_alt.upper() if requested_alt else None
    if requested_alt:
        alleles = [requested_alt.upper()]
    else:
        observed = []
        for base in BASES:
            allele_depth = counts[base]
            allele_freq = allele_depth / depth if depth else 0
            if base == ref or allele_freq >= min_alt_freq:
                if allele_depth > 0:
                    observed.append(base)
        alleles = observed or ([ref] if ref else [])
        if ref and ref not in alleles:
            alleles.insert(0, ref)
        if not alleles:
            alleles = [ref]

    rows = []
    for allele in alleles:
        if len(allele) > 1 and ref and allele.startswith(ref):
            continue
        allele_depth = counts.get(allele, 0)
        allele_freq = allele_depth / depth if depth else 0
        variant_type = "REF" if allele == ref else "SNV"
        note = f"site_query:{query['query_type']}={query['query_label']}"
        if depth == 0:
            note += "; no reads mapped at requested site"
        elif allele == ref:
            note += "; reference allele"
        elif allele_depth == 0:
            note += "; requested allele not observed"

        row = {
            "Gene name": gene_record.get("description") or gene,
            "Gene symbol": gene,
            "gene": gene,
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": allele,
            "qual": "",
            "filter": "PASS" if depth else "NO_MAPPING",
            "depth": depth,
            "ref_depth": ref_depth,
            "alt_depth": allele_depth,
            "allele_freq": f"{allele_freq:.6f}",
            "variant_type": variant_type,
            "source": "bam_site_query",
            "note": note,
        }
        rows.append(annotate_candidate_row(row, gene_db))

    insertion_alleles = []
    if requested_alt and ref and len(requested_alt) > 1:
        inserted_sequence = requested_alt[len(ref):] if requested_alt.startswith(ref) else requested_alt
        insertion_alleles = [inserted_sequence]
    elif not requested_alt:
        insertion_alleles = [
            sequence for sequence, count in sorted(insertions.items())
            if depth and count / depth >= min_alt_freq
        ]

    for inserted_sequence in insertion_alleles:
        if not inserted_sequence:
            continue
        allele_depth = insertions.get(inserted_sequence, 0)
        allele_freq = allele_depth / depth if depth else 0
        note = f"site_query:{query['query_type']}={query['query_label']}; insertion_sequence={inserted_sequence}"
        if depth == 0:
            note += "; no reads mapped at requested site"
        elif allele_depth == 0:
            note += "; requested insertion not observed"

        row = {
            "Gene name": gene_record.get("description") or gene,
            "Gene symbol": gene,
            "gene": gene,
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": f"{ref}{inserted_sequence}" if ref else inserted_sequence,
            "qual": "",
            "filter": "PASS" if depth else "NO_MAPPING",
            "depth": depth,
            "ref_depth": ref_depth,
            "alt_depth": allele_depth,
            "allele_freq": f"{allele_freq:.6f}",
            "variant_type": "INS",
            "source": "bam_site_query",
            "note": note,
        }
        rows.append(annotate_candidate_row(row, gene_db))
    return rows


def query_bam_sites(gene_db_path: Path,
                    bam_path: Path,
                    out_csv: Path,
                    gene: str,
                    genomic_pos: int | list[int] | None,
                    cds_pos: int | list[int] | None,
                    aa_pos: int | list[int] | None,
                   alt: str | None,
                   ref_fasta: Path | None,
                   min_mapq: int,
                   min_alt_freq: float,
                   force: bool,
                   logger):
    if out_csv.exists() and not force:
        logger.info(f"Site query table exists, skip: {out_csv}")
        return out_csv

    gene_db = load_gene_database(gene_db_path)
    gene_record = gene_db.get("genes_by_symbol", {}).get(gene)
    if not gene_record:
        raise ValueError(f"Gene not found in database: {gene}")

    if not bam_path.exists():
        raise FileNotFoundError(f"BAM not found: {bam_path}")
    ensure_bam_index(bam_path, logger)

    if ref_fasta is None:
        ref_fasta = resolve_project_path(gene_db["reference_fasta"])
    if not ref_fasta.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {ref_fasta}")

    provided_query_types = [
        bool(_as_position_list(value))
        for value in (genomic_pos, cds_pos, aa_pos)
    ].count(True)
    if provided_query_types != 1:
        raise ValueError("Provide exactly one of --genomic_pos, --cds_pos, or --aa_pos.")

    if _as_position_list(aa_pos):
        codon_queries = _resolve_aa_codon_queries(gene_record, _as_position_list(aa_pos))
        all_positions = [
            pos
            for query in codon_queries
            for pos in query["genomic_positions"]
        ]
        chrom = codon_queries[0]["chrom"]
        sam_lines = _run_samtools_view(
            bam_path=bam_path,
            chrom=chrom,
            start=min(all_positions),
            end=max(all_positions),
            min_mapq=min_mapq,
            logger=logger,
        )
        rows = _codon_rows_for_aa_queries(
            gene=gene,
            gene_record=gene_record,
            codon_queries=codon_queries,
            sam_lines=sam_lines,
            requested_alt=alt,
            min_alt_freq=min_alt_freq,
        )
    else:
        queries = _resolve_many_query_positions(gene_record, genomic_pos, cds_pos, aa_pos)
        start = min(query["pos"] for query in queries)
        end = max(query["pos"] for query in queries)
        chrom = queries[0]["chrom"]
        lines = _run_mpileup(
            bam_path=bam_path,
            ref_fasta=ref_fasta,
            chrom=chrom,
            start=start,
            end=end,
            min_mapq=min_mapq,
            logger=logger,
        )
        pileup_by_pos = _pileup_line_by_position(lines)

        rows = []
        for query in queries:
            rows.extend(_rows_for_site(
                gene=gene,
                gene_record=gene_record,
                gene_db=gene_db,
                query=query,
                pileup_parts=pileup_by_pos.get(query["pos"]),
                requested_alt=alt,
                min_alt_freq=min_alt_freq,
            ))

    ensure_dir(out_csv.parent)
    with open(out_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    logger.info(f"Wrote {len(rows)} queried site rows: {out_csv}")
    return out_csv
