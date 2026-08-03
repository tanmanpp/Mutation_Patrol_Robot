# -*- coding: utf-8 -*-

import csv
import json
import re
from pathlib import Path
from collections.abc import Iterable

from modules.annotate_mut import AA_NAMES, FIELDNAMES, annotate_candidate_row
from modules.bam_allele_freq import BASES, _count_pileup_events, _run_mpileup, ensure_bam_index
from modules.gene_database import (
    cds_index_to_genomic_pos,
    load_gene_database,
    reverse_complement,
    translate_dna,
)
from modules.haplotype_translate import (
    build_protein_haplotypes,
    write_protein_haplotypes,
)
from modules.utils import ensure_dir, resolve_project_path, run_cmd


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")

NO_CALL_REGION_FIELDS = [
    "gene", "chrom", "start", "end", "length", "reason",
    "mean_depth", "max_depth", "minimum_call_depth",
]


def _optional_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _biological_position_sort_key(row: dict):
    """Sort coding results from the protein N-terminus to C-terminus."""
    aa_pos = _optional_int(row.get("aa_pos"))
    cds_pos = _optional_int(row.get("cds_pos"))
    genomic_pos = _optional_int(row.get("pos"))
    if aa_pos is not None:
        return (
            str(row.get("gene") or row.get("Gene symbol") or ""),
            0,
            aa_pos,
            cds_pos if cds_pos is not None else float("inf"),
            genomic_pos if genomic_pos is not None else float("inf"),
            str(row.get("alt") or row.get("alt_codon") or ""),
        )
    return (
        str(row.get("gene") or row.get("Gene symbol") or ""),
        1,
        genomic_pos if genomic_pos is not None else float("inf"),
        float("inf"),
        float("inf"),
        str(row.get("alt") or row.get("alt_codon") or ""),
    )


def _site_filter(depth: int, min_depth: int):
    if depth == 0:
        return "NO_MAPPING"
    if depth < min_depth:
        return "LOW_DEPTH"
    return "PASS"


def _qualifying_alleles(counts: dict[str, int],
                        depth: int,
                        min_allele_count: int,
                        min_allele_freq: float):
    if depth <= 0:
        return []
    return [
        allele
        for allele, count in counts.items()
        if count >= min_allele_count and count / depth >= min_allele_freq
    ]


def _site_call_status(depth: int,
                      min_depth: int,
                      counts: dict[str, int],
                      ref: str,
                      min_allele_count: int,
                      min_allele_freq: float):
    if depth < min_depth:
        return "NO_CALL"
    qualifying = _qualifying_alleles(
        counts, depth, min_allele_count, min_allele_freq
    )
    if len(qualifying) >= 2:
        return "MIXED_SIGNAL"
    dominant = max(counts, key=counts.get) if counts else ""
    if dominant == ref:
        return "REFERENCE"
    return "VARIANT"


def _allele_status(depth: int, min_depth: int, allele_depth: int):
    if depth < min_depth:
        return "NO_CALL"
    return "OBSERVED" if allele_depth > 0 else "NOT_OBSERVED"


def _allele_spectrum(counts: dict[str, int], depth: int):
    if depth <= 0:
        return ""
    observed = [
        (allele, count)
        for allele, count in counts.items()
        if count > 0
    ]
    observed.sort(key=lambda item: (-item[1], item[0]))
    return "; ".join(
        f"{allele}:{count} ({count / depth:.1%})"
        for allele, count in observed
    )


def _site_qc_flags(depth: int,
                   min_depth: int,
                   call_status: str,
                   allele_depth: int | None = None,
                   min_allele_count: int = 1):
    flags = []
    if depth == 0:
        flags.append("NO_MAPPING")
    elif depth < min_depth:
        flags.append("LOW_DEPTH")
    if call_status == "MIXED_SIGNAL":
        flags.append("MULTIPLE_ALLELES")
    if (
        allele_depth is not None
        and 0 < allele_depth < min_allele_count
        and depth >= min_depth
    ):
        flags.append("LOW_ALLELE_SUPPORT")
    return ";".join(flags) or "PASS"


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
    try:
        completed = run_cmd(
            [
                "samtools", "view", "-q", str(min_mapq), "-F", "2308",
                str(bam_path), region,
            ],
            logger=logger,
            capture_output=True,
        )
    except Exception as exc:
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


def _read_events_at_positions(sam_line: str,
                              target_positions: set[int],
                              min_baseq: int = 20):
    parts = sam_line.rstrip("\n").split("\t")
    if len(parts) < 11:
        return {"bases": {}, "insertions": {}, "strand": "unknown"}

    flag = int(parts[1])
    ref_pos = int(parts[3])
    cigar = parts[5]
    sequence = parts[9].upper()
    qualities = parts[10]
    strand = "reverse" if flag & 16 else "forward"
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
                    base_quality = (
                        ord(qualities[base_index]) - 33
                        if base_index < len(qualities) and qualities != "*"
                        else 0
                    )
                    if base_index < len(sequence) and base_quality >= min_baseq:
                        bases[current_ref_pos] = sequence[base_index]
            ref_pos += length
            read_pos += length
        elif op == "I":
            anchor_pos = ref_pos - 1
            inserted_sequence = sequence[read_pos:read_pos + length]
            inserted_qualities = qualities[read_pos:read_pos + length]
            quality_ok = (
                qualities != "*"
                and len(inserted_qualities) == length
                and all(ord(value) - 33 >= min_baseq for value in inserted_qualities)
            )
            if anchor_pos in target_positions and inserted_sequence and quality_ok:
                insertions.setdefault(anchor_pos, []).append(inserted_sequence)
            read_pos += length
        elif op == "S":
            read_pos += length
        elif op in ("D", "N"):
            ref_pos += length
        elif op in ("H", "P"):
            continue

    return {"bases": bases, "insertions": insertions, "strand": strand}


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
                               min_alt_freq: float,
                               min_depth: int,
                               min_mapq: int,
                               min_baseq: int,
                               min_allele_count: int):
    target_positions = {
        pos
        for query in codon_queries
        for pos in query["genomic_positions"]
    }
    read_events = [
        _read_events_at_positions(line, target_positions, min_baseq)
        for line in sam_lines
    ]
    rows = []

    for query in codon_queries:
        codon_counts = {}
        codon_strand_counts = {}
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
                strands = codon_strand_counts.setdefault(
                    codon, {"forward": 0, "reverse": 0}
                )
                strands[read_event["strand"]] += 1

            for anchor_pos in query["genomic_positions"]:
                for inserted_sequence in read_event["insertions"].get(anchor_pos, []):
                    inserted_coding = _coding_base_from_genomic_base(inserted_sequence, gene_record["strand"])
                    if set(inserted_coding).issubset(set(BASES)):
                        insertion_counts[inserted_coding] = insertion_counts.get(inserted_coding, 0) + 1

        depth = sum(codon_counts.values())
        ref_codon = query["ref_codon"]
        ref_depth = codon_counts.get(ref_codon, 0)
        call_status = _site_call_status(
            depth,
            min_depth,
            codon_counts,
            ref_codon,
            min_allele_count,
            min_alt_freq,
        )
        spectrum = _allele_spectrum(codon_counts, depth)
        allele_count = sum(1 for count in codon_counts.values() if count > 0)

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
                if count >= min_allele_count and allele_freq >= min_alt_freq:
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
                note += "; no reads mapped across requested codon; unable to determine"
            elif depth < min_depth:
                note += f"; depth below minimum {min_depth}; unable to determine"
            elif codon == ref_codon:
                note += "; reference codon"
            elif alt_depth == 0:
                note += "; requested codon/amino acid not observed"

            if call_status == "NO_CALL":
                aa_change = "無法判斷"
                effect = "no_call"
            elif codon == ref_codon:
                effect = "reference"
            elif ref_aa == alt_aa:
                effect = "synonymous"
            elif alt_aa == "*":
                effect = "stop_gained"
            elif ref_aa == "*":
                effect = "stop_lost"
            else:
                effect = "missense"
            strand_counts = codon_strand_counts.get(
                codon, {"forward": 0, "reverse": 0}
            )

            rows.append({
                "Gene name": gene_record.get("description") or gene,
                "Gene symbol": gene,
                "gene": gene,
                "chrom": query["chrom"],
                "pos": pos_label,
                "ref": ref_codon,
                "alt": codon,
                "qual": "",
                "filter": _site_filter(depth, min_depth),
                "depth": depth,
                "raw_depth": depth,
                "ref_depth": ref_depth,
                "alt_depth": alt_depth,
                "forward_depth": strand_counts["forward"],
                "reverse_depth": strand_counts["reverse"],
                "allele_freq": f"{allele_freq:.6f}",
                "allele_count": allele_count,
                "allele_spectrum": spectrum,
                "call_status": call_status,
                "allele_status": _allele_status(depth, min_depth, alt_depth),
                "qc_flags": _site_qc_flags(
                    depth, min_depth, call_status, alt_depth, min_allele_count
                ),
                "min_depth": min_depth,
                "min_mapq": min_mapq,
                "min_baseq": min_baseq,
                "min_allele_count": min_allele_count,
                "variant_type": (
                    "NO_CALL"
                    if call_status == "NO_CALL"
                    else ("REF" if codon == ref_codon else "CODON")
                ),
                "region_type": "CDS",
                "effect": effect,
                "source": "bam_codon_site_query",
                "cds_pos": query["cds_start"],
                "codon_pos": "1-3",
                "ref_codon": ref_codon,
                "alt_codon": codon,
                "codon_change": f"{ref_codon}>{codon}",
                "aa_pos": aa_pos,
                "ref_aa": ref_aa,
                "ref_aa_name": AA_NAMES.get(ref_aa, "Unknown"),
                "alt_aa": alt_aa,
                "alt_aa_name": AA_NAMES.get(alt_aa, "Unknown"),
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
                if (
                    depth
                    and count >= min_allele_count
                    and count / depth >= min_alt_freq
                )
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
                effect = "inframe_insertion"
            else:
                aa_change = f"{query['ref_aa']}{query['aa_pos']}fs"
                alt_aa = "frameshift"
                note = f"site_query:aa_pos={query['aa_pos']}; codon-spanning reads only; frameshift insertion"
                effect = "frameshift"
            if depth == 0:
                note += "; no reads mapped across requested codon; unable to determine"
            elif depth < min_depth:
                note += f"; depth below minimum {min_depth}; unable to determine"
            elif alt_depth == 0:
                note += "; requested insertion not observed"

            insertion_counts_for_call = {
                ref_codon: max(depth - alt_depth, 0),
                f"+{inserted_sequence}": alt_depth,
            }
            insertion_call_status = _site_call_status(
                depth,
                min_depth,
                insertion_counts_for_call,
                ref_codon,
                min_allele_count,
                min_alt_freq,
            )
            if insertion_call_status == "NO_CALL":
                aa_change = "無法判斷"
                effect = "no_call"

            rows.append({
                "Gene name": gene_record.get("description") or gene,
                "Gene symbol": gene,
                "gene": gene,
                "chrom": query["chrom"],
                "pos": pos_label,
                "ref": ref_codon,
                "alt": f"{ref_codon}+{inserted_sequence}",
                "qual": "",
                "filter": _site_filter(depth, min_depth),
                "depth": depth,
                "raw_depth": depth,
                "ref_depth": ref_depth,
                "alt_depth": alt_depth,
                "forward_depth": "",
                "reverse_depth": "",
                "allele_freq": f"{allele_freq:.6f}",
                "allele_count": allele_count,
                "allele_spectrum": spectrum,
                "call_status": insertion_call_status,
                "allele_status": _allele_status(depth, min_depth, alt_depth),
                "qc_flags": _site_qc_flags(
                    depth,
                    min_depth,
                    insertion_call_status,
                    alt_depth,
                    min_allele_count,
                ),
                "min_depth": min_depth,
                "min_mapq": min_mapq,
                "min_baseq": min_baseq,
                "min_allele_count": min_allele_count,
                "variant_type": (
                    "NO_CALL" if insertion_call_status == "NO_CALL" else "INS"
                ),
                "region_type": "CDS",
                "effect": effect,
                "source": "bam_codon_site_query",
                "cds_pos": query["cds_start"],
                "codon_pos": "1-3",
                "ref_codon": ref_codon,
                "alt_codon": f"{ref_codon}+{inserted_sequence}",
                "codon_change": f"{ref_codon}>{ref_codon}+{inserted_sequence}",
                "aa_pos": query["aa_pos"],
                "ref_aa": query["ref_aa"],
                "ref_aa_name": AA_NAMES.get(query["ref_aa"], "Unknown"),
                "alt_aa": alt_aa,
                "alt_aa_name": (
                    ", ".join(
                        AA_NAMES.get(value, "Unknown") for value in alt_aa
                    )
                    if effect == "inframe_insertion"
                    else ("Frameshift" if effect == "frameshift" else "")
                ),
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


def _compress_no_call_regions(gene: str,
                              chrom: str,
                              positions: list[tuple[int, int, str]],
                              min_depth: int):
    if not positions:
        return []
    regions = []
    current = None
    for pos, depth, reason in positions:
        if (
            current
            and pos == current["end"] + 1
            and reason == current["reason"]
        ):
            current["end"] = pos
            current["depths"].append(depth)
            continue
        if current:
            regions.append(current)
        current = {
            "gene": gene,
            "chrom": chrom,
            "start": pos,
            "end": pos,
            "reason": reason,
            "depths": [depth],
        }
    if current:
        regions.append(current)

    output = []
    for region in regions:
        depths = region.pop("depths")
        output.append({
            **region,
            "length": region["end"] - region["start"] + 1,
            "mean_depth": f"{sum(depths) / len(depths):.2f}",
            "max_depth": max(depths),
            "minimum_call_depth": min_depth,
        })
    return output


def _apply_protein_haplotype_annotations(variant_rows: list[dict],
                                         haplotype_rows: list[dict]):
    """Replace provisional site-level protein effects with phased consequences."""
    by_cluster = {}
    for row in haplotype_rows:
        cluster_id = row.get("cluster_id")
        if cluster_id:
            by_cluster.setdefault(cluster_id, []).append(row)

    for cluster_id, cluster_rows in by_cluster.items():
        resolved_rows = [
            row for row in cluster_rows
            if row.get("phase_status") == "PHASED"
        ]
        call_status = cluster_rows[0].get("call_status") or "NO_CALL"
        aa_positions = [
            value
            for row in resolved_rows
            for value in (
                _optional_int(row.get("aa_start")),
                _optional_int(row.get("aa_end")),
            )
            if value is not None
        ]
        aa_start = min(aa_positions) if aa_positions else None
        aa_end = max(aa_positions) if aa_positions else None
        genomic_start = min(
            _optional_int(row.get("cluster_start")) or 0
            for row in cluster_rows
        )
        genomic_end = max(
            _optional_int(row.get("cluster_end")) or 0
            for row in cluster_rows
        )
        dominant = max(
            resolved_rows,
            key=lambda row: float(row.get("haplotype_frequency") or 0),
            default=None,
        )

        for variant_row in variant_rows:
            if variant_row.get("region_type") != "CDS":
                continue
            row_aa_pos = _optional_int(variant_row.get("aa_pos"))
            row_genomic_pos = _optional_int(variant_row.get("pos"))
            in_aa_range = (
                aa_start is not None
                and aa_end is not None
                and row_aa_pos is not None
                and aa_start <= row_aa_pos <= aa_end
            )
            is_cluster_indel = (
                variant_row.get("variant_type") in {"INS", "DEL"}
                and row_genomic_pos is not None
                and genomic_start <= row_genomic_pos <= genomic_end
            )
            if not in_aa_range and not is_cluster_indel:
                continue

            original_change = variant_row.get("aa_change") or ""
            original_effect = variant_row.get("effect") or ""
            note = variant_row.get("note") or ""
            note += (
                f"; site_level_effect={original_effect}; "
                f"site_level_aa_change={original_change}; "
                f"see protein haplotype {cluster_id}"
            )
            variant_row["note"] = note.lstrip("; ")
            qc_flags = [
                value for value in (variant_row.get("qc_flags") or "").split(";")
                if value and value != "PASS"
            ]

            if call_status == "MIXED_SIGNAL":
                variant_row["effect"] = "haplotype_mixed"
                variant_row["aa_change"] = (
                    f"See {cluster_id} protein haplotypes (MIXED_SIGNAL)"
                )
                qc_flags.append("MULTIPLE_HAPLOTYPES")
            elif dominant is None:
                variant_row["effect"] = "no_call"
                variant_row["aa_change"] = "Unable to determine (PHASE_UNRESOLVED)"
                qc_flags.append("PHASE_UNRESOLVED")
            else:
                variant_row["effect"] = dominant.get("effect") or original_effect
                variant_row["aa_change"] = (
                    dominant.get("aa_changes")
                    or dominant.get("combined_aa_change")
                    or original_change
                )
                qc_flags.append("HAPLOTYPE_RECONSTRUCTED")
            variant_row["qc_flags"] = ";".join(dict.fromkeys(qc_flags)) or "PASS"


def _rows_for_site(gene: str,
                   gene_record: dict,
                   gene_db: dict,
                   query: dict,
                   pileup_parts,
                   requested_alt: str | None,
                   min_alt_freq: float,
                   min_depth: int,
                   min_mapq: int,
                   min_baseq: int,
                   min_allele_count: int):
    if pileup_parts:
        chrom, pos, ref_base, raw_depth, read_bases = pileup_parts[:5]
        events = _count_pileup_events(read_bases, ref_base)
        counts = events["bases"]
        strand_counts = events["strand_counts"]
        insertions = events["insertions"]
        deletions = events["deletions"]
    else:
        chrom = query["chrom"]
        pos = str(query["pos"])
        ref_base = ""
        raw_depth = "0"
        counts = {base: 0 for base in BASES}
        strand_counts = {
            base: {"forward": 0, "reverse": 0}
            for base in BASES
        }
        insertions = {}
        deletions = {}

    ref = ref_base.upper()
    depth = sum(counts.values())
    ref_depth = counts.get(ref, 0) if ref else 0
    call_status = _site_call_status(
        depth,
        min_depth,
        counts,
        ref,
        min_allele_count,
        min_alt_freq,
    )
    spectrum = _allele_spectrum(counts, depth)
    allele_count = sum(1 for count in counts.values() if count > 0)

    requested_alt = requested_alt.upper() if requested_alt else None
    if requested_alt:
        alleles = [requested_alt.upper()]
    else:
        observed = []
        for base in BASES:
            allele_depth = counts[base]
            allele_freq = allele_depth / depth if depth else 0
            if (
                base == ref
                or (
                    allele_depth >= min_allele_count
                    and allele_freq >= min_alt_freq
                )
            ):
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
            note += "; no reads mapped at requested site; unable to determine"
        elif depth < min_depth:
            note += f"; depth below minimum {min_depth}; unable to determine"
        elif allele == ref:
            note += "; reference allele"
        elif allele_depth == 0:
            note += "; requested allele not observed"

        allele_strands = strand_counts.get(
            allele, {"forward": 0, "reverse": 0}
        )
        row = {
            "Gene name": gene_record.get("description") or gene,
            "Gene symbol": gene,
            "gene": gene,
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": allele,
            "qual": "",
            "filter": _site_filter(depth, min_depth),
            "depth": depth,
            "raw_depth": raw_depth,
            "ref_depth": ref_depth,
            "alt_depth": allele_depth,
            "forward_depth": allele_strands["forward"],
            "reverse_depth": allele_strands["reverse"],
            "allele_freq": f"{allele_freq:.6f}",
            "allele_count": allele_count,
            "allele_spectrum": spectrum,
            "call_status": call_status,
            "allele_status": _allele_status(depth, min_depth, allele_depth),
            "qc_flags": _site_qc_flags(
                depth, min_depth, call_status, allele_depth, min_allele_count
            ),
            "min_depth": min_depth,
            "min_mapq": min_mapq,
            "min_baseq": min_baseq,
            "min_allele_count": min_allele_count,
            "variant_type": "NO_CALL" if call_status == "NO_CALL" else variant_type,
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
            if (
                depth
                and count >= min_allele_count
                and count / depth >= min_alt_freq
            )
        ]

    for inserted_sequence in insertion_alleles:
        if not inserted_sequence:
            continue
        allele_depth = insertions.get(inserted_sequence, 0)
        allele_freq = allele_depth / depth if depth else 0
        note = f"site_query:{query['query_type']}={query['query_label']}; insertion_sequence={inserted_sequence}"
        if depth == 0:
            note += "; no reads mapped at requested site; unable to determine"
        elif depth < min_depth:
            note += f"; depth below minimum {min_depth}; unable to determine"
        elif allele_depth == 0:
            note += "; requested insertion not observed"

        insertion_counts_for_call = {
            ref: max(depth - allele_depth, 0),
            f"+{inserted_sequence}": allele_depth,
        }
        insertion_call_status = _site_call_status(
            depth,
            min_depth,
            insertion_counts_for_call,
            ref,
            min_allele_count,
            min_alt_freq,
        )
        row = {
            "Gene name": gene_record.get("description") or gene,
            "Gene symbol": gene,
            "gene": gene,
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": f"{ref}{inserted_sequence}" if ref else inserted_sequence,
            "qual": "",
            "filter": _site_filter(depth, min_depth),
            "depth": depth,
            "raw_depth": raw_depth,
            "ref_depth": ref_depth,
            "alt_depth": allele_depth,
            "forward_depth": "",
            "reverse_depth": "",
            "allele_freq": f"{allele_freq:.6f}",
            "allele_count": allele_count,
            "allele_spectrum": spectrum,
            "call_status": insertion_call_status,
            "allele_status": _allele_status(depth, min_depth, allele_depth),
            "qc_flags": _site_qc_flags(
                depth,
                min_depth,
                insertion_call_status,
                allele_depth,
                min_allele_count,
            ),
            "min_depth": min_depth,
            "min_mapq": min_mapq,
            "min_baseq": min_baseq,
            "min_allele_count": min_allele_count,
            "variant_type": (
                "NO_CALL" if insertion_call_status == "NO_CALL" else "INS"
            ),
            "source": "bam_site_query",
            "note": note,
        }
        rows.append(annotate_candidate_row(row, gene_db))

    if not requested_alt:
        deletion_alleles = [
            sequence for sequence, count in sorted(deletions.items())
            if (
                depth
                and count >= min_allele_count
                and count / depth >= min_alt_freq
            )
        ]
    else:
        deletion_alleles = []

    for deleted_sequence in deletion_alleles:
        allele_depth = deletions.get(deleted_sequence, 0)
        allele_freq = allele_depth / depth if depth else 0
        deletion_label = f"-{deleted_sequence}"
        deletion_counts_for_call = {
            ref: max(depth - allele_depth, 0),
            deletion_label: allele_depth,
        }
        deletion_call_status = _site_call_status(
            depth,
            min_depth,
            deletion_counts_for_call,
            ref,
            min_allele_count,
            min_alt_freq,
        )
        note = (
            f"site_query:{query['query_type']}={query['query_label']}; "
            f"deletion_sequence={deleted_sequence}"
        )
        row = {
            "Gene name": gene_record.get("description") or gene,
            "Gene symbol": gene,
            "gene": gene,
            "chrom": chrom,
            "pos": pos,
            "ref": f"{ref}{deleted_sequence}" if ref else deleted_sequence,
            "alt": ref,
            "qual": "",
            "filter": _site_filter(depth, min_depth),
            "depth": depth,
            "raw_depth": raw_depth,
            "ref_depth": ref_depth,
            "alt_depth": allele_depth,
            "forward_depth": "",
            "reverse_depth": "",
            "allele_freq": f"{allele_freq:.6f}",
            "allele_count": sum(
                count > 0 for count in deletion_counts_for_call.values()
            ),
            "allele_spectrum": _allele_spectrum(
                deletion_counts_for_call, depth
            ),
            "call_status": deletion_call_status,
            "allele_status": _allele_status(depth, min_depth, allele_depth),
            "qc_flags": _site_qc_flags(
                depth,
                min_depth,
                deletion_call_status,
                allele_depth,
                min_allele_count,
            ),
            "min_depth": min_depth,
            "min_mapq": min_mapq,
            "min_baseq": min_baseq,
            "min_allele_count": min_allele_count,
            "variant_type": (
                "NO_CALL" if deletion_call_status == "NO_CALL" else "DEL"
            ),
            "source": "bam_site_query",
            "note": note,
        }
        rows.append(annotate_candidate_row(row, gene_db))
    return rows


def scan_bam_gene_region(gene_db_path: Path,
                         bam_path: Path,
                         out_csv: Path,
                         gene: str,
                         ref_fasta: Path | None,
                         min_mapq: int,
                         min_alt_freq: float,
                         min_depth: int,
                         min_baseq: int,
                         min_allele_count: int,
                         force: bool,
                         logger):
    summary_path = out_csv.parent / "scan_summary.json"
    no_call_csv = out_csv.parent / "no_call_regions.csv"
    complete_table_csv = out_csv.parent / "complete_gene_table.csv"
    protein_haplotype_csv = out_csv.parent / "protein_haplotypes.csv"
    if (
        out_csv.exists()
        and summary_path.exists()
        and no_call_csv.exists()
        and complete_table_csv.exists()
        and protein_haplotype_csv.exists()
        and not force
    ):
        logger.info(f"Whole-gene scan exists, skip: {out_csv}")
        return {
            "table": out_csv,
            "summary": summary_path,
            "no_call_regions": no_call_csv,
            "complete_table": complete_table_csv,
            "protein_haplotypes": protein_haplotype_csv,
        }

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

    chrom = gene_record.get("resolved_chrom") or gene_record["chrom"]
    start = int(gene_record["gene_start"])
    end = int(gene_record["gene_end"])
    lines = _run_mpileup(
        bam_path=bam_path,
        ref_fasta=ref_fasta,
        chrom=chrom,
        start=start,
        end=end,
        min_mapq=min_mapq,
        logger=logger,
        min_baseq=min_baseq,
    )
    pileup_by_pos = _pileup_line_by_position(lines)

    variant_rows = []
    complete_rows = []
    no_call_positions = []
    callable_positions = 0
    reference_positions = 0
    mixed_positions = set()
    variant_positions = set()
    variant_type_counts = {}

    for pos in range(start, end + 1):
        query = {
            "chrom": chrom,
            "pos": pos,
            "query_type": "gene_region",
            "query_label": f"{start}-{end}",
        }
        rows = _rows_for_site(
            gene=gene,
            gene_record=gene_record,
            gene_db=gene_db,
            query=query,
            pileup_parts=pileup_by_pos.get(pos),
            requested_alt=None,
            min_alt_freq=min_alt_freq,
            min_depth=min_depth,
            min_mapq=min_mapq,
            min_baseq=min_baseq,
            min_allele_count=min_allele_count,
        )
        first_row = rows[0] if rows else {}
        complete_rows.extend(rows)
        depth = int(first_row.get("depth") or 0)
        call_status = first_row.get("call_status") or "NO_CALL"
        if call_status == "NO_CALL":
            no_call_positions.append(
                (pos, depth, _site_filter(depth, min_depth))
            )
            continue

        callable_positions += 1
        if call_status == "MIXED_SIGNAL":
            mixed_positions.add(pos)

        significant_rows = [
            row
            for row in rows
            if row.get("variant_type") not in {"", "REF", "NO_CALL"}
            and row.get("allele_status") == "OBSERVED"
            and int(row.get("alt_depth") or 0) >= min_allele_count
            and float(row.get("allele_freq") or 0) >= min_alt_freq
        ]
        if not significant_rows:
            reference_positions += 1
            continue

        variant_positions.add(pos)
        variant_rows.extend(significant_rows)
        for row in significant_rows:
            variant_type = row.get("variant_type") or "OTHER"
            variant_type_counts[variant_type] = (
                variant_type_counts.get(variant_type, 0) + 1
            )

    no_call_regions = _compress_no_call_regions(
        gene, chrom, no_call_positions, min_depth
    )
    protein_haplotype_rows = build_protein_haplotypes(
        gene=gene,
        gene_record=gene_record,
        variant_rows=variant_rows,
        bam_path=bam_path,
        min_mapq=min_mapq,
        min_alt_freq=min_alt_freq,
        min_depth=min_depth,
        min_baseq=min_baseq,
        min_allele_count=min_allele_count,
        logger=logger,
    )
    _apply_protein_haplotype_annotations(
        variant_rows,
        protein_haplotype_rows,
    )
    gene_length = end - start + 1
    summary = {
        "scan_type": "whole_gene",
        "gene": gene,
        "gene_name": gene_record.get("description") or gene,
        "chrom": chrom,
        "start": start,
        "end": end,
        "gene_length": gene_length,
        "callable_positions": callable_positions,
        "callable_percent": (
            round(callable_positions / gene_length * 100, 2)
            if gene_length else 0
        ),
        "no_call_positions": len(no_call_positions),
        "no_call_regions": len(no_call_regions),
        "reference_only_positions": reference_positions,
        "variant_sites": len(variant_positions),
        "variant_rows": len(variant_rows),
        "complete_table_rows": len(complete_rows),
        "mixed_signal_sites": len(mixed_positions),
        "variant_type_counts": variant_type_counts,
        "protein_haplotype_clusters": len({
            row.get("cluster_id") for row in protein_haplotype_rows
            if row.get("cluster_id")
        }),
        "protein_haplotype_rows": len(protein_haplotype_rows),
        "frame_restored_haplotypes": sum(
            row.get("frame_status") == "FRAME_RESTORED"
            for row in protein_haplotype_rows
        ),
        "persistent_frameshift_haplotypes": sum(
            row.get("frame_status") == "PERSISTENT_FRAMESHIFT"
            for row in protein_haplotype_rows
        ),
        "thresholds": {
            "min_depth": min_depth,
            "min_mapq": min_mapq,
            "min_baseq": min_baseq,
            "min_allele_count": min_allele_count,
            "min_alt_freq": min_alt_freq,
        },
        "interpretation": (
            "MIXED_SIGNAL reports multiple supported alleles only; "
            "it does not diagnose mixed infection."
        ),
    }

    ensure_dir(out_csv.parent)
    write_protein_haplotypes(protein_haplotype_rows, protein_haplotype_csv)
    variant_rows.sort(key=_biological_position_sort_key)
    with open(out_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(variant_rows)
    with open(no_call_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=NO_CALL_REGION_FIELDS
        )
        writer.writeheader()
        writer.writerows(no_call_regions)
    with open(
        complete_table_csv, "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(complete_rows)
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    logger.info(
        f"Whole-gene scan wrote {len(variant_rows)} variant rows, "
        f"{len(no_call_regions)} NO_CALL regions: {out_csv}"
    )
    return {
        "table": out_csv,
        "summary": summary_path,
        "no_call_regions": no_call_csv,
        "complete_table": complete_table_csv,
        "protein_haplotypes": protein_haplotype_csv,
    }


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
                   min_depth: int,
                   min_baseq: int,
                   min_allele_count: int,
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
            min_depth=min_depth,
            min_mapq=min_mapq,
            min_baseq=min_baseq,
            min_allele_count=min_allele_count,
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
            min_baseq=min_baseq,
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
                min_depth=min_depth,
                min_mapq=min_mapq,
                min_baseq=min_baseq,
                min_allele_count=min_allele_count,
            ))

    ensure_dir(out_csv.parent)
    rows.sort(key=_biological_position_sort_key)
    with open(out_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    logger.info(f"Wrote {len(rows)} queried site rows: {out_csv}")
    return out_csv
