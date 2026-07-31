# -*- coding: utf-8 -*-

import csv
from pathlib import Path

from modules.annotate_mut import FIELDNAMES, annotate_candidate_row
from modules.bam_utils import ensure_bam_index
from modules.gene_database import load_gene_database
from modules.utils import ensure_dir, run_cmd


BASES = ("A", "C", "G", "T")


def _format_allele_spectrum(counts: dict[str, int], depth: int):
    observed = sorted(
        ((allele, count) for allele, count in counts.items() if count > 0),
        key=lambda item: (-item[1], item[0]),
    )
    return "; ".join(
        f"{allele}:{count} ({count / depth:.1%})"
        for allele, count in observed
    ) if depth else ""


def _candidate_call_status(counts: dict[str, int],
                           depth: int,
                           min_allele_count: int,
                           min_allele_freq: float):
    qualifying = [
        allele
        for allele, count in counts.items()
        if count > 0
        and count >= min_allele_count
        and depth
        and count / depth >= min_allele_freq
    ]
    return "MIXED_SIGNAL" if len(qualifying) >= 2 else "VARIANT"


def _count_pileup_bases(read_bases: str, ref_base: str):
    return _count_pileup_events(read_bases, ref_base)["bases"]


def _count_pileup_events(read_bases: str, ref_base: str):
    counts = {base: 0 for base in BASES}
    strand_counts = {
        base: {"forward": 0, "reverse": 0}
        for base in BASES
    }
    insertions = {}
    deletions = {}
    i = 0
    while i < len(read_bases):
        char = read_bases[i]

        if char == "^":
            i += 2
            continue
        if char == "$":
            i += 1
            continue
        if char in "+-":
            event_type = char
            i += 1
            digits = []
            while i < len(read_bases) and read_bases[i].isdigit():
                digits.append(read_bases[i])
                i += 1
            skip = int("".join(digits)) if digits else 0
            sequence = read_bases[i:i + skip].upper()
            if sequence:
                if event_type == "+":
                    insertions[sequence] = insertions.get(sequence, 0) + 1
                else:
                    deletions[sequence] = deletions.get(sequence, 0) + 1
            i += skip
            continue
        if char in ".,":  # reference base on forward/reverse strand
            ref = ref_base.upper()
            if ref in counts:
                counts[ref] += 1
                strand = "forward" if char == "." else "reverse"
                strand_counts[ref][strand] += 1
            i += 1
            continue

        base = char.upper()
        if base in counts:
            counts[base] += 1
            strand = "forward" if char.isupper() else "reverse"
            strand_counts[base][strand] += 1
        i += 1

    return {
        "bases": counts,
        "strand_counts": strand_counts,
        "insertions": insertions,
        "deletions": deletions,
    }


def _run_mpileup(bam_path: Path,
                 ref_fasta: Path,
                 chrom: str,
                 start: int,
                 end: int,
                 min_mapq: int,
                 logger,
                 min_baseq: int = 20):
    region = f"{chrom}:{start}-{end}"
    try:
        completed = run_cmd(
            [
                "samtools", "mpileup", "-aa", "-q", str(min_mapq),
                "-Q", str(min_baseq),
                "-d", "1000000",
                "-f", str(ref_fasta), "-r", region, str(bam_path),
            ],
            logger=logger,
            capture_output=True,
        )
    except Exception as exc:
        raise RuntimeError(
            "samtools mpileup failed. Check that the BAM has an index, "
            "the queried contig exists in the BAM, and the BAM/reference use matching IDs."
        ) from exc
    return completed.stdout.splitlines()


def _candidate_rows_from_pileup(gene: str,
                                gene_record: dict,
                                pileup_lines: list[str],
                                gene_db: dict,
                                min_depth: int,
                                min_alt_count: int,
                                min_alt_freq: float,
                                min_mapq: int,
                                min_baseq: int):
    rows = []
    for line in pileup_lines:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 5:
            continue

        chrom, pos, ref_base, raw_depth, read_bases = parts[:5]
        events = _count_pileup_events(read_bases, ref_base)
        counts = events["bases"]
        strand_counts = events["strand_counts"]
        counted_depth = sum(counts.values())
        if counted_depth < min_depth:
            continue

        ref = ref_base.upper()
        ref_depth = counts.get(ref, 0)
        call_status = _candidate_call_status(
            counts, counted_depth, min_alt_count, min_alt_freq
        )
        spectrum = _format_allele_spectrum(counts, counted_depth)
        allele_count = sum(1 for count in counts.values() if count > 0)
        for alt in BASES:
            if alt == ref:
                continue
            alt_depth = counts[alt]
            if alt_depth < min_alt_count:
                continue
            allele_freq = alt_depth / counted_depth if counted_depth else 0
            if allele_freq < min_alt_freq:
                continue

            row = {
                "Gene name": gene_record.get("description") or gene,
                "Gene symbol": gene,
                "gene": gene,
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "qual": "",
                "filter": "PASS",
                "depth": counted_depth,
                "raw_depth": raw_depth,
                "ref_depth": ref_depth,
                "alt_depth": alt_depth,
                "forward_depth": strand_counts[alt]["forward"],
                "reverse_depth": strand_counts[alt]["reverse"],
                "allele_freq": f"{allele_freq:.6f}",
                "allele_count": allele_count,
                "allele_spectrum": spectrum,
                "call_status": call_status,
                "allele_status": "OBSERVED",
                "qc_flags": (
                    "MULTIPLE_ALLELES"
                    if call_status == "MIXED_SIGNAL" else "PASS"
                ),
                "min_depth": min_depth,
                "min_mapq": min_mapq,
                "min_baseq": min_baseq,
                "min_allele_count": min_alt_count,
                "variant_type": "SNV",
                "source": "bam_pileup",
            }
            rows.append(annotate_candidate_row(row, gene_db))

        for inserted_sequence, alt_depth in sorted(events["insertions"].items()):
            if alt_depth < min_alt_count:
                continue
            allele_freq = alt_depth / counted_depth if counted_depth else 0
            if allele_freq < min_alt_freq:
                continue

            insertion_call_counts = {
                ref: max(counted_depth - alt_depth, 0),
                f"+{inserted_sequence}": alt_depth,
            }
            insertion_call_status = _candidate_call_status(
                insertion_call_counts,
                counted_depth,
                min_alt_count,
                min_alt_freq,
            )
            row = {
                "Gene name": gene_record.get("description") or gene,
                "Gene symbol": gene,
                "gene": gene,
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": f"{ref}{inserted_sequence}",
                "qual": "",
                "filter": "PASS",
                "depth": counted_depth,
                "raw_depth": raw_depth,
                "ref_depth": ref_depth,
                "alt_depth": alt_depth,
                "forward_depth": "",
                "reverse_depth": "",
                "allele_freq": f"{allele_freq:.6f}",
                "allele_count": sum(
                    count > 0 for count in insertion_call_counts.values()
                ),
                "allele_spectrum": _format_allele_spectrum(
                    insertion_call_counts, counted_depth
                ),
                "call_status": insertion_call_status,
                "allele_status": "OBSERVED",
                "qc_flags": (
                    "MULTIPLE_ALLELES"
                    if insertion_call_status == "MIXED_SIGNAL" else "PASS"
                ),
                "min_depth": min_depth,
                "min_mapq": min_mapq,
                "min_baseq": min_baseq,
                "min_allele_count": min_alt_count,
                "variant_type": "INS",
                "source": "bam_pileup",
                "note": f"insertion_sequence={inserted_sequence}",
            }
            rows.append(annotate_candidate_row(row, gene_db))
    return rows


def mutation_candidates_from_bam(gene_bams: dict,
                                 gene_db_path: Path,
                                 ref_fasta: Path,
                                 out_dir: Path,
                                 min_depth: int,
                                 min_alt_count: int,
                                 min_alt_freq: float,
                                 min_mapq: int,
                                 dry_run: bool,
                                 force: bool,
                                 logger,
                                 min_baseq: int = 20):
    ensure_dir(out_dir)
    out_csv = out_dir / "mutation_candidates.csv"
    if out_csv.exists() and not force:
        logger.info(f"Mutation table exists, skip: {out_csv}")
        return out_csv

    gene_db = load_gene_database(gene_db_path)
    rows = []

    if dry_run:
        logger.info(f"Dry run: would call samtools mpileup and write {out_csv}")
    else:
        for gene, bam_path in gene_bams.items():
            gene_record = gene_db.get("genes_by_symbol", {}).get(gene)
            if not gene_record:
                logger.warning(f"Gene BAM has no database record, skip: {gene}")
                continue

            bam_path = Path(bam_path)
            if not bam_path.exists():
                raise FileNotFoundError(f"Gene BAM not found: {bam_path}")
            ensure_bam_index(bam_path, logger)

            chrom = gene_record.get("resolved_chrom") or gene_record["chrom"]
            start = int(gene_record["gene_start"])
            end = int(gene_record["gene_end"])
            pileup_lines = _run_mpileup(
                bam_path=bam_path,
                ref_fasta=ref_fasta,
                chrom=chrom,
                start=start,
                end=end,
                min_mapq=min_mapq,
                logger=logger,
                min_baseq=min_baseq,
            )
            rows.extend(_candidate_rows_from_pileup(
                gene=gene,
                gene_record=gene_record,
                pileup_lines=pileup_lines,
                gene_db=gene_db,
                min_depth=min_depth,
                min_alt_count=min_alt_count,
                min_alt_freq=min_alt_freq,
                min_mapq=min_mapq,
                min_baseq=min_baseq,
            ))

    with open(out_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    logger.info(f"Wrote {len(rows)} BAM allele rows: {out_csv}")
    return out_csv
