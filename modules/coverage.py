# -*- coding: utf-8 -*-

import csv
import hashlib
import json
import re
import statistics
from pathlib import Path

from modules.bam_utils import ensure_bam_index
from modules.gene_database import load_gene_database
from modules.utils import ensure_dir, run_cmd


GENE_COVERAGE_FIELDS = [
    "Gene name", "Gene symbol", "chrom", "start", "end", "length",
    "covered_bases", "callable_bases", "no_call_bases", "mean_depth",
    "median_depth", "min_depth_observed", "max_depth_observed",
    "pct_covered", "pct_callable", "minimum_call_depth", "qc_status",
]

REGION_COVERAGE_FIELDS = [
    "Gene name", "Gene symbol", "region_id", "region_type", "chrom",
    "start", "end", "length", "covered_bases", "callable_bases",
    "no_call_bases", "mean_depth", "median_depth", "min_depth_observed",
    "max_depth_observed", "pct_covered", "pct_callable",
    "minimum_call_depth", "qc_status",
]

POSITION_COVERAGE_FIELDS = [
    "gene", "chrom", "pos", "relative_pos", "depth", "call_status",
]


def _safe_file_stem(value: str):
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._")
    return cleaned or "gene"


def _coverage_stats(depths: list[int], minimum_call_depth: int):
    length = len(depths)
    covered_bases = sum(depth > 0 for depth in depths)
    callable_bases = sum(depth >= minimum_call_depth for depth in depths)
    no_call_bases = length - callable_bases
    pct_covered = covered_bases / length if length else 0
    pct_callable = callable_bases / length if length else 0
    if length == 0:
        qc_status = "NO_CALL"
    elif covered_bases == 0:
        qc_status = "NO_COVERAGE"
    elif callable_bases == 0:
        qc_status = "LOW_DEPTH"
    elif callable_bases < length:
        qc_status = "PARTIAL"
    else:
        qc_status = "PASS"
    return {
        "length": length,
        "covered_bases": covered_bases,
        "callable_bases": callable_bases,
        "no_call_bases": no_call_bases,
        "mean_depth": f"{statistics.fmean(depths):.2f}" if depths else "0.00",
        "median_depth": f"{statistics.median(depths):.2f}" if depths else "0.00",
        "min_depth_observed": min(depths) if depths else 0,
        "max_depth_observed": max(depths) if depths else 0,
        "pct_covered": f"{pct_covered:.6f}",
        "pct_callable": f"{pct_callable:.6f}",
        "minimum_call_depth": minimum_call_depth,
        "qc_status": qc_status,
    }


def _run_samtools_depth(bam_path: Path,
                        chrom: str,
                        start: int,
                        end: int,
                        min_mapq: int,
                        min_baseq: int,
                        logger):
    region = f"{chrom}:{start}-{end}"
    completed = run_cmd(
        [
            "samtools", "depth", "-aa",
            "-q", str(min_baseq),
            "-Q", str(min_mapq),
            "-r", region,
            str(bam_path),
        ],
        logger=logger,
        capture_output=True,
    )
    depths_by_position = {}
    for line in completed.stdout.splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        depths_by_position[int(parts[1])] = int(parts[2])
    return [
        depths_by_position.get(position, 0)
        for position in range(start, end + 1)
    ]


def _coverage_regions(gene_record: dict):
    amplicons = gene_record.get("amplicons") or []
    if amplicons:
        return [
            {
                "region_id": item.get("name") or f"amplicon_{index}",
                "region_type": "amplicon",
                "start": int(item["start"]),
                "end": int(item["end"]),
            }
            for index, item in enumerate(amplicons, start=1)
        ]
    return [
        {
            "region_id": f"exon_{exon.get('order') or index}",
            "region_type": "exon",
            "start": int(exon["start"]),
            "end": int(exon["end"]),
        }
        for index, exon in enumerate(gene_record.get("exons") or [], start=1)
    ]


def generate_coverage_tables(gene_db_path: Path,
                             bam_path: Path,
                             out_dir: Path,
                             minimum_call_depth: int,
                             min_mapq: int,
                             min_baseq: int,
                             force: bool,
                             logger):
    summary_csv = out_dir / "gene_coverage.csv"
    regions_csv = out_dir / "region_coverage.csv"
    index_json = out_dir / "coverage_index.json"
    if summary_csv.exists() and regions_csv.exists() and index_json.exists() and not force:
        logger.info(f"Coverage tables exist, skip: {out_dir}")
        return {
            "summary": summary_csv,
            "regions": regions_csv,
            "index": index_json,
        }

    ensure_dir(out_dir)
    positions_dir = out_dir / "positions"
    ensure_dir(positions_dir)
    ensure_bam_index(bam_path, logger)
    gene_db = load_gene_database(gene_db_path)

    summary_rows = []
    region_rows = []
    coverage_index = {
        "gene_db": str(gene_db_path),
        "minimum_call_depth": minimum_call_depth,
        "min_mapq": min_mapq,
        "min_baseq": min_baseq,
        "genes": {},
    }

    for gene_record in gene_db.get("genes") or []:
        gene = gene_record["symbol"]
        gene_name = gene_record.get("description") or gene
        chrom = gene_record.get("resolved_chrom") or gene_record.get("chrom")
        start = int(gene_record["gene_start"])
        end = int(gene_record["gene_end"])
        if not chrom or start < 1 or end < start:
            logger.warning(f"Coverage skipped for invalid gene interval: {gene}")
            continue

        depths = _run_samtools_depth(
            bam_path=bam_path,
            chrom=chrom,
            start=start,
            end=end,
            min_mapq=min_mapq,
            min_baseq=min_baseq,
            logger=logger,
        )
        stats = _coverage_stats(depths, minimum_call_depth)
        summary_rows.append({
            "Gene name": gene_name,
            "Gene symbol": gene,
            "chrom": chrom,
            "start": start,
            "end": end,
            **stats,
        })

        gene_hash = hashlib.sha1(gene.encode("utf-8")).hexdigest()[:10]
        file_name = f"{_safe_file_stem(gene)}_{gene_hash}.csv"
        position_path = positions_dir / file_name
        with open(position_path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=POSITION_COVERAGE_FIELDS)
            writer.writeheader()
            for offset, depth in enumerate(depths):
                writer.writerow({
                    "gene": gene,
                    "chrom": chrom,
                    "pos": start + offset,
                    "relative_pos": offset + 1,
                    "depth": depth,
                    "call_status": (
                        "CALLABLE" if depth >= minimum_call_depth else "NO_CALL"
                    ),
                })
        coverage_index["genes"][gene] = {
            "file": file_name,
            "chrom": chrom,
            "start": start,
            "end": end,
        }

        for region in _coverage_regions(gene_record):
            region_start = max(start, int(region["start"]))
            region_end = min(end, int(region["end"]))
            if region_start > region_end:
                continue
            left = region_start - start
            right = region_end - start + 1
            region_stats = _coverage_stats(
                depths[left:right], minimum_call_depth
            )
            region_rows.append({
                "Gene name": gene_name,
                "Gene symbol": gene,
                "region_id": region["region_id"],
                "region_type": region["region_type"],
                "chrom": chrom,
                "start": region_start,
                "end": region_end,
                **region_stats,
            })

    with open(summary_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=GENE_COVERAGE_FIELDS)
        writer.writeheader()
        writer.writerows(summary_rows)

    with open(regions_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=REGION_COVERAGE_FIELDS)
        writer.writeheader()
        writer.writerows(region_rows)

    index_json.write_text(
        json.dumps(coverage_index, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    logger.info(
        f"Wrote coverage for {len(summary_rows)} genes and "
        f"{len(region_rows)} regions: {out_dir}"
    )
    return {
        "summary": summary_csv,
        "regions": regions_csv,
        "index": index_json,
    }
