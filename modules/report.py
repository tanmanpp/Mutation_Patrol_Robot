# -*- coding: utf-8 -*-

import csv
import json
from pathlib import Path

from modules.gene_database import load_gene_database
from modules.utils import ensure_dir, make_json_safe, relative_path, resolve_project_path


SUMMARY_FIELDNAMES = [
    "Gene name", "Gene symbol", "Chromosome", "Start", "End",
    "Sequence status", "Variant count", "AA change count", "Synonymous count"
]


def _count_rows(csv_path: Path):
    if not csv_path or not csv_path.exists():
        return 0
    with open(csv_path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return sum(1 for _ in reader)


def _read_mutation_rows(csv_path: Path):
    if not csv_path or not csv_path.exists():
        return []
    with open(csv_path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_sample_summary(out_dir: Path, mutation_table: Path | None, extra: dict):
    gene_db_path = (extra or {}).get("gene_db")
    mutation_rows = _read_mutation_rows(mutation_table) if mutation_table else []
    rows_by_gene = {}
    for row in mutation_rows:
        gene_name = row.get("gene", "")
        rows_by_gene.setdefault(gene_name, []).append(row)

    summary_rows = []
    if gene_db_path:
        gene_db = load_gene_database(resolve_project_path(gene_db_path))
        for record in gene_db.get("genes", []):
            gene_symbol = record["symbol"]
            gene_rows = rows_by_gene.get(gene_symbol, [])
            aa_rows = [
                row for row in gene_rows
                if row.get("aa_change") and row.get("aa_change") != "NA"
            ]
            synonymous_rows = [
                row for row in aa_rows
                if "synonymous" in row.get("aa_change", "")
            ]
            summary_rows.append({
                "Gene name": record.get("description") or gene_symbol,
                "Gene symbol": gene_symbol,
                "Chromosome": record.get("resolved_chrom") or record.get("chrom", ""),
                "Start": record.get("gene_start", ""),
                "End": record.get("gene_end", ""),
                "Sequence status": record.get("sequence_status", ""),
                "Variant count": len(gene_rows),
                "AA change count": len(aa_rows),
                "Synonymous count": len(synonymous_rows),
            })
    else:
        for gene_symbol, gene_rows in sorted(rows_by_gene.items()):
            aa_rows = [
                row for row in gene_rows
                if row.get("aa_change") and row.get("aa_change") != "NA"
            ]
            synonymous_rows = [
                row for row in aa_rows
                if "synonymous" in row.get("aa_change", "")
            ]
            summary_rows.append({
                "Gene name": gene_symbol,
                "Gene symbol": gene_symbol,
                "Chromosome": gene_rows[0].get("chrom", "") if gene_rows else "",
                "Start": "",
                "End": "",
                "Sequence status": "",
                "Variant count": len(gene_rows),
                "AA change count": len(aa_rows),
                "Synonymous count": len(synonymous_rows),
            })

    summary_csv = out_dir / "sample_summary.csv"
    with open(summary_csv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDNAMES)
        writer.writeheader()
        writer.writerows(summary_rows)
    return summary_csv


def write_tables(out_dir: Path,
                 aa_table_path,
                 extra: dict,
                 logger):
    ensure_dir(out_dir)

    table_path = Path(aa_table_path) if aa_table_path else None
    sample_summary_path = _write_sample_summary(out_dir, table_path, extra or {})
    summary = {
        "mutation_table": relative_path(table_path) if table_path else None,
        "sample_summary_table": relative_path(sample_summary_path),
        "mutation_rows": _count_rows(table_path) if table_path else 0,
        "extra": make_json_safe(extra or {}),
    }

    summary_path = out_dir / "summary.json"
    with open(summary_path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)

    logger.info(f"Wrote report summary: {summary_path}")
    logger.info(f"Wrote sample summary table: {sample_summary_path}")
    return summary_path
