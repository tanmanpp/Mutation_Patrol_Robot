# -*- coding: utf-8 -*-

import csv
import json
from pathlib import Path

from modules.utils import ensure_dir, relative_path, run_cmd


CSV_FIELDNAMES = [
    "symbol", "chrom", "strand", "type", "order", "start", "end",
    "length", "cds_offset"
]


def _normalize_strand(orientation):
    if orientation is None:
        return "."
    if isinstance(orientation, int):
        return "+" if orientation >= 0 else "-"
    if isinstance(orientation, str):
        value = orientation.strip().lower()
        if value in ("plus", "+", "forward", "fwd", "1", "positive"):
            return "+"
        if value in ("minus", "-", "reverse", "rev", "-1", "negative"):
            return "-"
    return "."


def _parse_ncbi_product_report_line(line: str):
    data = json.loads(line)
    transcripts = data.get("transcripts") or []
    if not transcripts:
        return None

    transcript = transcripts[0]
    locations = transcript.get("genomicLocations") or []
    if not locations:
        return None

    location = locations[0]
    exons = location.get("exons") or []
    if not exons:
        return None

    symbol = data.get("symbol", "N/A")
    chrom = location.get("genomicAccessionVersion")
    genomic_range = location.get("genomicRange") or {}
    strand = _normalize_strand(genomic_range.get("orientation"))

    sorted_exons = sorted(exons, key=lambda exon: int(exon.get("order", 0)))
    exon_min = min(int(exon["begin"]) for exon in sorted_exons)
    exon_max = max(int(exon["end"]) for exon in sorted_exons)

    gr_begin = genomic_range.get("begin")
    gr_end = genomic_range.get("end")
    gene_start = exon_min if gr_begin is None else min(int(gr_begin), exon_min)
    gene_end = exon_max if gr_end is None else max(int(gr_end), exon_max)

    rows = [{
        "symbol": symbol,
        "chrom": chrom,
        "strand": strand,
        "type": "gene",
        "order": 0,
        "start": gene_start,
        "end": gene_end,
        "length": gene_end - gene_start + 1,
        "cds_offset": 0,
    }]

    offset = 0
    for exon in sorted_exons:
        start = int(exon["begin"])
        end = int(exon["end"])
        length = end - start + 1
        rows.append({
            "symbol": symbol,
            "chrom": chrom,
            "strand": strand,
            "type": "exon",
            "order": int(exon.get("order", 0)),
            "start": start,
            "end": end,
            "length": length,
            "cds_offset": offset,
        })
        offset += length

    bed = {
        "chrom": chrom,
        "start0": gene_start - 1,
        "end0": gene_end,
        "name": symbol,
        "score": 0,
        "strand": strand,
    }
    return rows, bed


def _convert_jsonl_to_tables(annotation_path: Path,
                             bed_out: Path,
                             csv_out: Path,
                             genes: list[str] | None,
                             force: bool,
                             logger):
    if bed_out.exists() and csv_out.exists() and not force:
        logger.info(f"Annotation tables exist, skip: {bed_out}, {csv_out}")
        return

    selected = set(genes or [])
    seen_bed = set()
    parsed_count = 0

    with open(annotation_path, "r", encoding="utf-8") as infile, \
            open(csv_out, "w", newline="", encoding="utf-8") as out_csv, \
            open(bed_out, "w", newline="", encoding="utf-8") as out_bed:
        writer = csv.DictWriter(out_csv, fieldnames=CSV_FIELDNAMES)
        writer.writeheader()

        for line in infile:
            if not line.strip():
                continue
            parsed = _parse_ncbi_product_report_line(line)
            if not parsed:
                continue
            rows, bed = parsed
            if selected and bed["name"] not in selected:
                continue

            for row in rows:
                writer.writerow(row)

            key = (bed["chrom"], bed["start0"], bed["end0"], bed["name"], bed["strand"])
            if key not in seen_bed:
                seen_bed.add(key)
                out_bed.write(
                    f'{bed["chrom"]}\t{bed["start0"]}\t{bed["end0"]}\t'
                    f'{bed["name"]}\t{bed["score"]}\t{bed["strand"]}\n'
                )
            parsed_count += 1

    if parsed_count == 0:
        raise ValueError(f"No matching gene records parsed from {annotation_path}")

    logger.info(f"Parsed {parsed_count} gene records into {bed_out}")


def _copy_filtered_bed(source_bed: Path,
                       dest_bed: Path,
                       genes: list[str] | None,
                       force: bool,
                       logger):
    if not genes:
        return source_bed

    if dest_bed.exists() and not force:
        logger.info(f"Filtered BED exists, skip: {dest_bed}")
        return dest_bed

    selected = set(genes)
    kept = 0
    with open(source_bed, "r", encoding="utf-8") as src, \
            open(dest_bed, "w", encoding="utf-8", newline="") as dst:
        for line in src:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                raise ValueError("BED must have at least 4 columns: chrom start end gene")
            if parts[3] in selected:
                dst.write(line)
                kept += 1

    if kept == 0:
        raise ValueError(f"No requested genes found in {source_bed}: {', '.join(genes)}")
    logger.info(f"Filtered BED to {kept} regions: {dest_bed}")
    return dest_bed


def prepare_annotation(out_dir: Path,
                       ref_id: str | None,
                       annotation: str | None,
                       genes: list[str] | None,
                       gene_bed: str | None,
                       dry_run: bool,
                       force: bool,
                       logger):
    """
    Return dict of paths:
      {"bed": ".../gene_regions.bed", "annotation_jsonl": "...", "annotation_csv": "..."}

    Supported in this first version:
      - Existing BED.
      - Existing NCBI product_report JSONL.
      - A ref_id only when a project-local scripts/json_spider.py is provided.
    """
    ensure_dir(out_dir)

    if gene_bed:
        source_bed = Path(gene_bed).resolve()
        if not source_bed.exists():
            raise FileNotFoundError(f"gene_bed not found: {source_bed}")
        bed_path = _copy_filtered_bed(
            source_bed, out_dir / "gene_regions.filtered.bed", genes, force, logger
        )
        return {"bed": relative_path(bed_path), "annotation_jsonl": None, "annotation_csv": None}

    if annotation:
        annotation_path = Path(annotation).resolve()
        if not annotation_path.exists():
            raise FileNotFoundError(f"annotation not found: {annotation_path}")
        if annotation_path.suffix.lower() != ".jsonl":
            raise NotImplementedError("Initial version supports NCBI product_report .jsonl annotations only.")

        bed_out = out_dir / "gene_regions.bed"
        csv_out = out_dir / "gene_regions.csv"
        _convert_jsonl_to_tables(annotation_path, bed_out, csv_out, genes, force, logger)
        return {
            "bed": relative_path(bed_out),
            "annotation_jsonl": relative_path(annotation_path),
            "annotation_csv": relative_path(csv_out),
        }

    if ref_id:
        project_root = Path(__file__).resolve().parents[1]
        spider = project_root / "scripts" / "json_spider.py"
        if not spider.exists():
            legacy_spider = project_root / "test_script" / "json_spider.py"
            if legacy_spider.exists():
                spider = legacy_spider

        if not spider.exists():
            raise FileNotFoundError("No scripts/json_spider.py or test_script/json_spider.py found.")

        annotation_path = out_dir / f"{ref_id}.jsonl"
        cmd = ["python", str(spider), "--ref_id", ref_id, "--out", str(annotation_path)]
        run_cmd(cmd, logger=logger, dry_run=dry_run)
        if dry_run:
            logger.warning("Dry run: annotation JSONL was not created, so BED conversion is skipped.")
            return {
                "bed": relative_path(out_dir / "gene_regions.bed"),
                "annotation_jsonl": relative_path(annotation_path),
                "annotation_csv": relative_path(out_dir / "gene_regions.csv"),
            }

        bed_out = out_dir / "gene_regions.bed"
        csv_out = out_dir / "gene_regions.csv"
        _convert_jsonl_to_tables(annotation_path, bed_out, csv_out, genes, force, logger)
        return {
            "bed": relative_path(bed_out),
            "annotation_jsonl": relative_path(annotation_path),
            "annotation_csv": relative_path(csv_out),
        }

    raise ValueError("Need one of: --gene_bed, --annotation, or --ref_id")
