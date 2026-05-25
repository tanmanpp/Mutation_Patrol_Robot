# -*- coding: utf-8 -*-

import csv
import json
import shutil
from pathlib import Path

from fastapi import UploadFile

from modules.bam_allele_freq import mutation_candidates_from_bam
from modules.gene_database import build_gene_database_from_annotations, load_gene_database
from modules.mapping import map_reads_to_ref
from modules.report import write_tables
from modules.roi_extract import extract_gene_bams
from modules.site_query import query_bam_sites
from modules.utils import ensure_dir, make_json_safe, relative_path, resolve_project_path, setup_logger, write_run_metadata
from modules.variant_call import call_variants_per_gene
from modules.annotate_mut import variants_to_aa_table

from .config import DATABASES_DIR, LEGACY_DATABASES_DIR, RESULTS_DIR, UPLOADS_DIR


def safe_name(value: str):
    cleaned = "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in value.strip())
    return cleaned or "run"


async def save_upload(upload: UploadFile, dest: Path):
    ensure_dir(dest.parent)
    with open(dest, "wb") as handle:
        while True:
            chunk = await upload.read(1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)
    await upload.close()
    return dest


def read_csv_records(path: Path, limit: int | None = None):
    if not path.exists():
        return []
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = []
        for row in reader:
            rows.append(row)
            if limit and len(rows) >= limit:
                break
        return rows


def list_databases():
    databases = []
    seen = set()
    db_paths = list(sorted(DATABASES_DIR.glob("*/gene_database.json")))
    if LEGACY_DATABASES_DIR.exists():
        db_paths.extend(sorted(LEGACY_DATABASES_DIR.glob("*/gene_database.json")))

    for db_json in db_paths:
        if db_json.parent.name in seen:
            continue
        seen.add(db_json.parent.name)
        payload = load_gene_database(db_json)
        databases.append({
            "name": db_json.parent.name,
            "path": relative_path(db_json),
            "gene_count": payload.get("gene_count", len(payload.get("genes", []))),
            "warning_count": payload.get("warning_count", 0),
            "genes": [
                {
                    "symbol": gene.get("symbol", ""),
                    "description": gene.get("description", ""),
                    "sequence_status": gene.get("sequence_status", ""),
                }
                for gene in payload.get("genes", [])
            ],
        })
    return databases


def get_database_path(name: str):
    db_name = safe_name(name)
    candidates = [
        DATABASES_DIR / db_name / "gene_database.json",
        LEGACY_DATABASES_DIR / db_name / "gene_database.json",
    ]
    for db_path in candidates:
        if db_path.exists():
            return db_path
    raise FileNotFoundError(f"Database not found: {name}")


def get_database_genes(name: str):
    payload = load_gene_database(get_database_path(name))
    return [
        {
            "symbol": gene.get("symbol", ""),
            "description": gene.get("description", ""),
            "chrom": gene.get("resolved_chrom") or gene.get("chrom", ""),
            "gene_start": gene.get("gene_start", ""),
            "gene_end": gene.get("gene_end", ""),
            "sequence_status": gene.get("sequence_status", ""),
        }
        for gene in payload.get("genes", [])
    ]


def build_database(name: str, ref_fasta: Path, annotations: list[Path], genes: list[str] | None, force: bool):
    db_name = safe_name(name)
    out_dir = DATABASES_DIR / db_name
    ensure_dir(out_dir)
    logger = setup_logger(out_dir / "build_gene_db.log")
    outputs = build_gene_database_from_annotations(
        annotation_paths=annotations,
        ref_fasta=ref_fasta,
        out_dir=out_dir,
        genes=genes,
        force=force,
    )
    for label, path in outputs.items():
        logger.info(f"{label}: {path}")

    db_payload = load_gene_database(outputs["json"])
    return {
        "name": db_name,
        "outputs": {key: relative_path(value) for key, value in outputs.items()},
        "gene_count": db_payload.get("gene_count", 0),
        "warning_count": db_payload.get("warning_count", 0),
        "genes": get_database_genes(db_name),
        "summary_rows": read_csv_records(outputs["csv"]),
    }


def analyze_sample(
    run_name: str,
    db_name: str,
    fastq_path: Path | None,
    bam_path: Path | None,
    ref_fasta: Path | None,
    min_depth: int,
    threads: int,
    min_mapq: int,
    min_alt_count: int,
    min_alt_freq: float,
    candidate_source: str,
    force: bool,
):
    run_id = safe_name(run_name)
    out_dir = RESULTS_DIR / run_id
    ensure_dir(out_dir)
    logger = setup_logger(out_dir / "analyze_sample.log")

    db_path = get_database_path(db_name)
    gene_db = load_gene_database(db_path)
    ref_path = Path(ref_fasta) if ref_fasta else resolve_project_path(gene_db["reference_fasta"])

    work = out_dir / "work"
    bam_dir = work / "bam"
    roi_dir = work / "roi_bam"
    vcf_dir = work / "vcf"
    tables_dir = out_dir / "tables"
    for directory in [work, bam_dir, roi_dir, vcf_dir, tables_dir]:
        ensure_dir(directory)

    bed_path = db_path.parent / "gene_regions.bed"
    if bam_path:
        source_bam = bam_path
        logger.info(f"Using uploaded BAM: {source_bam}")
    elif fastq_path:
        source_bam = map_reads_to_ref(
            raw_dir=None,
            out_dir=bam_dir,
            ref_fasta=ref_path,
            threads=threads,
            min_mapq=min_mapq,
            dry_run=False,
            force=force,
            logger=logger,
            read_files=[fastq_path],
        )
    else:
        raise ValueError("Need FASTQ or BAM for sample analysis.")

    gene_bams = extract_gene_bams(
        bam_path=source_bam,
        bed_path=bed_path,
        out_dir=roi_dir,
        threads=threads,
        dry_run=False,
        force=force,
        logger=logger,
    )

    if candidate_source == "bam":
        mutation_table = mutation_candidates_from_bam(
            gene_bams=gene_bams,
            gene_db_path=db_path,
            ref_fasta=ref_path,
            out_dir=tables_dir,
            min_depth=min_depth,
            min_alt_count=min_alt_count,
            min_alt_freq=min_alt_freq,
            min_mapq=min_mapq,
            dry_run=False,
            force=force,
            logger=logger,
        )
    else:
        vcf_paths = call_variants_per_gene(
            gene_bams=gene_bams,
            ref_fasta=ref_path,
            out_dir=vcf_dir,
            caller="bcftools",
            threads=threads,
            min_depth=min_depth,
            dry_run=False,
            force=force,
            logger=logger,
        )
        mutation_table = variants_to_aa_table(
            vcf_paths=vcf_paths,
            annotation_path=None,
            bed_path=bed_path,
            out_dir=tables_dir,
            dry_run=False,
            force=force,
            logger=logger,
            gene_db_path=db_path,
        )

    write_tables(
        out_dir=tables_dir,
        aa_table_path=mutation_table,
        extra={
            "gene_db": relative_path(db_path),
            "gene_count": gene_db.get("gene_count", 0),
            "candidate_source": candidate_source,
            "threads": threads,
        },
        logger=logger,
    )

    return get_result(run_id)


def query_sites(
    analysis_run: str,
    db_name: str,
    gene: str,
    query_type: str,
    query_value: str,
    alt: str | None,
    min_alt_freq: float,
    min_mapq: int,
    force: bool,
):
    analysis_id = safe_name(analysis_run)
    query_id = safe_name(gene)
    analysis_dir = RESULTS_DIR / analysis_id
    bam_path = analysis_dir / "work" / "bam" / "merged.sorted.bam"
    if not bam_path.exists():
        raise FileNotFoundError(f"Analysis BAM not found: {bam_path}")

    out_dir = analysis_dir / "site_query" / query_id
    ensure_dir(out_dir)
    logger = setup_logger(out_dir / "query_sites.log")
    out_csv = out_dir / "site_query.csv"
    kwargs = {"genomic_pos": None, "cds_pos": None, "aa_pos": None}
    values = parse_query_values(query_value)
    kwargs[query_type] = values
    query_bam_sites(
        gene_db_path=get_database_path(db_name),
        bam_path=bam_path,
        out_csv=out_csv,
        gene=gene,
        alt=alt,
        ref_fasta=None,
        min_mapq=min_mapq,
        min_alt_freq=min_alt_freq,
        force=force,
        logger=logger,
        **kwargs,
    )
    return {
        "run_id": analysis_id,
        "query_id": query_id,
        "gene": gene,
        "site_query_table": relative_path(out_csv),
        "rows": read_csv_records(out_csv),
    }


def parse_query_values(value: str):
    values = []
    for item in value.split(","):
        item = item.strip()
        if not item:
            continue
        try:
            parsed = int(item)
        except ValueError as exc:
            raise ValueError(f"Query positions must be comma-separated integers: {value}") from exc
        if parsed < 1:
            raise ValueError("Query positions must be 1-based and >= 1.")
        values.append(parsed)
    if not values:
        raise ValueError("Provide at least one query position.")
    return values


def clear_uploaded_samples():
    samples_dir = UPLOADS_DIR / "samples"
    before = sample_upload_stats()
    if samples_dir.exists():
        shutil.rmtree(samples_dir)
    ensure_dir(samples_dir)
    return {
        "cleared": True,
        "path": relative_path(samples_dir),
        "deleted_files": before["file_count"],
        "deleted_dirs": before["dir_count"],
        "deleted_bytes": before["total_bytes"],
        "stats": sample_upload_stats(),
    }


def sample_upload_stats():
    samples_dir = UPLOADS_DIR / "samples"
    total_bytes = 0
    file_count = 0
    dir_count = 0
    if samples_dir.exists():
        for path in samples_dir.rglob("*"):
            if path.is_file():
                file_count += 1
                total_bytes += path.stat().st_size
            elif path.is_dir():
                dir_count += 1
    return {
        "path": relative_path(samples_dir),
        "total_bytes": total_bytes,
        "file_count": file_count,
        "dir_count": dir_count,
    }


def get_result(run_id: str):
    result_dir = RESULTS_DIR / safe_name(run_id)
    tables = result_dir / "tables"
    summary_path = tables / "summary.json"
    summary = {}
    if summary_path.exists():
        summary = make_json_safe(json.loads(summary_path.read_text(encoding="utf-8")))
    return {
        "run_id": result_dir.name,
        "summary": summary,
        "sample_summary": read_csv_records(tables / "sample_summary.csv"),
        "mutation_candidates": read_csv_records(tables / "mutation_candidates.csv", limit=500),
        "site_query": read_latest_site_query(result_dir),
        "site_queries": read_site_queries(result_dir),
    }


def read_latest_site_query(result_dir: Path):
    site_root = result_dir / "site_query"
    if not site_root.exists():
        return []
    query_dirs = [path for path in site_root.iterdir() if path.is_dir()]
    if not query_dirs:
        return []
    latest = max(query_dirs, key=lambda path: path.stat().st_mtime)
    return read_csv_records(latest / "site_query.csv", limit=500)


def read_site_queries(result_dir: Path):
    site_root = result_dir / "site_query"
    if not site_root.exists():
        return []
    items = []
    query_dirs = sorted(
        [path for path in site_root.iterdir() if path.is_dir()],
        key=lambda path: path.stat().st_mtime,
    )
    for query_dir in query_dirs:
        csv_path = query_dir / "site_query.csv"
        if not csv_path.exists():
            continue
        items.append({
            "query_id": query_dir.name,
            "path": relative_path(csv_path),
            "updated_at": query_dir.stat().st_mtime,
            "rows": read_csv_records(csv_path, limit=500),
        })
    return items


def list_results():
    runs = []
    for result_dir in sorted(RESULTS_DIR.iterdir()) if RESULTS_DIR.exists() else []:
        if result_dir.is_dir():
            runs.append(get_result(result_dir.name))
    return runs


def copy_uploaded_bam_to_run(upload_path: Path, run_name: str):
    run_dir = UPLOADS_DIR / safe_name(run_name)
    ensure_dir(run_dir)
    dest = run_dir / upload_path.name
    if upload_path.resolve() != dest.resolve():
        shutil.copy2(upload_path, dest)
    return dest
