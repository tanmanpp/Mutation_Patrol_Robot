#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
from modules.utils import (
    ensure_dir, setup_logger, write_run_metadata
)
from modules.json_spider_adapter import prepare_annotation
from modules.mapping import map_reads_to_ref
from modules.roi_extract import extract_gene_bams
from modules.variant_call import call_variants_per_gene
from modules.annotate_mut import variants_to_aa_table
from modules.report import write_tables


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Raw data -> mapping -> gene ROI BAM -> variants -> AA mutation table",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # -------- Input / Output --------
    p.add_argument("--raw_dir", default=None, help="raw data folder (fastq/fastq.gz or pod5/fast5)")
    p.add_argument("--out_dir", required=True, help="output folder")
    p.add_argument("--sample_sheet", default=None,
                   help="CSV/TSV mapping sample_id to raw files/barcodes (optional)")
    p.add_argument("--input_bam", default=None,
                   help="Existing sorted/indexed BAM. If set, skip minimap2 mapping and use this BAM.")

    # -------- Reference / Annotation --------
    p.add_argument("--ref_fasta", required=True, help="reference genome fasta")
    p.add_argument("--ref_id", default=None,
                   help="NCBI accession or identifier for pulling annotation via json_spider.py (optional)")
    p.add_argument("--annotation", default=None,
                   help="Existing annotation file (jsonl/gff/gbk). If provided, skip json_spider.")

    # Genes / Regions selection
    p.add_argument("--genes", nargs="+", default=None,
                   help="gene names to extract (space-separated). If not set, extract all from annotation.")
    p.add_argument("--gene_bed", default=None,
                   help="BED file listing ROI regions. If provided, use this and ignore --genes parsing.")

    # -------- Workflow switches --------
    p.add_argument("--skip_basecall", action="store_true", help="skip basecalling step")
    p.add_argument("--skip_demux", action="store_true", help="skip demux step")
    p.add_argument("--skip_porechop", action="store_true", help="skip porechop step")
    p.add_argument("--skip_mapping", action="store_true", help="skip mapping step")
    p.add_argument("--skip_roi", action="store_true", help="skip ROI extraction step")
    p.add_argument("--skip_variant", action="store_true", help="skip variant calling step")
    p.add_argument("--skip_aa", action="store_true", help="skip AA annotation step")

    # -------- Tools / threads --------
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--min_mapq", type=int, default=20)
    p.add_argument("--min_depth", type=int, default=10)

    # Variant caller selection (placeholder)
    p.add_argument("--caller", default="bcftools",
                   choices=["bcftools", "freebayes", "longshot", "medaka"],
                   help="variant caller backend; initial version implements bcftools")

    # General
    p.add_argument("--dry_run", action="store_true", help="print commands without running")
    p.add_argument("--force", action="store_true", help="overwrite existing outputs")

    return p


def main():
    args = build_parser().parse_args()

    if not args.input_bam and not args.skip_mapping and not args.raw_dir:
        raise ValueError("Need --raw_dir unless --input_bam or --skip_mapping is set.")

    out_dir = Path(args.out_dir).resolve()
    ensure_dir(out_dir)
    log = setup_logger(out_dir / "pipeline.log")
    write_run_metadata(out_dir / "run_metadata.json", args)

    # Subfolders
    work = out_dir / "work"
    ensure_dir(work)
    fastq_dir = work / "fastq"
    bam_dir = work / "bam"
    roi_dir = work / "roi_bam"
    vcf_dir = work / "vcf"
    tab_dir = out_dir / "tables"
    ensure_dir(tab_dir)

    # ------------------------------------------------------------
    # STEP 0) Annotation preparation (from json_spider or input file)
    # ------------------------------------------------------------
    anno_path = prepare_annotation(
        out_dir=work,
        ref_id=args.ref_id,
        annotation=args.annotation,
        genes=args.genes,
        gene_bed=args.gene_bed,
        dry_run=args.dry_run,
        force=args.force,
        logger=log
    )
    # anno_path should point to a normalized "gene_regions.bed" (preferred) and/or "annotation.jsonl"
    # We'll assume adapter returns a dict of paths
    bed_path = anno_path["bed"]  # required later

    # ------------------------------------------------------------
    # STEP 1) Basecalling / demux / porechop (placeholders)
    # 你說先從 basecallers 開始、demux、porechop
    # 但你也說目前先不要長度篩選，先 mapping
    # -> 這裡保留骨架，預設你可以 --skip_* 直接略過
    # ------------------------------------------------------------
    if not args.skip_basecall:
        log.warning("TODO: basecalling module not implemented. Use --skip_basecall for now.")
        # from modules.basecalling import basecall
        # fastq_dir = basecall(...)
    else:
        ensure_dir(fastq_dir)

    if not args.skip_demux:
        log.warning("TODO: demux module not implemented. Use --skip_demux for now.")
        # from modules.demux import demux
        # fastq_dir = demux(...)
    else:
        ensure_dir(fastq_dir)

    if not args.skip_porechop:
        log.warning("TODO: porechop module not implemented. Use --skip_porechop for now.")
        # from modules.qc_trim import porechop
        # fastq_dir = porechop(...)
    else:
        ensure_dir(fastq_dir)

    # ------------------------------------------------------------
    # STEP 2) Mapping (minimap2 -> sort/index BAM)
    # ------------------------------------------------------------
    if args.input_bam:
        bam_path = Path(args.input_bam).resolve()
        if not bam_path.exists() and not args.dry_run:
            raise FileNotFoundError(f"input_bam not found: {bam_path}")
        log.info(f"Using existing BAM: {bam_path}")
    elif not args.skip_mapping:
        bam_path = map_reads_to_ref(
            raw_dir=Path(args.raw_dir),
            out_dir=bam_dir,
            ref_fasta=Path(args.ref_fasta),
            threads=args.threads,
            min_mapq=args.min_mapq,
            dry_run=args.dry_run,
            force=args.force,
            logger=log
        )
    else:
        log.warning("skip mapping enabled; assume BAM exists in work/bam.")
        bam_path = bam_dir / "merged.sorted.bam"

    # ------------------------------------------------------------
    # STEP 3) ROI extraction per gene (samtools view -L bed)
    # ------------------------------------------------------------
    if not args.skip_roi:
        gene_bams = extract_gene_bams(
            bam_path=bam_path,
            bed_path=Path(bed_path),
            out_dir=roi_dir,
            threads=args.threads,
            dry_run=args.dry_run,
            force=args.force,
            logger=log
        )
        # gene_bams: dict {gene: bam}
    else:
        log.warning("skip ROI enabled; discovering existing BAMs in work/roi_bam.")
        gene_bams = {p.stem: p for p in roi_dir.glob("*.bam")}

    # ------------------------------------------------------------
    # STEP 4) Variant calling per gene
    # ------------------------------------------------------------
    if not args.skip_variant:
        vcf_paths = call_variants_per_gene(
            gene_bams=gene_bams,
            ref_fasta=Path(args.ref_fasta),
            out_dir=vcf_dir,
            caller=args.caller,
            threads=args.threads,
            min_depth=args.min_depth,
            dry_run=args.dry_run,
            force=args.force,
            logger=log
        )
    else:
        log.warning("skip variant enabled; discovering existing VCFs in work/vcf.")
        vcf_paths = {p.name.split(".vcf")[0]: p for p in vcf_dir.glob("*.vcf*")}

    # ------------------------------------------------------------
    # STEP 5) Variants -> mutation candidate table
    # ------------------------------------------------------------
    if not args.skip_aa:
        aa_table = variants_to_aa_table(
            vcf_paths=vcf_paths,
            annotation_path=anno_path.get("annotation_jsonl", None),
            bed_path=Path(bed_path),
            out_dir=tab_dir,
            dry_run=args.dry_run,
            force=args.force,
            logger=log
        )
    else:
        log.warning("skip AA enabled.")
        aa_table = None

    # ------------------------------------------------------------
    # STEP 6) Report/export
    # ------------------------------------------------------------
    write_tables(
        out_dir=tab_dir,
        aa_table_path=aa_table,
        extra={},
        logger=log
    )

    log.info("Done.")


if __name__ == "__main__":
    main()
