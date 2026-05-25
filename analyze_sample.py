#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

from modules.annotate_mut import variants_to_aa_table
from modules.bam_allele_freq import mutation_candidates_from_bam
from modules.gene_database import load_gene_database
from modules.mapping import map_reads_to_ref
from modules.report import write_tables
from modules.roi_extract import extract_gene_bams
from modules.utils import ensure_dir, setup_logger, write_run_metadata
from modules.variant_call import call_variants_per_gene


def build_parser():
    parser = argparse.ArgumentParser(
        description="Analyze one sequencing sample against a prebuilt AMR gene database.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--gene_db", required=True, help="gene_database.json from build_gene_db.py")
    parser.add_argument("--out_dir", required=True, help="Sample analysis output folder")
    parser.add_argument("--raw_dir", default=None, help="FASTQ folder. Required unless --fastq or --input_bam is used.")
    parser.add_argument("--fastq", default=None, help="Single FASTQ/FASTQ.GZ file for one sample")
    parser.add_argument("--input_bam", default=None, help="Existing sorted/indexed sample BAM")
    parser.add_argument("--ref_fasta", default=None, help="Reference FASTA. Defaults to value stored in gene_db.")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--min_mapq", type=int, default=20)
    parser.add_argument("--min_depth", type=int, default=10)
    parser.add_argument("--min_alt_count", type=int, default=1)
    parser.add_argument("--min_alt_freq", type=float, default=0.05)
    parser.add_argument("--candidate_source", choices=["bam", "vcf"], default="bam")
    parser.add_argument("--caller", default="bcftools", choices=["bcftools", "freebayes", "longshot", "medaka"])
    parser.add_argument("--skip_mapping", action="store_true", help="Use work/bam/merged.sorted.bam")
    parser.add_argument("--skip_roi", action="store_true", help="Use existing work/roi_bam/*.bam")
    parser.add_argument("--skip_variant", action="store_true", help="Use existing work/vcf/*.vcf*")
    parser.add_argument("--dry_run", action="store_true")
    parser.add_argument("--force", action="store_true")
    return parser


def main():
    args = build_parser().parse_args()
    gene_db_path = Path(args.gene_db).resolve()
    gene_db = load_gene_database(gene_db_path)
    ref_fasta = Path(args.ref_fasta or gene_db["reference_fasta"]).resolve()

    if args.fastq and args.raw_dir:
        raise ValueError("Use either --fastq or --raw_dir, not both.")
    if not args.input_bam and not args.skip_mapping and not args.raw_dir and not args.fastq:
        raise ValueError("Need --fastq or --raw_dir unless --input_bam or --skip_mapping is set.")

    out_dir = Path(args.out_dir).resolve()
    ensure_dir(out_dir)
    logger = setup_logger(out_dir / "analyze_sample.log")
    write_run_metadata(out_dir / "analyze_sample_metadata.json", args)

    work = out_dir / "work"
    bam_dir = work / "bam"
    roi_dir = work / "roi_bam"
    vcf_dir = work / "vcf"
    tables_dir = out_dir / "tables"
    for directory in [work, bam_dir, roi_dir, vcf_dir, tables_dir]:
        ensure_dir(directory)

    bed_path = Path(gene_db_path.parent / "gene_regions.bed")
    if not bed_path.exists():
        raise FileNotFoundError(f"Expected BED next to gene database: {bed_path}")

    if args.input_bam:
        bam_path = Path(args.input_bam).resolve()
        if not bam_path.exists() and not args.dry_run:
            raise FileNotFoundError(f"input_bam not found: {bam_path}")
        logger.info(f"Using existing BAM: {bam_path}")
    elif args.skip_mapping:
        bam_path = bam_dir / "merged.sorted.bam"
        logger.warning(f"skip_mapping enabled; using expected BAM path: {bam_path}")
    else:
        bam_path = map_reads_to_ref(
            raw_dir=Path(args.raw_dir) if args.raw_dir else None,
            out_dir=bam_dir,
            ref_fasta=ref_fasta,
            threads=args.threads,
            min_mapq=args.min_mapq,
            dry_run=args.dry_run,
            force=args.force,
            logger=logger,
            read_files=[Path(args.fastq)] if args.fastq else None,
        )

    if args.skip_roi:
        gene_bams = {path.stem: path for path in roi_dir.glob("*.bam")}
        logger.warning(f"skip_roi enabled; found {len(gene_bams)} existing ROI BAMs")
    else:
        gene_bams = extract_gene_bams(
            bam_path=bam_path,
            bed_path=bed_path,
            out_dir=roi_dir,
            threads=args.threads,
            dry_run=args.dry_run,
            force=args.force,
            logger=logger,
        )

    if args.candidate_source == "bam":
        mutation_table = mutation_candidates_from_bam(
            gene_bams=gene_bams,
            gene_db_path=gene_db_path,
            ref_fasta=ref_fasta,
            out_dir=tables_dir,
            min_depth=args.min_depth,
            min_alt_count=args.min_alt_count,
            min_alt_freq=args.min_alt_freq,
            min_mapq=args.min_mapq,
            dry_run=args.dry_run,
            force=args.force,
            logger=logger,
        )
    else:
        if args.skip_variant:
            vcf_paths = {path.name.split(".vcf")[0]: path for path in vcf_dir.glob("*.vcf*")}
            logger.warning(f"skip_variant enabled; found {len(vcf_paths)} existing VCFs")
        else:
            vcf_paths = call_variants_per_gene(
                gene_bams=gene_bams,
                ref_fasta=ref_fasta,
                out_dir=vcf_dir,
                caller=args.caller,
                threads=args.threads,
                min_depth=args.min_depth,
                dry_run=args.dry_run,
                force=args.force,
                logger=logger,
            )

        mutation_table = variants_to_aa_table(
            vcf_paths=vcf_paths,
            annotation_path=None,
            bed_path=bed_path,
            out_dir=tables_dir,
            dry_run=args.dry_run,
            force=args.force,
            logger=logger,
            gene_db_path=gene_db_path,
        )
    write_tables(
        out_dir=tables_dir,
        aa_table_path=mutation_table,
        extra={
            "gene_db": str(gene_db_path),
            "gene_count": gene_db["gene_count"],
            "candidate_source": args.candidate_source,
        },
        logger=logger,
    )
    logger.info("Sample analysis complete.")


if __name__ == "__main__":
    main()
