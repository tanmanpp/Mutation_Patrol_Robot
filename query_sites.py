#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

from modules.site_query import query_bam_sites
from modules.utils import ensure_dir, setup_logger, write_run_metadata


def build_parser():
    parser = argparse.ArgumentParser(
        description="Query user-selected AMR gene positions from a BAM and output mutation_candidates-like rows.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--gene_db", required=True, help="gene_database.json from build_gene_db.py")
    parser.add_argument("--bam", required=True, help="Sorted/indexed sample BAM")
    parser.add_argument("--out_csv", required=True, help="Output CSV path")
    parser.add_argument("--gene", required=True, help="Gene symbol in the database")
    parser.add_argument("--genomic_pos", default=None, help="1-based genomic coordinate, or comma-separated coordinates")
    parser.add_argument("--cds_pos", default=None, help="1-based CDS nucleotide position(s), comma-separated")
    parser.add_argument("--aa_pos", default=None, help="1-based amino-acid position(s), comma-separated; outputs codon bases")
    parser.add_argument("--alt", default=None, help="Optional allele to report even when not observed")
    parser.add_argument("--ref_fasta", default=None, help="Reference FASTA. Defaults to value stored in gene_db.")
    parser.add_argument("--min_mapq", type=int, default=20)
    parser.add_argument("--min_alt_freq", type=float, default=0.05,
                        help="Minimum non-reference allele frequency to report when --alt is not set")
    parser.add_argument("--force", action="store_true")
    return parser


def main():
    args = build_parser().parse_args()
    out_csv = Path(args.out_csv).resolve()
    ensure_dir(out_csv.parent)
    logger = setup_logger(out_csv.parent / "query_sites.log")
    write_run_metadata(out_csv.parent / "query_sites_metadata.json", args)

    query_bam_sites(
        gene_db_path=Path(args.gene_db).resolve(),
        bam_path=Path(args.bam).resolve(),
        out_csv=out_csv,
        gene=args.gene,
        genomic_pos=parse_positions(args.genomic_pos),
        cds_pos=parse_positions(args.cds_pos),
        aa_pos=parse_positions(args.aa_pos),
        alt=args.alt,
        ref_fasta=Path(args.ref_fasta).resolve() if args.ref_fasta else None,
        min_mapq=args.min_mapq,
        min_alt_freq=args.min_alt_freq,
        force=args.force,
        logger=logger,
    )
    logger.info("Site query complete.")


def parse_positions(value):
    if value is None:
        return None
    positions = []
    for item in str(value).split(","):
        item = item.strip()
        if not item:
            continue
        positions.append(int(item))
    return positions or None


if __name__ == "__main__":
    main()
