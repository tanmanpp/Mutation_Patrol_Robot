#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

from modules.gene_database import build_gene_database_from_annotations
from modules.utils import ensure_dir, setup_logger, write_run_metadata


def build_parser():
    parser = argparse.ArgumentParser(
        description="Build an AMR gene database from a reference genome and NCBI product_report JSONL.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--ref_fasta", required=True, help="Reference genome FASTA")
    parser.add_argument(
        "--annotation",
        nargs="+",
        default=None,
        help="One or more NCBI product_report JSONL files",
    )
    parser.add_argument(
        "--annotation_dir",
        default=None,
        help="Folder containing per-gene .jsonl files to merge into one database",
    )
    parser.add_argument("--out_dir", required=True, help="Output database folder")
    parser.add_argument("--genes", nargs="+", default=None, help="Optional gene symbols to include")
    parser.add_argument("--force", action="store_true", help="Overwrite existing database outputs")
    return parser


def collect_annotation_paths(annotation_args, annotation_dir):
    paths = []
    for raw_path in annotation_args or []:
        path = Path(raw_path)
        if not path.exists():
            raise FileNotFoundError(f"annotation not found: {path}")
        paths.append(path)

    if annotation_dir:
        directory = Path(annotation_dir)
        if not directory.exists():
            raise FileNotFoundError(f"annotation_dir not found: {directory}")
        if not directory.is_dir():
            raise NotADirectoryError(str(directory))
        paths.extend(sorted(directory.glob("*.jsonl")))

    if not paths:
        raise ValueError("Need --annotation and/or --annotation_dir")

    unique_paths = []
    seen = set()
    for path in paths:
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        unique_paths.append(resolved)
    return unique_paths


def main():
    args = build_parser().parse_args()
    out_dir = Path(args.out_dir).resolve()
    ensure_dir(out_dir)
    logger = setup_logger(out_dir / "build_gene_db.log")
    write_run_metadata(out_dir / "build_gene_db_metadata.json", args)

    annotation_paths = collect_annotation_paths(args.annotation, args.annotation_dir)
    logger.info(f"Building database from {len(annotation_paths)} annotation file(s).")

    outputs = build_gene_database_from_annotations(
        annotation_paths=annotation_paths,
        ref_fasta=Path(args.ref_fasta),
        out_dir=out_dir,
        genes=args.genes,
        force=args.force,
    )

    for label, path in outputs.items():
        logger.info(f"{label}: {path}")
    warnings_path = outputs.get("warnings")
    if warnings_path and Path(warnings_path).exists():
        warnings = [
            line.strip()
            for line in Path(warnings_path).read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        for warning in warnings:
            logger.warning(warning)
        if warnings:
            logger.warning(f"Gene database built with {len(warnings)} warning(s).")
    logger.info("Gene database build complete.")


if __name__ == "__main__":
    main()
