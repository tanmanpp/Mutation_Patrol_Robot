# -*- coding: utf-8 -*-

from pathlib import Path

from modules.utils import ensure_dir, run_cmd


def bam_index_candidates(bam_path: Path):
    bam_path = Path(bam_path)
    return [
        Path(str(bam_path) + ".bai"),
        bam_path.with_suffix(".bai"),
    ]


def fasta_index_path(fasta_path: Path):
    return Path(str(Path(fasta_path)) + ".fai")


def ensure_fasta_index(fasta_path: Path, logger, dry_run: bool = False):
    fasta_path = Path(fasta_path)
    index_path = fasta_index_path(fasta_path)
    if (
        index_path.exists()
        and index_path.stat().st_mtime >= fasta_path.stat().st_mtime
    ):
        return index_path
    run_cmd(
        ["samtools", "faidx", str(fasta_path)],
        logger=logger,
        dry_run=dry_run,
    )
    return index_path


def ensure_bam_index(bam_path: Path, logger, dry_run: bool = False):
    bam_path = Path(bam_path)
    if any(path.exists() for path in bam_index_candidates(bam_path)):
        return
    run_cmd(
        ["samtools", "index", str(bam_path)],
        logger=logger,
        dry_run=dry_run,
    )


def prepare_sorted_bam(source_bam: Path,
                       dest_bam: Path,
                       threads: int,
                       force: bool,
                       logger,
                       dry_run: bool = False):
    """Copy/sort an input BAM into the canonical run location and index it."""
    source_bam = Path(source_bam)
    dest_bam = Path(dest_bam)
    ensure_dir(dest_bam.parent)

    if not source_bam.exists() and not dry_run:
        raise FileNotFoundError(f"Input BAM not found: {source_bam}")

    same_path = source_bam.resolve() == dest_bam.resolve()
    if same_path:
        ensure_bam_index(dest_bam, logger=logger, dry_run=dry_run)
        return dest_bam

    if dest_bam.exists() and not force:
        ensure_bam_index(dest_bam, logger=logger, dry_run=dry_run)
        return dest_bam

    temp_bam = dest_bam.with_name(f"{dest_bam.stem}.preparing.bam")
    run_cmd(
        [
            "samtools", "sort",
            "-@", str(threads),
            "-o", str(temp_bam),
            str(source_bam),
        ],
        logger=logger,
        dry_run=dry_run,
    )
    if not dry_run:
        temp_bam.replace(dest_bam)
    ensure_bam_index(dest_bam, logger=logger, dry_run=dry_run)
    return dest_bam
