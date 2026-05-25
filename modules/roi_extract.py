# -*- coding: utf-8 -*-

from pathlib import Path

from modules.utils import ensure_dir, run_cmd


def _load_bed(bed_path: Path):
    """
    Expect BED4+:
      chrom  start0  end0  gene
    """
    items = []
    with open(bed_path, "r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                raise ValueError(f"BED line {line_no} must have at least 4 columns")
            chrom, start0, end0, gene = parts[0], int(parts[1]), int(parts[2]), parts[3]
            items.append((chrom, start0, end0, gene))
    return items


def extract_gene_bams(bam_path: Path,
                      bed_path: Path,
                      out_dir: Path,
                      threads: int,
                      dry_run: bool,
                      force: bool,
                      logger):
    ensure_dir(out_dir)

    if not bam_path.exists() and not dry_run:
        raise FileNotFoundError(str(bam_path))
    if not bed_path.exists():
        raise FileNotFoundError(str(bed_path))

    regions = _load_bed(bed_path)
    gene_bams = {}

    for chrom, start0, end0, gene in regions:
        out_bam = out_dir / f"{gene}.bam"
        if out_bam.exists() and not force:
            logger.info(f"ROI BAM exists, skip: {out_bam}")
            gene_bams[gene] = out_bam
            continue

        # samtools region strings are 1-based inclusive; BED is 0-based half-open.
        region_str = f"{chrom}:{start0 + 1}-{end0}"
        cmd = [
            "bash", "-lc",
            (
                f"samtools view -@ {threads} -b {bam_path} {region_str} "
                f"| samtools sort -@ {threads} -o {out_bam} - && "
                f"samtools index {out_bam}"
            )
        ]
        run_cmd(cmd, logger=logger, dry_run=dry_run)
        gene_bams[gene] = out_bam

    return gene_bams
