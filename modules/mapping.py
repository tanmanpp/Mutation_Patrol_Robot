# -*- coding: utf-8 -*-

from pathlib import Path

from modules.utils import ensure_dir, run_cmd


def _find_fastq_files(raw_dir: Path):
    exts = ["*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"]
    files = []
    for ext in exts:
        files.extend(raw_dir.rglob(ext))
    return sorted(files)


def _validate_fastq_files(read_files):
    files = [Path(path) for path in read_files]
    for path in files:
        if not path.exists():
            raise FileNotFoundError(f"FASTQ file not found: {path}")
        if not path.is_file():
            raise ValueError(f"FASTQ input is not a file: {path}")
    return files


def map_reads_to_ref(raw_dir: Path,
                     out_dir: Path,
                     ref_fasta: Path,
                     threads: int,
                     min_mapq: int,
                     dry_run: bool,
                     force: bool,
                     logger,
                     read_files=None):
    """
    Map FASTQ reads to a reference with minimap2 and produce a sorted BAM.

    Initial assumptions:
      - Nanopore reads, using minimap2 preset map-ont.
      - All FASTQ inputs belong to one sample/run.
    """
    ensure_dir(out_dir)

    if not ref_fasta.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {ref_fasta}")

    if read_files:
        fq_files = _validate_fastq_files(read_files)
    else:
        if raw_dir is None:
            raise ValueError("Need raw_dir or read_files for mapping.")
        if not raw_dir.exists():
            raise FileNotFoundError(f"Raw read folder not found: {raw_dir}")
        fq_files = _find_fastq_files(raw_dir)

    if not fq_files:
        raise FileNotFoundError(f"No FASTQ input found.")

    bam_out = out_dir / "merged.sorted.bam"
    if bam_out.exists() and not force:
        logger.info(f"BAM exists, skip: {bam_out}")
        return bam_out

    read_inputs = " ".join(str(p) for p in fq_files)
    cmd = [
        "bash", "-lc",
        (
            f"minimap2 -t {threads} -a -x map-ont {ref_fasta} {read_inputs} "
            f"| samtools view -b -q {min_mapq} - "
            f"| samtools sort -@ {threads} -o {bam_out} - && "
            f"samtools index {bam_out}"
        )
    ]
    run_cmd(cmd, logger=logger, dry_run=dry_run)
    return bam_out
