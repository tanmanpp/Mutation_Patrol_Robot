# -*- coding: utf-8 -*-

from pathlib import Path

from modules.utils import ensure_dir, run_cmd


def call_variants_per_gene(gene_bams: dict,
                           ref_fasta: Path,
                           out_dir: Path,
                           caller: str,
                           threads: int,
                           min_depth: int,
                           dry_run: bool,
                           force: bool,
                           logger):
    """
    Call variants for each gene-level BAM.

    First implemented backend:
      bcftools mpileup + call + filter, producing compressed VCFs.
    """
    ensure_dir(out_dir)

    if caller != "bcftools":
        raise NotImplementedError("Initial version implements only --caller bcftools.")
    if not ref_fasta.exists() and not dry_run:
        raise FileNotFoundError(f"Reference FASTA not found: {ref_fasta}")

    vcf_paths = {}
    for gene, bam_path in gene_bams.items():
        bam_path = Path(bam_path)
        if not bam_path.exists() and not dry_run:
            raise FileNotFoundError(f"Gene BAM not found: {bam_path}")

        vcf_out = out_dir / f"{gene}.vcf.gz"
        if vcf_out.exists() and not force:
            logger.info(f"VCF exists, skip: {vcf_out}")
            vcf_paths[gene] = vcf_out
            continue

        raw_vcf = out_dir / f"{gene}.raw.vcf.gz"
        cmd = [
            "bash", "-lc",
            (
                f"bcftools mpileup --threads {threads} -Ou -f {ref_fasta} {bam_path} "
                f"| bcftools call --threads {threads} -mv -Oz -o {raw_vcf} && "
                f"bcftools filter -i 'DP>={min_depth}' -Oz -o {vcf_out} {raw_vcf} && "
                f"bcftools index -t {vcf_out}"
            )
        ]
        run_cmd(cmd, logger=logger, dry_run=dry_run)
        vcf_paths[gene] = vcf_out

    return vcf_paths
