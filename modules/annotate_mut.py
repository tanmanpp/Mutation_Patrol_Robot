# -*- coding: utf-8 -*-

import csv
import gzip
from pathlib import Path

from modules.gene_database import (
    coding_alt_base,
    genomic_pos_to_cds_index,
    load_gene_database,
    reverse_complement,
    translate_dna,
)
from modules.utils import ensure_dir


FIELDNAMES = [
    "Gene name", "Gene symbol", "gene", "chrom", "pos", "ref", "alt",
    "qual", "filter", "depth", "ref_depth", "alt_depth", "allele_freq",
    "variant_type", "source", "cds_pos", "codon_pos", "ref_codon",
    "alt_codon", "aa_pos", "ref_aa", "alt_aa", "aa_change", "note"
]


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def _info_value(info: str, key: str):
    prefix = f"{key}="
    for item in info.split(";"):
        if item.startswith(prefix):
            return item[len(prefix):]
    return ""


def _variant_type(ref: str, alt: str):
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) < len(alt):
        return "INS"
    if len(ref) > len(alt):
        return "DEL"
    return "MNV"


def _gene_display_name(gene: str, gene_db: dict | None):
    if not gene_db:
        return gene
    gene_record = gene_db.get("genes_by_symbol", {}).get(gene)
    if not gene_record:
        return gene
    return gene_record.get("description") or gene


def _annotate_record(row: dict, gene_db: dict | None):
    empty = {
        "cds_pos": "",
        "codon_pos": "",
        "ref_codon": "",
        "alt_codon": "",
        "aa_pos": "",
        "ref_aa": "",
        "alt_aa": "",
        "aa_change": "NA",
        "note": "",
    }

    if not gene_db:
        empty["note"] = "No gene database provided"
        return empty

    gene_record = gene_db.get("genes_by_symbol", {}).get(row["gene"])
    if not gene_record:
        empty["note"] = "Gene not found in database"
        return empty

    if row["variant_type"] not in ("SNV", "REF", "INS"):
        empty["note"] = "AA annotation currently supports SNVs, insertions, and reference alleles only"
        return empty

    cds_index = genomic_pos_to_cds_index(gene_record, row["chrom"], int(row["pos"]))
    if cds_index is None:
        empty["note"] = "Variant is outside annotated CDS exons"
        return empty

    cds_sequence = gene_record["cds_sequence"]
    if cds_index >= len(cds_sequence):
        empty["note"] = "CDS index is outside stored CDS sequence"
        return empty

    codon_start = (cds_index // 3) * 3
    ref_codon = cds_sequence[codon_start:codon_start + 3]
    if len(ref_codon) != 3:
        empty["note"] = "Incomplete codon"
        return empty

    if row["variant_type"] == "INS":
        inserted_sequence = _inserted_sequence(row["ref"], row["alt"], gene_record["strand"])
        ref_aa = translate_dna(ref_codon)
        aa_pos = (cds_index // 3) + 1
        codon_pos = (cds_index % 3) + 1
        inserted_aa = translate_dna(inserted_sequence) if len(inserted_sequence) % 3 == 0 else ""
        if len(inserted_sequence) % 3 == 0:
            aa_change = f"{ref_aa}{aa_pos}_ins{inserted_aa}"
            note = f"In-frame insertion from gene database; inserted_cds={inserted_sequence}"
        else:
            aa_change = f"{ref_aa}{aa_pos}fs"
            note = f"Frameshift insertion from gene database; inserted_cds={inserted_sequence}"

        return {
            "cds_pos": cds_index + 1,
            "codon_pos": codon_pos,
            "ref_codon": ref_codon,
            "alt_codon": f"{ref_codon}+{inserted_sequence}",
            "aa_pos": aa_pos,
            "ref_aa": ref_aa,
            "alt_aa": inserted_aa or "frameshift",
            "aa_change": aa_change,
            "note": note,
        }

    alt_codon_list = list(ref_codon)
    alt_base = coding_alt_base(row["ref"], row["alt"], gene_record["strand"])
    alt_codon_list[cds_index % 3] = alt_base
    alt_codon = "".join(alt_codon_list)

    ref_aa = translate_dna(ref_codon)
    alt_aa = translate_dna(alt_codon)
    aa_pos = (cds_index // 3) + 1
    aa_change = f"{ref_aa}{aa_pos}{alt_aa}"
    if ref_aa == alt_aa:
        aa_change = f"{ref_aa}{aa_pos}{ref_aa} (synonymous)"

    return {
        "cds_pos": cds_index + 1,
        "codon_pos": (cds_index % 3) + 1,
        "ref_codon": ref_codon,
        "alt_codon": alt_codon,
        "aa_pos": aa_pos,
        "ref_aa": ref_aa,
        "alt_aa": alt_aa,
        "aa_change": aa_change,
        "note": "AA annotated from gene database",
    }


def _inserted_sequence(ref: str, alt: str, strand: str):
    ref = (ref or "").upper()
    alt = (alt or "").upper()
    if alt.startswith(ref):
        inserted = alt[len(ref):]
    elif alt.startswith("+"):
        inserted = alt[1:]
    else:
        inserted = alt
    if strand == "-":
        return reverse_complement(inserted)
    return inserted


def annotate_candidate_row(row: dict, gene_db: dict | None):
    gene = row["gene"]
    row.setdefault("Gene name", _gene_display_name(gene, gene_db))
    row.setdefault("Gene symbol", gene)
    row.setdefault("qual", "")
    row.setdefault("filter", "")
    row.setdefault("ref_depth", "")
    row.setdefault("alt_depth", "")
    row.setdefault("allele_freq", "")
    row.setdefault("source", "")
    existing_note = row.get("note", "")
    annotation = _annotate_record(row, gene_db)
    if existing_note and annotation.get("note"):
        annotation["note"] = f"{existing_note}; {annotation['note']}"
    elif existing_note:
        annotation["note"] = existing_note
    row.update(annotation)
    return row


def _read_vcf_records(gene: str, vcf_path: Path, gene_db: dict | None):
    with _open_text(vcf_path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _id, ref, alts, qual, filt, info = parts[:8]
            depth = _info_value(info, "DP")
            for alt in alts.split(","):
                row = {
                    "Gene name": _gene_display_name(gene, gene_db),
                    "Gene symbol": gene,
                    "gene": gene,
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "qual": qual,
                    "filter": filt,
                    "depth": depth,
                    "ref_depth": "",
                    "alt_depth": "",
                    "allele_freq": "",
                    "variant_type": _variant_type(ref, alt),
                    "source": "vcf",
                }
                yield annotate_candidate_row(row, gene_db)


def variants_to_aa_table(vcf_paths: dict,
                         annotation_path: str | None,
                         bed_path: Path,
                         out_dir: Path,
                         dry_run: bool,
                         force: bool,
                         logger,
                         gene_db_path: str | Path | None = None):
    """
    Convert per-gene VCFs into a mutation candidate table.

    The table is intentionally shaped for later amino-acid annotation. This first
    version reports nucleotide variants and leaves aa_change as NA.
    """
    ensure_dir(out_dir)
    out_csv = out_dir / "mutation_candidates.csv"
    if out_csv.exists() and not force:
        logger.info(f"Mutation table exists, skip: {out_csv}")
        return out_csv

    if dry_run:
        logger.info(f"Dry run: would write mutation table to {out_csv}")
        return out_csv

    gene_db = load_gene_database(Path(gene_db_path)) if gene_db_path else None
    rows = []
    for gene, vcf_path in vcf_paths.items():
        vcf_path = Path(vcf_path)
        if not vcf_path.exists():
            logger.warning(f"VCF missing, skip table rows for {gene}: {vcf_path}")
            continue
        rows.extend(_read_vcf_records(gene, vcf_path, gene_db))

    with open(out_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    logger.info(f"Wrote {len(rows)} variant rows: {out_csv}")
    return out_csv
