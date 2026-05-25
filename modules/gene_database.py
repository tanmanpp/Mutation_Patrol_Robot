# -*- coding: utf-8 -*-

import csv
import json
from datetime import datetime
from pathlib import Path

from modules.utils import ensure_dir, relative_path


SUMMARY_FIELDNAMES = [
    "symbol", "gene_id", "description", "chrom", "resolved_chrom", "strand",
    "gene_start", "gene_end", "cds_length", "aa_length", "exon_count",
    "sequence_status"
]

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def load_fasta(fasta_path: Path):
    sequences = {}
    name = None
    chunks = []

    with open(fasta_path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    sequences[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

    if name:
        sequences[name] = "".join(chunks).upper()

    if not sequences:
        raise ValueError(f"No sequences found in FASTA: {fasta_path}")
    return sequences


def reverse_complement(sequence: str):
    return sequence.translate(COMPLEMENT)[::-1].upper()


def translate_dna(sequence: str):
    amino_acids = []
    upper = sequence.upper()
    for i in range(0, len(upper) - 2, 3):
        codon = upper[i:i + 3]
        amino_acids.append(CODON_TABLE.get(codon, "X"))
    return "".join(amino_acids)


def wrap_sequence(sequence: str, width: int = 60):
    return "\n".join(sequence[i:i + width] for i in range(0, len(sequence), width))


def normalize_strand(orientation):
    if orientation is None:
        return "."
    if isinstance(orientation, int):
        return "+" if orientation >= 0 else "-"
    if isinstance(orientation, str):
        value = orientation.strip().lower()
        if value in ("plus", "+", "forward", "fwd", "1", "positive"):
            return "+"
        if value in ("minus", "-", "reverse", "rev", "-1", "negative"):
            return "-"
    return "."


def _strip_accession_version(accession: str):
    if not accession:
        return accession
    return accession.rsplit(".", 1)[0]


def _resolve_reference_chrom(reference: dict, chrom: str):
    if chrom in reference:
        return chrom, None

    chrom_base = _strip_accession_version(chrom)
    matches = [
        ref_name for ref_name in reference
        if _strip_accession_version(ref_name) == chrom_base
    ]
    if len(matches) == 1:
        warning = (
            f"Chromosome ID '{chrom}' not found exactly; using reference "
            f"contig '{matches[0]}' by version-insensitive match."
        )
        return matches[0], warning

    if len(matches) > 1:
        warning = (
            f"Chromosome ID '{chrom}' matched multiple reference contigs by "
            f"version-insensitive match: {', '.join(matches)}. Sequence was not extracted."
        )
        return None, warning

    warning = (
        f"Chromosome ID '{chrom}' was not found in reference FASTA. "
        "Coordinates were kept, but CDS/protein sequence was not extracted."
    )
    return None, warning


def _extract_reference_slice(reference: dict, resolved_chrom: str, start: int, end: int):
    if start < 1 or end > len(reference[resolved_chrom]) or start > end:
        raise ValueError(f"Invalid interval for {resolved_chrom}: {start}-{end}")
    return reference[resolved_chrom][start - 1:end]


def _parse_ncbi_product_report_line(line: str, reference: dict):
    data = json.loads(line)
    transcripts = data.get("transcripts") or []
    if not transcripts:
        return None

    transcript = transcripts[0]
    locations = transcript.get("genomicLocations") or []
    if not locations:
        return None

    location = locations[0]
    raw_exons = location.get("exons") or []
    if not raw_exons:
        return None

    genomic_range = location.get("genomicRange") or {}
    symbol = data.get("symbol", "N/A")
    chrom = location.get("genomicAccessionVersion")
    strand = normalize_strand(genomic_range.get("orientation"))
    warnings = []
    resolved_chrom, chrom_warning = _resolve_reference_chrom(reference, chrom)
    if chrom_warning:
        warnings.append(chrom_warning)

    exons = []
    for exon in sorted(raw_exons, key=lambda item: int(item.get("order", 0))):
        start = int(exon["begin"])
        end = int(exon["end"])
        exons.append({
            "order": int(exon.get("order", 0)),
            "start": start,
            "end": end,
            "length": end - start + 1,
        })

    exon_min = min(exon["start"] for exon in exons)
    exon_max = max(exon["end"] for exon in exons)
    gr_begin = genomic_range.get("begin")
    gr_end = genomic_range.get("end")
    gene_start = exon_min if gr_begin is None else min(int(gr_begin), exon_min)
    gene_end = exon_max if gr_end is None else max(int(gr_end), exon_max)

    offset = 0
    for exon in exons:
        exon["cds_offset"] = offset
        offset += exon["length"]

    cds_chunks = []
    sequence_status = "ok"
    if resolved_chrom:
        for exon in exons:
            try:
                seq = _extract_reference_slice(reference, resolved_chrom, exon["start"], exon["end"])
            except ValueError as exc:
                warnings.append(f"{symbol}: {exc}. Sequence was not extracted.")
                cds_chunks = []
                sequence_status = "sequence_unavailable"
                break
            if strand == "-":
                seq = reverse_complement(seq)
            cds_chunks.append(seq)
    else:
        sequence_status = "sequence_unavailable"

    cds_sequence = "".join(cds_chunks).upper()
    aa_sequence = translate_dna(cds_sequence)
    if not cds_sequence and sequence_status == "ok":
        sequence_status = "sequence_unavailable"

    return {
        "symbol": symbol,
        "gene_id": data.get("geneId", ""),
        "description": data.get("description", ""),
        "tax_id": data.get("taxId", ""),
        "tax_name": data.get("taxname", ""),
        "transcript_accession": transcript.get("accessionVersion", ""),
        "protein_accession": (transcript.get("protein") or {}).get("accessionVersion", ""),
        "chrom": chrom,
        "resolved_chrom": resolved_chrom or "",
        "strand": strand,
        "gene_start": gene_start,
        "gene_end": gene_end,
        "bed_start0": gene_start - 1,
        "bed_end0": gene_end,
        "exons": exons,
        "cds_length": len(cds_sequence),
        "aa_length": len(aa_sequence.rstrip("*")),
        "cds_sequence": cds_sequence,
        "aa_sequence": aa_sequence,
        "sequence_status": sequence_status,
        "warnings": warnings,
    }


def parse_ncbi_product_report(annotation_path: Path,
                              ref_fasta: Path,
                              genes: list[str] | None = None):
    reference = load_fasta(ref_fasta)
    selected = set(genes or [])
    records = []

    with open(annotation_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            record = _parse_ncbi_product_report_line(line, reference)
            if not record:
                continue
            if selected and record["symbol"] not in selected:
                continue
            records.append(record)

    if not records:
        raise ValueError(f"No matching genes parsed from {annotation_path}")
    return records


def parse_ncbi_product_reports(annotation_paths: list[Path],
                               ref_fasta: Path,
                               genes: list[str] | None = None):
    reference = load_fasta(ref_fasta)
    selected = set(genes or [])
    records = []
    seen_symbols = {}

    for annotation_path in annotation_paths:
        with open(annotation_path, "r", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                record = _parse_ncbi_product_report_line(line, reference)
                if not record:
                    continue
                if selected and record["symbol"] not in selected:
                    continue

                symbol = record["symbol"]
                if symbol in seen_symbols:
                    first_path = seen_symbols[symbol]
                    raise ValueError(
                        f"Duplicate gene symbol '{symbol}' found in both "
                        f"{first_path} and {annotation_path}"
                    )
                seen_symbols[symbol] = annotation_path
                record["source_annotation"] = relative_path(annotation_path)
                records.append(record)

    if not records:
        joined = ", ".join(str(path) for path in annotation_paths)
        raise ValueError(f"No matching genes parsed from annotation inputs: {joined}")
    return records


def write_gene_database(records: list[dict],
                        ref_fasta: Path,
                        annotation_paths,
                        out_dir: Path,
                        force: bool = False):
    ensure_dir(out_dir)

    json_path = out_dir / "gene_database.json"
    csv_path = out_dir / "gene_database.csv"
    bed_path = out_dir / "gene_regions.bed"
    cds_fasta_path = out_dir / "gene_cds.fna"
    aa_fasta_path = out_dir / "gene_proteins.faa"
    warnings_path = out_dir / "gene_database_warnings.txt"

    output_paths = [json_path, csv_path, bed_path, cds_fasta_path, aa_fasta_path, warnings_path]
    if any(path.exists() for path in output_paths) and not force:
        existing = ", ".join(str(path) for path in output_paths if path.exists())
        raise FileExistsError(f"Output exists; use --force to overwrite: {existing}")

    database_warnings = []
    for record in records:
        for warning in record.get("warnings", []):
            database_warnings.append(f'{record["symbol"]}: {warning}')

    if isinstance(annotation_paths, (str, Path)):
        annotation_list = [Path(annotation_paths)]
    else:
        annotation_list = [Path(path) for path in annotation_paths]

    payload = {
        "created_at": datetime.now().isoformat(),
        "reference_fasta": relative_path(ref_fasta),
        "annotations": [relative_path(path) for path in annotation_list],
        "gene_count": len(records),
        "warning_count": len(database_warnings),
        "warnings": database_warnings,
        "genes": records,
    }

    with open(json_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)

    with open(csv_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDNAMES)
        writer.writeheader()
        for record in records:
            writer.writerow({
                "symbol": record["symbol"],
                "gene_id": record["gene_id"],
                "description": record["description"],
                "chrom": record["chrom"],
                "resolved_chrom": record.get("resolved_chrom", ""),
                "strand": record["strand"],
                "gene_start": record["gene_start"],
                "gene_end": record["gene_end"],
                "cds_length": record["cds_length"],
                "aa_length": record["aa_length"],
                "exon_count": len(record["exons"]),
                "sequence_status": record.get("sequence_status", ""),
            })

    with open(bed_path, "w", encoding="utf-8", newline="") as handle:
        for record in records:
            bed_chrom = record.get("resolved_chrom") or record["chrom"]
            handle.write(
                f'{bed_chrom}\t{record["bed_start0"]}\t{record["bed_end0"]}\t'
                f'{record["symbol"]}\t0\t{record["strand"]}\n'
            )

    with open(cds_fasta_path, "w", encoding="utf-8") as handle:
        for record in records:
            header = f'{record["symbol"]}|{record["chrom"]}:{record["gene_start"]}-{record["gene_end"]}({record["strand"]})'
            handle.write(f">{header}\n{wrap_sequence(record['cds_sequence'])}\n")

    with open(aa_fasta_path, "w", encoding="utf-8") as handle:
        for record in records:
            header = f'{record["symbol"]}|{record.get("protein_accession", "")}'
            handle.write(f">{header}\n{wrap_sequence(record['aa_sequence'])}\n")

    with open(warnings_path, "w", encoding="utf-8") as handle:
        for warning in database_warnings:
            handle.write(warning + "\n")

    return {
        "json": json_path,
        "csv": csv_path,
        "bed": bed_path,
        "cds_fasta": cds_fasta_path,
        "protein_fasta": aa_fasta_path,
        "warnings": warnings_path,
    }


def build_gene_database(annotation_path: Path,
                        ref_fasta: Path,
                        out_dir: Path,
                        genes: list[str] | None = None,
                        force: bool = False):
    return build_gene_database_from_annotations(
        annotation_paths=[annotation_path],
        ref_fasta=ref_fasta,
        out_dir=out_dir,
        genes=genes,
        force=force,
    )


def build_gene_database_from_annotations(annotation_paths: list[Path],
                                         ref_fasta: Path,
                                         out_dir: Path,
                                         genes: list[str] | None = None,
                                         force: bool = False):
    records = parse_ncbi_product_reports(annotation_paths, ref_fasta, genes)
    return write_gene_database(records, ref_fasta, annotation_paths, out_dir, force)


def load_gene_database(db_path: Path):
    with open(db_path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    genes = payload.get("genes") or []
    payload["genes_by_symbol"] = {record["symbol"]: record for record in genes}
    return payload


def genomic_pos_to_cds_index(gene_record: dict, chrom: str, pos: int):
    accepted_chroms = {gene_record["chrom"]}
    if gene_record.get("resolved_chrom"):
        accepted_chroms.add(gene_record["resolved_chrom"])
    if chrom not in accepted_chroms:
        return None

    for exon in gene_record["exons"]:
        start = int(exon["start"])
        end = int(exon["end"])
        if start <= pos <= end:
            offset = int(exon["cds_offset"])
            if gene_record["strand"] == "-":
                return offset + (end - pos)
            return offset + (pos - start)
    return None


def cds_index_to_genomic_pos(gene_record: dict, cds_index: int):
    if cds_index < 0:
        return None

    for exon in gene_record["exons"]:
        start = int(exon["start"])
        end = int(exon["end"])
        offset = int(exon["cds_offset"])
        length = int(exon["length"])
        if offset <= cds_index < offset + length:
            exon_offset = cds_index - offset
            if gene_record["strand"] == "-":
                return end - exon_offset
            return start + exon_offset
    return None


def coding_alt_base(ref_base: str, alt_base: str, strand: str):
    if strand == "-":
        return reverse_complement(alt_base)
    return alt_base.upper()
