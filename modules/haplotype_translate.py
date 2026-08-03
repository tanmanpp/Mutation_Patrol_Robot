# -*- coding: utf-8 -*-

import csv
from collections import Counter, defaultdict
from pathlib import Path

from modules.gene_database import (
    cds_index_to_genomic_pos,
    reverse_complement,
    translate_dna,
)
from modules.utils import run_cmd


HAPLOTYPE_FIELDNAMES = [
    "gene",
    "gene_name",
    "chrom",
    "cluster_id",
    "cluster_start",
    "cluster_end",
    "cluster_events",
    "haplotype_id",
    "spanning_read_count",
    "haplotype_read_count",
    "haplotype_frequency",
    "forward_depth",
    "reverse_depth",
    "haplotype_spectrum",
    "call_status",
    "phase_status",
    "qc_flags",
    "frame_status",
    "effect",
    "net_nt_change",
    "cds_window_start",
    "cds_window_end",
    "genomic_window_start",
    "genomic_window_end",
    "reference_cds_haplotype",
    "sample_cds_haplotype",
    "raw_alignment_events",
    "reference_peptide",
    "sample_peptide",
    "aa_start",
    "aa_end",
    "aa_changes",
    "combined_aa_change",
    "premature_stop",
    "note",
]

CIGAR_OPERATIONS = {"M", "I", "D", "N", "S", "H", "P", "=", "X"}


def _optional_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _parse_cigar(cigar: str):
    number = []
    operations = []
    for char in cigar:
        if char.isdigit():
            number.append(char)
            continue
        if char not in CIGAR_OPERATIONS or not number:
            return []
        operations.append((int("".join(number)), char))
        number = []
    return operations if not number else []


def _coding_indel_clusters(variant_rows: list[dict], max_cds_gap: int = 12):
    candidates = []
    for row in variant_rows:
        if row.get("variant_type") not in {"INS", "DEL"}:
            continue
        if row.get("region_type") != "CDS":
            continue
        cds_pos = _optional_int(row.get("cds_pos"))
        genomic_pos = _optional_int(row.get("pos"))
        if cds_pos is None or genomic_pos is None:
            continue
        length_change = abs(len(row.get("alt") or "") - len(row.get("ref") or ""))
        candidates.append({
            "row": row,
            "cds_index": cds_pos - 1,
            "genomic_pos": genomic_pos,
            "length_change": max(1, length_change),
        })

    candidates.sort(key=lambda item: item["cds_index"])
    clusters = []
    for candidate in candidates:
        if (
            not clusters
            or candidate["cds_index"] - clusters[-1][-1]["cds_index"] > max_cds_gap
        ):
            clusters.append([candidate])
        else:
            clusters[-1].append(candidate)
    return clusters


def _run_samtools_view(bam_path: Path,
                       chrom: str,
                       start: int,
                       end: int,
                       min_mapq: int,
                       logger):
    region = f"{chrom}:{start}-{end}"
    try:
        completed = run_cmd(
            [
                "samtools",
                "view",
                "-q",
                str(min_mapq),
                "-F",
                "2308",
                str(bam_path),
                region,
            ],
            logger=logger,
            capture_output=True,
        )
    except Exception as exc:
        raise RuntimeError(
            "samtools view failed during protein haplotype reconstruction. "
            "Check the BAM index and contig identifiers."
        ) from exc
    return completed.stdout.splitlines()


def _base_quality(qualities: str, index: int):
    if qualities == "*" or index < 0 or index >= len(qualities):
        return 0
    return ord(qualities[index]) - 33


def _read_coding_haplotype(sam_line: str,
                           coding_positions: list[int],
                           strand: str,
                           min_baseq: int):
    parts = sam_line.rstrip("\n").split("\t")
    if len(parts) < 11:
        return None

    flag = int(parts[1])
    ref_pos = int(parts[3])
    cigar = parts[5]
    sequence = parts[9].upper()
    qualities = parts[10]
    operations = _parse_cigar(cigar)
    if not operations or not sequence or sequence == "*":
        return None

    target = set(coding_positions)
    target_min = min(target)
    target_max = max(target)
    read_pos = 0
    bases = {}
    deletions = set()
    insertions = defaultdict(list)
    raw_events = []

    for length, operation in operations:
        if operation in {"M", "=", "X"}:
            for offset in range(length):
                current_ref = ref_pos + offset
                if current_ref not in target:
                    continue
                current_read = read_pos + offset
                if current_read >= len(sequence):
                    return None
                if _base_quality(qualities, current_read) < min_baseq:
                    return None
                bases[current_ref] = sequence[current_read]
            ref_pos += length
            read_pos += length
        elif operation == "I":
            anchor = ref_pos - 1
            inserted = sequence[read_pos:read_pos + length]
            relevant = anchor in target or anchor + 1 in target
            if relevant:
                inserted_qualities = [
                    _base_quality(qualities, index)
                    for index in range(read_pos, read_pos + length)
                ]
                if (
                    len(inserted) != length
                    or any(value < min_baseq for value in inserted_qualities)
                ):
                    return None
                insertions[anchor].append(inserted)
                raw_events.append(f"{anchor}:+{inserted}")
            read_pos += length
        elif operation == "D":
            deleted_positions = set(range(ref_pos, ref_pos + length))
            relevant_positions = deleted_positions & target
            if relevant_positions:
                deletions.update(relevant_positions)
                raw_events.append(f"{ref_pos - 1}:-{length}")
            ref_pos += length
        elif operation == "N":
            ref_pos += length
        elif operation == "S":
            read_pos += length

    if any(pos not in bases and pos not in deletions for pos in target):
        return None

    genomic_sequence = []
    for pos in sorted(coding_positions):
        if pos not in deletions:
            genomic_sequence.append(bases[pos])
        for inserted in insertions.get(pos, []):
            genomic_sequence.append(inserted)

    coding_sequence = "".join(genomic_sequence)
    if strand == "-":
        coding_sequence = reverse_complement(coding_sequence)
    if not coding_sequence or not set(coding_sequence).issubset({"A", "C", "G", "T"}):
        return None

    return {
        "sequence": coding_sequence,
        "events": ";".join(raw_events) or "NO_INDEL",
        "strand": "reverse" if flag & 16 else "forward",
        "target_min": target_min,
        "target_max": target_max,
    }


def _common_prefix_length(reference: str, sample: str):
    length = 0
    for ref_value, sample_value in zip(reference, sample):
        if ref_value != sample_value:
            break
        length += 1
    return length


def _common_suffix_length(reference: str, sample: str, prefix_length: int):
    maximum = min(len(reference), len(sample)) - prefix_length
    length = 0
    while length < maximum and reference[-(length + 1)] == sample[-(length + 1)]:
        length += 1
    return length


def _protein_consequence(reference_cds: str,
                         sample_cds: str,
                         has_alignment_indel: bool):
    reference_protein = translate_dna(reference_cds)
    sample_protein = translate_dna(sample_cds)
    prefix = _common_prefix_length(reference_protein, sample_protein)
    net_nt_change = len(sample_cds) - len(reference_cds)

    if reference_protein == sample_protein:
        peptide_start = max(0, prefix - 3)
        peptide_end = min(len(reference_protein), peptide_start + 12)
        return {
            "frame_status": "FRAME_RESTORED" if has_alignment_indel else "IN_FRAME",
            "effect": "synonymous_or_alignment_equivalent",
            "net_nt_change": net_nt_change,
            "reference_peptide": reference_protein[peptide_start:peptide_end],
            "sample_peptide": sample_protein[peptide_start:peptide_end],
            "aa_start": "",
            "aa_end": "",
            "aa_changes": "No amino-acid change",
            "combined_aa_change": "No amino-acid change",
            "premature_stop": "",
        }

    first_position = prefix + 1
    peptide_start = max(0, prefix - 2)

    if net_nt_change % 3 != 0:
        reference_aa = reference_protein[prefix] if prefix < len(reference_protein) else "?"
        stop_index = sample_protein.find("*", prefix)
        premature_stop = stop_index + 1 if stop_index >= 0 else ""
        stop_suffix = (
            f"*{stop_index - prefix + 1}"
            if stop_index >= 0 else ""
        )
        peptide_end = (
            min(len(sample_protein), stop_index + 1)
            if stop_index >= 0 else min(len(sample_protein), prefix + 20)
        )
        return {
            "frame_status": "PERSISTENT_FRAMESHIFT",
            "effect": "frameshift",
            "net_nt_change": net_nt_change,
            "reference_peptide": reference_protein[
                peptide_start:min(len(reference_protein), prefix + 20)
            ],
            "sample_peptide": sample_protein[peptide_start:peptide_end],
            "aa_start": first_position,
            "aa_end": "",
            "aa_changes": f"{reference_aa}{first_position}fs{stop_suffix}",
            "combined_aa_change": f"{reference_aa}{first_position}fs{stop_suffix}",
            "premature_stop": premature_stop,
        }

    suffix = _common_suffix_length(reference_protein, sample_protein, prefix)
    ref_end = len(reference_protein) - suffix if suffix else len(reference_protein)
    alt_end = len(sample_protein) - suffix if suffix else len(sample_protein)
    ref_segment = reference_protein[prefix:ref_end]
    alt_segment = sample_protein[prefix:alt_end]
    aa_end = prefix + len(ref_segment)

    changes = []
    if len(ref_segment) == len(alt_segment):
        for offset, (ref_aa, alt_aa) in enumerate(zip(ref_segment, alt_segment)):
            if ref_aa != alt_aa:
                changes.append(f"{ref_aa}{first_position + offset}{alt_aa}")

    if not ref_segment:
        left_pos = max(1, first_position - 1)
        left_aa = reference_protein[left_pos - 1] if reference_protein else "?"
        right_aa = reference_protein[prefix] if prefix < len(reference_protein) else "?"
        combined = f"{left_aa}{left_pos}_{right_aa}{first_position}ins{alt_segment}"
        effect = "inframe_insertion"
    elif not alt_segment:
        combined = (
            f"{ref_segment[0]}{first_position}_"
            f"{ref_segment[-1]}{aa_end}del"
        )
        effect = "inframe_deletion"
    elif len(ref_segment) == 1 and len(alt_segment) == 1:
        combined = f"{ref_segment}{first_position}{alt_segment}"
        effect = "missense"
    else:
        combined = (
            f"{ref_segment[0]}{first_position}_"
            f"{ref_segment[-1]}{aa_end}delins{alt_segment}"
        )
        effect = "frame_restored_complex" if has_alignment_indel else "complex_substitution"

    peptide_end = min(
        max(len(reference_protein), len(sample_protein)),
        max(ref_end, alt_end) + 8,
    )
    return {
        "frame_status": "FRAME_RESTORED" if has_alignment_indel else "IN_FRAME",
        "effect": effect,
        "net_nt_change": net_nt_change,
        "reference_peptide": reference_protein[peptide_start:peptide_end],
        "sample_peptide": sample_protein[peptide_start:peptide_end],
        "aa_start": first_position,
        "aa_end": aa_end,
        "aa_changes": "; ".join(changes) or combined,
        "combined_aa_change": combined,
        "premature_stop": "",
    }


def _no_call_row(gene: str,
                 gene_record: dict,
                 cluster_id: str,
                 cluster: list[dict],
                 window_start: int,
                 window_end: int,
                 genomic_positions: list[int],
                 spanning_depth: int,
                 min_depth: int):
    cluster_positions = [item["genomic_pos"] for item in cluster]
    return {
        "gene": gene,
        "gene_name": gene_record.get("description") or gene,
        "chrom": gene_record.get("resolved_chrom") or gene_record["chrom"],
        "cluster_id": cluster_id,
        "cluster_start": min(cluster_positions),
        "cluster_end": max(cluster_positions),
        "cluster_events": "; ".join(
            f"{item['row'].get('pos')}:{item['row'].get('ref')}>{item['row'].get('alt')}"
            for item in cluster
        ),
        "haplotype_id": "",
        "spanning_read_count": spanning_depth,
        "haplotype_read_count": 0,
        "haplotype_frequency": "0.000000",
        "forward_depth": 0,
        "reverse_depth": 0,
        "haplotype_spectrum": "",
        "call_status": "NO_CALL",
        "phase_status": "PHASE_UNRESOLVED",
        "qc_flags": "LOW_SPANNING_DEPTH",
        "frame_status": "UNRESOLVED",
        "effect": "no_call",
        "net_nt_change": "",
        "cds_window_start": window_start + 1,
        "cds_window_end": window_end + 1,
        "genomic_window_start": min(genomic_positions),
        "genomic_window_end": max(genomic_positions),
        "reference_cds_haplotype": gene_record["cds_sequence"][window_start:window_end + 1],
        "sample_cds_haplotype": "",
        "raw_alignment_events": "",
        "reference_peptide": "",
        "sample_peptide": "",
        "aa_start": "",
        "aa_end": "",
        "aa_changes": "Unable to determine",
        "combined_aa_change": "Unable to determine",
        "premature_stop": "",
        "note": (
            f"Only {spanning_depth} quality-filtered reads reconstructed the complete "
            f"cluster window; minimum required depth is {min_depth}."
        ),
    }


def _haplotype_rows_for_cluster(gene: str,
                                gene_record: dict,
                                cluster_id: str,
                                cluster: list[dict],
                                sam_lines: list[str],
                                min_alt_freq: float,
                                min_depth: int,
                                min_baseq: int,
                                min_allele_count: int,
                                flank_codons: int = 2):
    cds_sequence = gene_record.get("cds_sequence", "").upper()
    if not cds_sequence:
        return []

    affected_indices = []
    for item in cluster:
        index = item["cds_index"]
        length_change = item["length_change"]
        affected_indices.extend(range(
            max(0, index - length_change),
            min(len(cds_sequence), index + length_change + 1),
        ))
    event_start = min(affected_indices)
    event_end = max(affected_indices)
    window_start = max(0, (event_start // 3) * 3 - flank_codons * 3)
    window_end = min(
        len(cds_sequence) - 1,
        (event_end // 3) * 3 + 2 + flank_codons * 3,
    )
    coding_positions = [
        cds_index_to_genomic_pos(gene_record, index)
        for index in range(window_start, window_end + 1)
    ]
    if any(position is None for position in coding_positions):
        return []

    sequence_counts = Counter()
    sequence_strands = defaultdict(Counter)
    sequence_events = defaultdict(Counter)
    for sam_line in sam_lines:
        result = _read_coding_haplotype(
            sam_line,
            coding_positions,
            gene_record.get("strand", "+"),
            min_baseq,
        )
        if result is None:
            continue
        sequence = result["sequence"]
        sequence_counts[sequence] += 1
        sequence_strands[sequence][result["strand"]] += 1
        sequence_events[sequence][result["events"]] += 1

    spanning_depth = sum(sequence_counts.values())
    if spanning_depth < min_depth:
        return [_no_call_row(
            gene,
            gene_record,
            cluster_id,
            cluster,
            window_start,
            window_end,
            coding_positions,
            spanning_depth,
            min_depth,
        )]

    qualifying = [
        (sequence, count)
        for sequence, count in sequence_counts.most_common()
        if count >= min_allele_count and count / spanning_depth >= min_alt_freq
    ]
    if not qualifying:
        qualifying = [sequence_counts.most_common(1)[0]]

    reference_local = cds_sequence[window_start:window_end + 1]
    if len(qualifying) >= 2:
        call_status = "MIXED_SIGNAL"
    elif qualifying[0][0] == reference_local:
        call_status = "REFERENCE"
    else:
        call_status = "VARIANT"

    spectrum = "; ".join(
        f"H{index}:{count} ({count / spanning_depth:.1%})"
        for index, (_, count) in enumerate(qualifying, start=1)
    )
    cluster_positions = [item["genomic_pos"] for item in cluster]
    cluster_events = "; ".join(
        f"{item['row'].get('pos')}:{item['row'].get('ref')}>{item['row'].get('alt')}"
        for item in cluster
    )
    rows = []
    for index, (sample_local, count) in enumerate(qualifying, start=1):
        full_sample_cds = (
            cds_sequence[:window_start]
            + sample_local
            + cds_sequence[window_end + 1:]
        )
        event_signature = sequence_events[sample_local].most_common(1)[0][0]
        consequence = _protein_consequence(
            cds_sequence,
            full_sample_cds,
            event_signature != "NO_INDEL",
        )
        strands = sequence_strands[sample_local]
        is_reference = sample_local == reference_local
        qc_flags = []
        if call_status == "MIXED_SIGNAL":
            qc_flags.append("MULTIPLE_HAPLOTYPES")
        if strands.get("forward", 0) == 0 or strands.get("reverse", 0) == 0:
            qc_flags.append("SINGLE_STRAND_SUPPORT")
        rows.append({
            "gene": gene,
            "gene_name": gene_record.get("description") or gene,
            "chrom": gene_record.get("resolved_chrom") or gene_record["chrom"],
            "cluster_id": cluster_id,
            "cluster_start": min(cluster_positions),
            "cluster_end": max(cluster_positions),
            "cluster_events": cluster_events,
            "haplotype_id": f"{cluster_id}_H{index}",
            "spanning_read_count": spanning_depth,
            "haplotype_read_count": count,
            "haplotype_frequency": f"{count / spanning_depth:.6f}",
            "forward_depth": strands.get("forward", 0),
            "reverse_depth": strands.get("reverse", 0),
            "haplotype_spectrum": spectrum,
            "call_status": call_status,
            "phase_status": "PHASED",
            "qc_flags": ";".join(qc_flags) or "PASS",
            "frame_status": (
                "REFERENCE" if is_reference else consequence["frame_status"]
            ),
            "effect": "reference" if is_reference else consequence["effect"],
            "net_nt_change": consequence["net_nt_change"],
            "cds_window_start": window_start + 1,
            "cds_window_end": window_end + 1,
            "genomic_window_start": min(coding_positions),
            "genomic_window_end": max(coding_positions),
            "reference_cds_haplotype": reference_local,
            "sample_cds_haplotype": sample_local,
            "raw_alignment_events": event_signature,
            "reference_peptide": consequence["reference_peptide"],
            "sample_peptide": consequence["sample_peptide"],
            "aa_start": consequence["aa_start"],
            "aa_end": consequence["aa_end"],
            "aa_changes": (
                "Reference haplotype" if is_reference else consequence["aa_changes"]
            ),
            "combined_aa_change": (
                "Reference haplotype"
                if is_reference else consequence["combined_aa_change"]
            ),
            "premature_stop": consequence["premature_stop"],
            "note": (
                "Protein consequence reconstructed from reads spanning the complete "
                "indel cluster; raw alignment events are retained separately."
            ),
        })
    return rows


def build_protein_haplotypes(gene: str,
                             gene_record: dict,
                             variant_rows: list[dict],
                             bam_path: Path,
                             min_mapq: int,
                             min_alt_freq: float,
                             min_depth: int,
                             min_baseq: int,
                             min_allele_count: int,
                             logger):
    rows = []
    chrom = gene_record.get("resolved_chrom") or gene_record["chrom"]
    for index, cluster in enumerate(_coding_indel_clusters(variant_rows), start=1):
        cluster_id = f"C{index}"
        cds_indices = [item["cds_index"] for item in cluster]
        maximum_length_change = max(
            item["length_change"] for item in cluster
        )
        event_start = max(
            0,
            min(cds_indices) - maximum_length_change - 9,
        )
        event_end = min(
            len(gene_record.get("cds_sequence", "")) - 1,
            max(cds_indices) + maximum_length_change + 9,
        )
        genomic_positions = [
            cds_index_to_genomic_pos(gene_record, cds_index)
            for cds_index in range(event_start, event_end + 1)
        ]
        genomic_positions = [pos for pos in genomic_positions if pos is not None]
        if not genomic_positions:
            continue
        sam_lines = _run_samtools_view(
            bam_path,
            chrom,
            min(genomic_positions),
            max(genomic_positions),
            min_mapq,
            logger,
        )
        rows.extend(_haplotype_rows_for_cluster(
            gene=gene,
            gene_record=gene_record,
            cluster_id=cluster_id,
            cluster=cluster,
            sam_lines=sam_lines,
            min_alt_freq=min_alt_freq,
            min_depth=min_depth,
            min_baseq=min_baseq,
            min_allele_count=min_allele_count,
        ))
    return rows


def write_protein_haplotypes(rows: list[dict], out_csv: Path):
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=HAPLOTYPE_FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    return out_csv
