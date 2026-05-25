import json
import csv
import argparse
import os


def normalize_strand(orientation):
    """
    NCBI orientation 可能是:
      - "plus" / "minus"
      - "+" / "-"
      - 1 / -1
      - True/False (少見)
    轉成 BED 的 strand: '+' / '-' / '.'
    """
    if orientation is None:
        return '.'

    # 數字型
    if isinstance(orientation, int):
        return '+' if orientation >= 0 else '-'

    # 字串型
    if isinstance(orientation, str):
        o = orientation.strip().lower()
        if o in ("plus", "+", "forward", "fwd", "1", "positive"):
            return '+'
        if o in ("minus", "-", "reverse", "rev", "-1", "negative"):
            return '-'

    return '.'


def parse_ncbi_gene_json(json_str):
    """解析單行 JSON 字串並提取：gene range + exons（for CSV）"""
    data = json.loads(json_str)

    if not data.get('transcripts'):
        return None

    transcript = data['transcripts'][0]
    location = transcript['genomicLocations'][0]
    genomic_range = location.get("genomicRange", {})

    symbol = data.get("symbol", "N/A")
    chrom = location.get("genomicAccessionVersion")
    orientation = genomic_range.get("orientation")
    strand = normalize_strand(orientation)

    if not location.get("exons"):
        return None

    sorted_exons = sorted(location['exons'], key=lambda x: x['order'])

    # exon min/max（保險）
    exon_min = min(int(ex["begin"]) for ex in sorted_exons)
    exon_max = max(int(ex["end"]) for ex in sorted_exons)

    gr_begin = genomic_range.get("begin")
    gr_end = genomic_range.get("end")

    gene_start_1based = exon_min if gr_begin is None else min(int(gr_begin), exon_min)
    gene_end_1based = exon_max if gr_end is None else max(int(gr_end), exon_max)

    # 給 CSV 用（1-based, inclusive）
    rows = []
    rows.append({
        "symbol": symbol,
        "chrom": chrom,
        "strand": strand,
        "type": "gene",
        "order": 0,
        "start": gene_start_1based,
        "end": gene_end_1based,
        "length": gene_end_1based - gene_start_1based + 1,
        "cds_offset": 0
    })

    accumulated_len = 0
    for ex in sorted_exons:
        begin = int(ex['begin'])
        end = int(ex['end'])
        length = end - begin + 1

        rows.append({
            "symbol": symbol,
            "chrom": chrom,
            "strand": strand,
            "type": "exon",
            "order": ex['order'],
            "start": begin,
            "end": end,
            "length": length,
            "cds_offset": accumulated_len
        })
        accumulated_len += length

    # 給 BED 用（0-based, half-open）
    # BED start = gene_start_1based - 1
    # BED end   = gene_end_1based （因為 half-open，不用 +1）
    bed = {
        "chrom": chrom,
        "start0": gene_start_1based - 1,
        "end0": gene_end_1based,
        "name": symbol,
        "score": 0,
        "strand": strand
    }

    return rows, bed


def main():
    parser = argparse.ArgumentParser(
        description="NCBI product_report.jsonl -> CSV (gene+exons) + BED (gene ranges for samtools -L)"
    )
    parser.add_argument("-i", "--input", required=True, help="Input NCBI product_report.jsonl file")
    parser.add_argument("--csv", required=True, help="Output CSV file path")
    parser.add_argument("--bed", required=True, help="Output BED file path (for samtools -L)")
    args = parser.parse_args()

    csv_fieldnames = ["symbol", "chrom", "strand", "type", "order", "start", "end", "length", "cds_offset"]

    try:
        with open(args.input, 'r', encoding='utf-8') as infile, \
             open(args.csv, 'w', newline='', encoding='utf-8') as out_csv, \
             open(args.bed, 'w', newline='', encoding='utf-8') as out_bed:

            writer = csv.DictWriter(out_csv, fieldnames=csv_fieldnames)
            writer.writeheader()

            # 避免重複 gene（同 symbol 多筆時只寫一次 bed），你也可拿掉這段
            seen_bed = set()

            for line in infile:
                if not line.strip():
                    continue

                parsed = parse_ncbi_gene_json(line)
                if not parsed:
                    continue

                rows, bed = parsed

                # CSV：gene + exons 全寫
                for r in rows:
                    writer.writerow(r)

                # BED：只寫 gene range
                bed_key = (bed["chrom"], bed["start0"], bed["end0"], bed["name"], bed["strand"])
                if bed_key not in seen_bed:
                    seen_bed.add(bed_key)
                    # BED6: chrom start end name score strand
                    out_bed.write(
                        f'{bed["chrom"]}\t{bed["start0"]}\t{bed["end0"]}\t{bed["name"]}\t{bed["score"]}\t{bed["strand"]}\n'
                    )

                print(f"Processed gene: {rows[0]['symbol']}")

        print(f"Done. We get the gene info.\nCSV saved: {args.csv}\nBED saved: {args.bed}")

    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()
