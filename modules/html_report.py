# -*- coding: utf-8 -*-

import csv
import html
import json
from collections import Counter
from datetime import datetime
from pathlib import Path

from modules.utils import ensure_dir, relative_path


REPORT_NAME = "final_report.html"


def read_csv_records(path: Path):
    if not path.exists():
        return []
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def esc(value):
    return html.escape(str(value if value is not None else ""))


def number(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def table_html(rows: list[dict], empty: str, page_size: int = 20):
    if not rows:
        return f'<div class="empty">{esc(empty)}</div>'
    columns = list(rows[0].keys())
    head = "".join(f"<th>{esc(column)}</th>" for column in columns)
    body_rows = []
    for row in rows:
        cells = "".join(f"<td>{esc(row.get(column, ''))}</td>" for column in columns)
        body_rows.append(f"<tr>{cells}</tr>")
    return (
        f'<div class="paged-table" data-page-size="{page_size}">'
        '<div class="table-pager">'
        '<span class="page-status">Showing rows</span>'
        '<div class="pager-buttons">'
        '<button type="button" data-page-first>First</button>'
        '<button type="button" data-page-prev>Previous</button>'
        '<button type="button" data-page-next>Next</button>'
        '<button type="button" data-page-last>Last</button>'
        '</div>'
        '</div>'
        '<div class="table-wrap"><table>'
        f"<thead><tr>{head}</tr></thead>"
        f"<tbody>{''.join(body_rows)}</tbody>"
        "</table></div></div>"
    )


def metric_cards(metrics: dict):
    cards = []
    for label, value in metrics.items():
        cards.append(
            '<div class="metric">'
            f"<span>{esc(label)}</span>"
            f"<strong>{esc(value)}</strong>"
            "</div>"
        )
    return f'<div class="metrics">{"".join(cards)}</div>'


def bar_svg(items: list[tuple[str, float]], title: str, color: str = "#2f6f9f"):
    if not items:
        return '<div class="empty">No chart data.</div>'
    width = 920
    height = 260
    left = 56
    bottom = 58
    top = 28
    chart_w = width - left - 24
    chart_h = height - top - bottom
    max_value = max(value for _, value in items) or 1
    bar_gap = 8
    bar_w = max(8, (chart_w - bar_gap * (len(items) - 1)) / len(items))
    parts = [
        f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="{esc(title)}">',
        f'<text x="{left}" y="18" class="chart-title">{esc(title)}</text>',
        f'<line x1="{left}" y1="{top + chart_h}" x2="{width - 20}" y2="{top + chart_h}" class="axis"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + chart_h}" class="axis"/>',
    ]
    for index, (label, value) in enumerate(items):
        x = left + index * (bar_w + bar_gap)
        h = (value / max_value) * chart_h
        y = top + chart_h - h
        parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{h:.1f}" fill="{color}"/>')
        parts.append(f'<text x="{x + bar_w / 2:.1f}" y="{y - 4:.1f}" text-anchor="middle" class="value">{value:g}</text>')
        parts.append(
            f'<text x="{x + bar_w / 2:.1f}" y="{top + chart_h + 16}" text-anchor="end" '
            f'transform="rotate(-35 {x + bar_w / 2:.1f},{top + chart_h + 16})" class="label">{esc(label)}</text>'
        )
    parts.append("</svg>")
    return "".join(parts)


def site_query_chart(rows: list[dict], title: str):
    if not rows:
        return '<div class="empty">No site query chart data.</div>'
    items = []
    for row in rows:
        gene_name = row.get("Gene name") or row.get("Gene symbol") or row.get("gene") or ""
        position = row.get("aa_pos") or row.get("pos") or ""
        allele = row.get("alt_codon") or row.get("alt") or ""
        label = f"{gene_name} {position} {allele}".strip()
        depth = number(row.get("depth"))
        freq = number(row.get("allele_freq"))
        items.append((label, depth, freq, row.get("variant_type", "")))

    max_depth = max(depth for _, depth, _, _ in items) or 1
    width = 920
    height = 300
    left = 56
    bottom = 72
    top = 30
    chart_w = width - left - 42
    chart_h = height - top - bottom
    gap = 8
    bar_w = max(8, (chart_w - gap * (len(items) - 1)) / len(items))
    parts = [
        f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="{esc(title)}">',
        f'<text x="{left}" y="18" class="chart-title">{esc(title)}</text>',
        f'<line x1="{left}" y1="{top + chart_h}" x2="{width - 20}" y2="{top + chart_h}" class="axis"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + chart_h}" class="axis"/>',
    ]
    points = []
    for index, (label, depth, freq, variant_type) in enumerate(items):
        x = left + index * (bar_w + gap)
        h = (depth / max_depth) * chart_h
        y = top + chart_h - h
        color = "#8fa1a6" if depth == 0 else ("#79a99b" if variant_type == "REF" else "#b85c38")
        parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{h:.1f}" fill="{color}"/>')
        point_x = x + bar_w / 2
        point_y = top + chart_h - (freq * chart_h)
        points.append((point_x, point_y))
        parts.append(
            f'<text x="{point_x:.1f}" y="{top + chart_h + 16}" text-anchor="end" '
            f'transform="rotate(-35 {point_x:.1f},{top + chart_h + 16})" class="label">{esc(label)}</text>'
        )
    if len(points) > 1:
        path = " ".join(("M" if idx == 0 else "L") + f"{x:.1f},{y:.1f}" for idx, (x, y) in enumerate(points))
        parts.append(f'<path d="{path}" fill="none" stroke="#2f6f9f" stroke-width="2"/>')
    for x, y in points:
        parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="#2f6f9f"/>')
    parts.append("</svg>")
    return "".join(parts)


def load_summary(result_dir: Path):
    summary_path = result_dir / "tables" / "summary.json"
    if not summary_path.exists():
        return {}
    return json.loads(summary_path.read_text(encoding="utf-8"))


def load_site_queries(result_dir: Path):
    site_root = result_dir / "site_query"
    if not site_root.exists():
        return []
    queries = []
    for query_dir in sorted([path for path in site_root.iterdir() if path.is_dir()], key=lambda path: path.stat().st_mtime):
        csv_path = query_dir / "site_query.csv"
        if not csv_path.exists():
            continue
        rows = read_csv_records(csv_path)
        gene_name = rows[0].get("Gene name", query_dir.name) if rows else query_dir.name
        queries.append({
            "query_id": query_dir.name,
            "gene_name": gene_name,
            "path": csv_path,
            "rows": rows,
        })
    return queries


def render_html_report(result_dir: Path, out_html: Path | None = None):
    result_dir = Path(result_dir)
    if out_html is None:
        out_html = result_dir / REPORT_NAME
    ensure_dir(out_html.parent)

    tables = result_dir / "tables"
    sample_summary = read_csv_records(tables / "sample_summary.csv")
    mutation_candidates = read_csv_records(tables / "mutation_candidates.csv")
    summary = load_summary(result_dir)
    site_queries = load_site_queries(result_dir)

    variant_counts = Counter(row.get("variant_type", "NA") or "NA" for row in mutation_candidates)
    site_query_count = sum(len(item["rows"]) for item in site_queries)
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    sections = [
        "<!doctype html><html><head><meta charset=\"utf-8\">",
        f"<title>Mutation Patrol Robot Report - {esc(result_dir.name)}</title>",
        STYLE,
        "</head><body>",
        "<main>",
        f"<h1>Mutation Patrol Robot Report</h1>",
        f"<p class=\"subtle\">Sample/output folder: <strong>{esc(result_dir.name)}</strong></p>",
        f"<p class=\"subtle\">Generated at: {esc(generated_at)}</p>",
        metric_cards({
            "Sample summary rows": len(sample_summary),
            "Mutation candidates": len(mutation_candidates),
            "Site query tables": len(site_queries),
            "Site query rows": site_query_count,
            "Output folder": relative_path(result_dir),
        }),
        "<h2>Analysis Summary</h2>",
        table_html(sample_summary, "No sample summary found."),
        "<h2>Mutation Candidate Overview</h2>",
        bar_svg(list(variant_counts.items()), "Mutation candidates by variant type", "#b85c38"),
        "<h2>Mutation Candidates</h2>",
        table_html(mutation_candidates, "No mutation candidates found."),
        "<h2>Site Query</h2>",
    ]

    if site_queries:
        for item in site_queries:
            sections.extend([
                f"<section class=\"site-query\"><h3>{esc(item['gene_name'])}</h3>",
                f"<p class=\"subtle\">Source: {esc(relative_path(item['path']))}</p>",
                site_query_chart(item["rows"], f"Depth and allele frequency - {item['gene_name']}"),
                table_html(item["rows"], "No site query rows found."),
                "</section>",
            ])
    else:
        sections.append('<div class="empty">No site query outputs found.</div>')

    if summary:
        sections.extend([
            "<h2>Run Metadata</h2>",
            f"<pre>{esc(json.dumps(summary, indent=2, ensure_ascii=False))}</pre>",
        ])

    sections.extend(["</main></body></html>"])
    sections.insert(-1, SCRIPT)
    out_html.write_text("\n".join(sections), encoding="utf-8")
    return out_html


STYLE = """
<style>
:root{color:#1d2a31;font-family:Inter,Arial,sans-serif;background:#eef3f1}
body{margin:0}
main{max-width:1180px;margin:0 auto;padding:28px}
h1{margin:0 0 8px;font-size:30px}
h2{margin:32px 0 12px;font-size:22px}
h3{margin:22px 0 8px;font-size:18px}
.subtle{color:#586c72}
.metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:10px;margin:20px 0}
.metric{display:grid;gap:6px;padding:12px;border:1px solid #d1ded9;border-radius:8px;background:white}
.metric span{font-size:12px;color:#586c72;font-weight:700}
.metric strong{font-size:20px;color:#163238;word-break:break-word}
.table-wrap{overflow:auto;border:1px solid #d1ded9;border-radius:8px;background:white}
.paged-table{display:grid;gap:8px}
.table-pager{display:flex;align-items:center;justify-content:space-between;gap:12px;padding:8px 10px;border:1px solid #d1ded9;border-radius:8px;background:white;color:#40545a;font-size:12px;font-weight:700}
.pager-buttons{display:flex;gap:6px}
.pager-buttons button{min-height:28px;border:1px solid #c8d6d2;border-radius:6px;background:#fff;color:#1d2a31;cursor:pointer}
.pager-buttons button:disabled{cursor:not-allowed;opacity:.45}
table{width:100%;border-collapse:collapse;font-size:12px}
th,td{padding:8px 10px;border-bottom:1px solid #e3ece8;text-align:left;white-space:nowrap}
th{position:sticky;top:0;background:#f4faf7;color:#2f444a}
.empty{padding:16px;border:1px dashed #b7c8c3;border-radius:8px;background:white;color:#586c72}
svg{width:100%;height:auto;margin:8px 0 14px;border:1px solid #d1ded9;border-radius:8px;background:white}
.axis{stroke:#8fa1a6;stroke-width:1}
.chart-title{font-size:14px;font-weight:700;fill:#163238}
.label{font-size:10px;fill:#40545a}
.value{font-size:10px;fill:#40545a}
pre{overflow:auto;padding:14px;border:1px solid #d1ded9;border-radius:8px;background:white;font-size:12px}
.site-query{margin-top:12px}
</style>
"""


SCRIPT = """
<script>
document.querySelectorAll(".paged-table").forEach(function(container) {
  const pageSize = Number(container.dataset.pageSize || 20);
  const rows = Array.from(container.querySelectorAll("tbody tr"));
  const status = container.querySelector(".page-status");
  const first = container.querySelector("[data-page-first]");
  const prev = container.querySelector("[data-page-prev]");
  const next = container.querySelector("[data-page-next]");
  const last = container.querySelector("[data-page-last]");
  const totalPages = Math.max(1, Math.ceil(rows.length / pageSize));
  let page = 1;

  function render() {
    const start = (page - 1) * pageSize;
    const end = Math.min(start + pageSize, rows.length);
    rows.forEach(function(row, index) {
      row.style.display = index >= start && index < end ? "" : "none";
    });
    status.textContent = rows.length
      ? `Showing ${start + 1}-${end} of ${rows.length} rows (page ${page}/${totalPages})`
      : "No rows";
    first.disabled = page === 1;
    prev.disabled = page === 1;
    next.disabled = page === totalPages;
    last.disabled = page === totalPages;
  }

  first.addEventListener("click", function() { page = 1; render(); });
  prev.addEventListener("click", function() { page = Math.max(1, page - 1); render(); });
  next.addEventListener("click", function() { page = Math.min(totalPages, page + 1); render(); });
  last.addEventListener("click", function() { page = totalPages; render(); });
  render();
});
</script>
"""
