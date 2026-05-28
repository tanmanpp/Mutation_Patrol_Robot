#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

from modules.html_report import render_html_report


def build_parser():
    parser = argparse.ArgumentParser(
        description="Export one Mutation Patrol Robot output folder as a final HTML report.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--out_dir", required=True, help="Analysis output folder, e.g. app_data/results/sample_01")
    parser.add_argument("--html", default=None, help="Optional output HTML path. Defaults to <out_dir>/final_report.html")
    return parser


def main():
    args = build_parser().parse_args()
    out_dir = Path(args.out_dir)
    html_path = Path(args.html) if args.html else None
    report_path = render_html_report(out_dir, html_path)
    print(report_path)


if __name__ == "__main__":
    main()
