# -*- coding: utf-8 -*-

import time
import unittest
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

from backend import jobs, services
from modules.bam_allele_freq import _count_pileup_events
from modules.coverage import (
    _coverage_stats,
    _run_samtools_depth,
    generate_coverage_tables,
)
from modules.gene_database import (
    cds_index_to_genomic_pos,
    genomic_pos_to_cds_index,
    load_gene_database,
)
from modules.site_query import (
    _apply_protein_haplotype_annotations,
    _biological_position_sort_key,
    _codon_rows_for_aa_queries,
    _compress_no_call_regions,
    _rows_for_site,
    _site_call_status,
    _site_filter,
    scan_bam_gene_region,
)
from modules.utils import run_bash_pipeline


PROJECT_ROOT = Path(__file__).resolve().parents[1]


class CoreTests(unittest.TestCase):
    def test_phased_haplotype_supersedes_provisional_site_effects(self):
        variants = [
            {
                "region_type": "CDS",
                "variant_type": "INS",
                "pos": "403618",
                "aa_pos": "74",
                "effect": "frameshift",
                "aa_change": "M74fs",
                "qc_flags": "PASS",
                "note": "raw insertion",
            },
            {
                "region_type": "CDS",
                "variant_type": "SNV",
                "pos": "403625",
                "aa_pos": "76",
                "effect": "missense",
                "aa_change": "K76T",
                "qc_flags": "PASS",
                "note": "raw substitution",
            },
        ]
        haplotypes = [{
            "cluster_id": "C1",
            "cluster_start": "403618",
            "cluster_end": "403622",
            "phase_status": "PHASED",
            "call_status": "VARIANT",
            "haplotype_frequency": "0.987667",
            "aa_start": "74",
            "aa_end": "76",
            "effect": "frame_restored_complex",
            "aa_changes": "M74I; N75E; K76T",
            "combined_aa_change": "M74_K76delinsIET",
        }]

        _apply_protein_haplotype_annotations(variants, haplotypes)

        self.assertTrue(all(
            row["aa_change"] == "M74I; N75E; K76T"
            for row in variants
        ))
        self.assertTrue(all(
            "HAPLOTYPE_RECONSTRUCTED" in row["qc_flags"]
            for row in variants
        ))
        self.assertIn("site_level_aa_change=M74fs", variants[0]["note"])

    def test_igv_config_opens_whole_gene_reference_and_alignment(self):
        with TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            result_dir = root / "results" / "run1"
            query_dir = result_dir / "site_query" / "gene1_whole_gene"
            query_dir.mkdir(parents=True)
            (result_dir / "tables").mkdir()
            bam_path = result_dir / "work" / "bam" / "merged.sorted.bam"
            bam_path.parent.mkdir(parents=True)
            bam_path.write_bytes(b"BAM")
            bam_index = Path(str(bam_path) + ".bai")
            bam_index.write_bytes(b"BAI")

            reference_path = root / "reference.fa"
            reference_path.write_text(
                ">chr1\n" + ("A" * 300) + "\n",
                encoding="utf-8",
            )
            reference_index = Path(str(reference_path) + ".fai")
            reference_index.write_text(
                "chr1\t300\t6\t300\t301\n",
                encoding="utf-8",
            )
            database_path = root / "gene_database.json"
            database_path.write_text(json.dumps({
                "reference_fasta": str(reference_path),
                "genes": [{
                    "symbol": "gene1",
                    "chrom": "chr1",
                    "resolved_chrom": "chr1",
                    "strand": "-",
                    "gene_start": 100,
                    "gene_end": 200,
                }],
            }), encoding="utf-8")
            (query_dir / "igv_context.json").write_text(json.dumps({
                "query_type": "gene_region",
                "database": "test_database",
                "gene_database": str(database_path),
                "reference_fasta": str(reference_path),
                "bam": str(bam_path),
                "gene": "gene1",
                "chrom": "chr1",
                "start": 100,
                "end": 200,
            }), encoding="utf-8")

            with patch.object(services, "RESULTS_DIR", root / "results"):
                config = services.get_igv_config(
                    "run1",
                    "gene1_whole_gene",
                )
                alignment_path, media_type = (
                    services.get_igv_resource_path(
                        "run1",
                        "gene1_whole_gene",
                        "alignment",
                    )
                )

            self.assertEqual(config["locus"], "chr1:50-250")
            self.assertEqual(config["gene"], "gene1")
            self.assertEqual(config["strand"], "-")
            self.assertTrue(
                config["reference"]["fasta_url"].endswith(
                    "/igv/reference"
                )
            )
            self.assertTrue(
                config["alignment"]["index_url"].endswith(
                    "/igv/alignment-index"
                )
            )
            self.assertEqual(alignment_path, bam_path)
            self.assertEqual(media_type, "application/octet-stream")

    def test_amino_acid_rows_sort_from_n_to_c_terminus(self):
        rows = [
            {"gene": "minus", "aa_pos": "10", "cds_pos": "30", "pos": "100"},
            {"gene": "minus", "aa_pos": "2", "cds_pos": "6", "pos": "200"},
            {"gene": "minus", "aa_pos": "2", "cds_pos": "4", "pos": "202"},
        ]

        ordered = sorted(rows, key=_biological_position_sort_key)

        self.assertEqual(
            [(row["aa_pos"], row["cds_pos"]) for row in ordered],
            [("2", "4"), ("2", "6"), ("10", "30")],
        )

    def test_safe_names_cannot_escape_storage(self):
        run_name = services.safe_name("../../sample 01")
        self.assertNotIn("/", run_name)
        self.assertNotIn("\\", run_name)
        self.assertNotIn(run_name, {".", ".."})
        self.assertEqual(
            services.safe_upload_filename(r"..\folder\reads 01.fastq.gz"),
            "reads_01.fastq.gz",
        )

    def test_query_positions_are_validated(self):
        self.assertEqual(services.parse_query_values("1, 2,10"), [1, 2, 10])
        for invalid in ("", "0", "-1", "one"):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    services.parse_query_values(invalid)

    def test_pileup_parser_counts_bases_and_indels(self):
        events = _count_pileup_events(".,Aa+2tt-1c^F.$", "G")
        self.assertEqual(events["bases"], {"A": 2, "C": 0, "G": 3, "T": 0})
        self.assertEqual(events["insertions"], {"TT": 1})
        self.assertEqual(events["deletions"], {"C": 1})

    def test_fixture_gene_coordinates_round_trip(self):
        db_path = PROJECT_ROOT / "test_data" / "3D7_database" / "gene_database.json"
        database = load_gene_database(db_path)
        self.assertEqual(database["gene_count"], 6)

        for gene in database["genes"]:
            chrom = gene.get("resolved_chrom") or gene["chrom"]
            cds_indexes = (0, int(gene["cds_length"]) // 2, int(gene["cds_length"]) - 1)
            for cds_index in cds_indexes:
                genomic_pos = cds_index_to_genomic_pos(gene, cds_index)
                self.assertEqual(
                    genomic_pos_to_cds_index(gene, chrom, genomic_pos),
                    cds_index,
                )

    def test_background_job_reports_success_and_failure(self):
        successful = jobs.submit_job("test", "success", lambda: {"ok": True})
        failed = jobs.submit_job("test", "failure", lambda: 1 / 0)

        deadline = time.time() + 3
        while time.time() < deadline:
            successful = jobs.get_job(successful["job_id"])
            failed = jobs.get_job(failed["job_id"])
            if successful["status"] == "completed" and failed["status"] == "failed":
                break
            time.sleep(0.02)

        self.assertEqual(successful["status"], "completed")
        self.assertEqual(successful["result"], {"ok": True})
        self.assertEqual(failed["status"], "failed")
        self.assertEqual(failed["error_type"], "ZeroDivisionError")

    def test_pipeline_arguments_are_shell_quoted(self):
        with patch("modules.utils.run_cmd", return_value=0) as mocked:
            run_bash_pipeline(
                [["tool", "sample name;touch unsafe", "reference.fa"]],
                logger=object(),
            )

        shell_command = mocked.call_args.args[0]
        self.assertEqual(shell_command[:2], ["bash", "-lc"])
        self.assertIn("'sample name;touch unsafe'", shell_command[2])

    def test_site_call_requires_minimum_depth(self):
        self.assertEqual(_site_filter(0, 10), "NO_MAPPING")
        self.assertEqual(_site_filter(9, 10), "LOW_DEPTH")
        self.assertEqual(_site_filter(10, 10), "PASS")
        self.assertEqual(
            _site_call_status(9, 10, {"A": 9}, "A", 2, 0.05),
            "NO_CALL",
        )
        self.assertEqual(
            _site_call_status(10, 10, {"A": 10}, "A", 2, 0.05),
            "REFERENCE",
        )
        self.assertEqual(
            _site_call_status(20, 10, {"A": 14, "G": 6}, "A", 2, 0.05),
            "MIXED_SIGNAL",
        )
        self.assertEqual(
            _site_call_status(20, 10, {"A": 0, "G": 20}, "A", 2, 0.05),
            "VARIANT",
        )

    def test_uncovered_codon_is_not_reported_as_reference(self):
        rows = _codon_rows_for_aa_queries(
            gene="example",
            gene_record={"description": "Example gene", "strand": "+"},
            codon_queries=[{
                "chrom": "chr1",
                "aa_pos": 1,
                "cds_start": 1,
                "genomic_positions": [1, 2, 3],
                "ref_codon": "ATG",
                "ref_aa": "M",
            }],
            sam_lines=[],
            requested_alt=None,
            min_alt_freq=0.05,
            min_depth=10,
            min_mapq=20,
            min_baseq=20,
            min_allele_count=2,
        )

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["call_status"], "NO_CALL")
        self.assertEqual(rows[0]["variant_type"], "NO_CALL")
        self.assertEqual(rows[0]["filter"], "NO_MAPPING")
        self.assertEqual(rows[0]["allele_status"], "NO_CALL")

    def test_coverage_stats_distinguish_partial_and_no_call_bases(self):
        stats = _coverage_stats([0, 5, 10, 20], minimum_call_depth=10)
        self.assertEqual(stats["length"], 4)
        self.assertEqual(stats["covered_bases"], 3)
        self.assertEqual(stats["callable_bases"], 2)
        self.assertEqual(stats["no_call_bases"], 2)
        self.assertEqual(stats["pct_callable"], "0.500000")
        self.assertEqual(stats["qc_status"], "PARTIAL")

    def test_mixed_signal_reports_allele_spectrum_and_strands(self):
        with patch(
            "modules.site_query.annotate_candidate_row",
            side_effect=lambda row, _database: row,
        ):
            rows = _rows_for_site(
                gene="example",
                gene_record={"description": "Example"},
                gene_db={},
                query={
                    "chrom": "chr1",
                    "pos": 10,
                    "query_type": "genomic_pos",
                    "query_label": "10",
                },
                pileup_parts=["chr1", "10", "A", "20", "." * 14 + "g" * 6],
                requested_alt=None,
                min_alt_freq=0.05,
                min_depth=10,
                min_mapq=20,
                min_baseq=20,
                min_allele_count=2,
            )

        self.assertEqual({row["alt"] for row in rows}, {"A", "G"})
        self.assertTrue(all(row["call_status"] == "MIXED_SIGNAL" for row in rows))
        self.assertTrue(all(row["qc_flags"] == "MULTIPLE_ALLELES" for row in rows))
        self.assertEqual(rows[0]["allele_spectrum"], "A:14 (70.0%); G:6 (30.0%)")
        g_row = next(row for row in rows if row["alt"] == "G")
        self.assertEqual(g_row["forward_depth"], 0)
        self.assertEqual(g_row["reverse_depth"], 6)

    def test_coverage_command_uses_base_and_mapping_quality(self):
        completed = SimpleNamespace(stdout="chr1\t2\t7\nchr1\t3\t12\n")
        with patch("modules.coverage.run_cmd", return_value=completed) as mocked:
            depths = _run_samtools_depth(
                bam_path=Path("sample.bam"),
                chrom="chr1",
                start=1,
                end=3,
                min_mapq=25,
                min_baseq=18,
                logger=object(),
            )

        self.assertEqual(depths, [0, 7, 12])
        command = mocked.call_args.args[0]
        self.assertEqual(command[:2], ["samtools", "depth"])
        self.assertIn(["-q", "18"], [command[index:index + 2] for index in range(len(command) - 1)])
        self.assertIn(["-Q", "25"], [command[index:index + 2] for index in range(len(command) - 1)])

    def test_coverage_tables_include_gene_region_and_positions(self):
        with TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            database_path = root / "gene_database.json"
            database_path.write_text(json.dumps({
                "reference_fasta": "reference.fa",
                "genes": [{
                    "symbol": "gene/one",
                    "description": "Example gene",
                    "chrom": "chr1",
                    "resolved_chrom": "chr1",
                    "gene_start": 1,
                    "gene_end": 3,
                    "exons": [{
                        "order": 1,
                        "start": 1,
                        "end": 3,
                        "length": 3,
                        "cds_offset": 0,
                    }],
                }],
            }), encoding="utf-8")
            bam_path = root / "sample.bam"
            bam_path.touch()
            out_dir = root / "coverage"

            with (
                patch("modules.coverage.ensure_bam_index"),
                patch(
                    "modules.coverage._run_samtools_depth",
                    return_value=[0, 10, 20],
                ),
            ):
                generate_coverage_tables(
                    gene_db_path=database_path,
                    bam_path=bam_path,
                    out_dir=out_dir,
                    minimum_call_depth=10,
                    min_mapq=20,
                    min_baseq=20,
                    force=True,
                    logger=MagicMock(),
                )

            self.assertTrue((out_dir / "gene_coverage.csv").exists())
            self.assertTrue((out_dir / "region_coverage.csv").exists())
            index = json.loads(
                (out_dir / "coverage_index.json").read_text(encoding="utf-8")
            )
            position_file = index["genes"]["gene/one"]["file"]
            self.assertNotIn("/", position_file)
            position_rows = services.read_csv_records(
                out_dir / "positions" / position_file
            )
            self.assertEqual(
                [row["call_status"] for row in position_rows],
                ["NO_CALL", "CALLABLE", "CALLABLE"],
            )

    def test_no_call_positions_are_compressed_by_reason(self):
        regions = _compress_no_call_regions(
            "gene1",
            "chr1",
            [
                (10, 0, "NO_MAPPING"),
                (11, 0, "NO_MAPPING"),
                (12, 3, "LOW_DEPTH"),
                (13, 4, "LOW_DEPTH"),
            ],
            min_depth=10,
        )
        self.assertEqual(len(regions), 2)
        self.assertEqual((regions[0]["start"], regions[0]["end"]), (10, 11))
        self.assertEqual(regions[1]["reason"], "LOW_DEPTH")

    def test_whole_gene_scan_reports_variants_and_no_call_summary(self):
        with TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            database_path = root / "gene_database.json"
            reference_path = root / "reference.fa"
            reference_path.write_text(">chr1\nAAAA\n", encoding="utf-8")
            database_path.write_text(json.dumps({
                "reference_fasta": str(reference_path),
                "genes": [{
                    "symbol": "gene1",
                    "description": "Example gene",
                    "chrom": "chr1",
                    "resolved_chrom": "chr1",
                    "strand": "+",
                    "gene_start": 1,
                    "gene_end": 4,
                    "cds_sequence": "AAAA",
                    "exons": [{
                        "order": 1,
                        "start": 1,
                        "end": 4,
                        "length": 4,
                        "cds_offset": 0,
                    }],
                }],
            }), encoding="utf-8")
            bam_path = root / "sample.bam"
            bam_path.touch()
            out_csv = root / "scan" / "site_query.csv"
            pileup = [
                "chr1\t1\tA\t10\t..........",
                "chr1\t2\tA\t10\t.....ggggg",
                "chr1\t3\tA\t1\t.",
                "chr1\t4\tA\t0\t*",
            ]

            with (
                patch("modules.site_query.ensure_bam_index"),
                patch("modules.site_query._run_mpileup", return_value=pileup),
            ):
                outputs = scan_bam_gene_region(
                    gene_db_path=database_path,
                    bam_path=bam_path,
                    out_csv=out_csv,
                    gene="gene1",
                    ref_fasta=reference_path,
                    min_mapq=20,
                    min_alt_freq=0.05,
                    min_depth=10,
                    min_baseq=20,
                    min_allele_count=2,
                    force=True,
                    logger=MagicMock(),
                )

            variant_rows = services.read_csv_records(outputs["table"])
            self.assertEqual(len(variant_rows), 1)
            self.assertEqual(variant_rows[0]["pos"], "2")
            self.assertEqual(variant_rows[0]["alt"], "G")
            self.assertEqual(variant_rows[0]["call_status"], "MIXED_SIGNAL")
            self.assertEqual(variant_rows[0]["region_type"], "CDS")
            self.assertEqual(variant_rows[0]["effect"], "missense")
            self.assertEqual(variant_rows[0]["ref_aa_name"], "Lysine")
            self.assertEqual(variant_rows[0]["alt_aa_name"], "Arginine")
            summary = json.loads(
                outputs["summary"].read_text(encoding="utf-8")
            )
            self.assertEqual(summary["variant_sites"], 1)
            self.assertEqual(summary["mixed_signal_sites"], 1)
            self.assertEqual(summary["no_call_positions"], 2)
            self.assertEqual(summary["callable_positions"], 2)
            complete_rows = services.read_csv_records(
                outputs["complete_table"]
            )
            self.assertEqual(len(complete_rows), 5)
            self.assertEqual(summary["complete_table_rows"], 5)
            self.assertEqual(summary["protein_haplotype_rows"], 0)
            self.assertTrue(outputs["protein_haplotypes"].exists())
            self.assertEqual(
                {row["call_status"] for row in complete_rows},
                {"REFERENCE", "MIXED_SIGNAL", "NO_CALL"},
            )


if __name__ == "__main__":
    unittest.main()
