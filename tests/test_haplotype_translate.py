import unittest

from modules.haplotype_translate import (
    _coding_indel_clusters,
    _haplotype_rows_for_cluster,
    _protein_consequence,
    _read_coding_haplotype,
)


class ProteinHaplotypeTests(unittest.TestCase):
    def setUp(self):
        prefix = "GCT" * 71
        crt_window = "TGTGTAATGAATAAAATTTTTGCTAAA"
        self.gene_record = {
            "symbol": "PF3D7_0709000",
            "description": "chloroquine resistance transporter",
            "chrom": "chr7",
            "resolved_chrom": "chr7",
            "strand": "+",
            "gene_start": 403399,
            "gene_end": 403399 + len(prefix + crt_window) - 1,
            "cds_sequence": prefix + crt_window,
            "exons": [{
                "order": 1,
                "start": 403399,
                "end": 403399 + len(prefix + crt_window) - 1,
                "length": len(prefix + crt_window),
                "cds_offset": 0,
            }],
        }
        self.variant_rows = [
            {
                "variant_type": "INS",
                "region_type": "CDS",
                "cds_pos": "220",
                "pos": "403618",
                "ref": "A",
                "alt": "AT",
            },
            {
                "variant_type": "DEL",
                "region_type": "CDS",
                "cds_pos": "224",
                "pos": "403622",
                "ref": "AT",
                "alt": "A",
            },
        ]

    @staticmethod
    def _sam_line(name, flag, start, cigar, sequence):
        return "\t".join([
            name,
            str(flag),
            "chr7",
            str(start),
            "60",
            cigar,
            "*",
            "0",
            "0",
            sequence,
            "I" * len(sequence),
        ])

    def test_negative_strand_read_is_reconstructed_in_coding_orientation(self):
        sam_line = self._sam_line("negative-gene", 0, 101, "3M", "CAT")
        result = _read_coding_haplotype(
            sam_line,
            coding_positions=[103, 102, 101],
            strand="-",
            min_baseq=20,
        )

        self.assertEqual(result["sequence"], "ATG")

    def test_compensating_indels_translate_crt_haplotype(self):
        cluster = _coding_indel_clusters(self.variant_rows)[0]
        reference_local = self.gene_record["cds_sequence"][210:231]

        sample_local = reference_local[:10] + "T" + reference_local[10:]
        sample_local = sample_local[:15] + sample_local[16:]
        sample_local = sample_local[:16] + "C" + sample_local[17:]

        sam_lines = [
            self._sam_line(
                f"variant-{index}",
                0 if index % 2 == 0 else 16,
                403609,
                "10M1I4M1D6M",
                sample_local,
            )
            for index in range(10)
        ]
        sam_lines.append(self._sam_line(
            "reference-read",
            0,
            403609,
            "21M",
            reference_local,
        ))

        rows = _haplotype_rows_for_cluster(
            gene="PF3D7_0709000",
            gene_record=self.gene_record,
            cluster_id="C1",
            cluster=cluster,
            sam_lines=sam_lines,
            min_alt_freq=0.05,
            min_depth=10,
            min_baseq=20,
            min_allele_count=2,
        )

        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row["call_status"], "VARIANT")
        self.assertEqual(row["phase_status"], "PHASED")
        self.assertEqual(row["frame_status"], "FRAME_RESTORED")
        self.assertEqual(row["aa_changes"], "M74I; N75E; K76T")
        self.assertEqual(row["combined_aa_change"], "M74_K76delinsIET")
        self.assertIn("CVMNKIFAK", row["reference_peptide"])
        self.assertIn("CVIETIFAK", row["sample_peptide"])
        self.assertEqual(row["raw_alignment_events"], "403618:+T;403622:-1")
        self.assertEqual(row["forward_depth"], 5)
        self.assertEqual(row["reverse_depth"], 5)

    def test_unspanned_indel_cluster_is_no_call(self):
        cluster = _coding_indel_clusters(self.variant_rows)[0]
        rows = _haplotype_rows_for_cluster(
            gene="PF3D7_0709000",
            gene_record=self.gene_record,
            cluster_id="C1",
            cluster=cluster,
            sam_lines=[],
            min_alt_freq=0.05,
            min_depth=10,
            min_baseq=20,
            min_allele_count=2,
        )

        self.assertEqual(rows[0]["call_status"], "NO_CALL")
        self.assertEqual(rows[0]["phase_status"], "PHASE_UNRESOLVED")
        self.assertEqual(rows[0]["aa_changes"], "Unable to determine")

    def test_multiple_protein_haplotypes_are_mixed_signal(self):
        cluster = _coding_indel_clusters(self.variant_rows)[0]
        reference_local = self.gene_record["cds_sequence"][210:231]
        sample_local = reference_local[:10] + "T" + reference_local[10:]
        sample_local = sample_local[:15] + sample_local[16:]
        sample_local = sample_local[:16] + "C" + sample_local[17:]
        sam_lines = [
            self._sam_line(
                f"variant-{index}",
                index % 2 * 16,
                403609,
                "10M1I4M1D6M",
                sample_local,
            )
            for index in range(6)
        ] + [
            self._sam_line(
                f"reference-{index}",
                index % 2 * 16,
                403609,
                "21M",
                reference_local,
            )
            for index in range(4)
        ]

        rows = _haplotype_rows_for_cluster(
            gene="PF3D7_0709000",
            gene_record=self.gene_record,
            cluster_id="C1",
            cluster=cluster,
            sam_lines=sam_lines,
            min_alt_freq=0.05,
            min_depth=10,
            min_baseq=20,
            min_allele_count=2,
        )

        self.assertEqual(len(rows), 2)
        self.assertTrue(all(row["call_status"] == "MIXED_SIGNAL" for row in rows))
        self.assertEqual(rows[0]["haplotype_frequency"], "0.600000")
        self.assertEqual(rows[1]["haplotype_frequency"], "0.400000")
        self.assertTrue(all("MULTIPLE_HAPLOTYPES" in row["qc_flags"] for row in rows))

    def test_uncompensated_insertion_remains_frameshift(self):
        reference = "ATGAAAACCTTTGCTGACTAA"
        sample = reference[:3] + "A" + reference[3:]
        consequence = _protein_consequence(
            reference,
            sample,
            has_alignment_indel=True,
        )

        self.assertEqual(consequence["frame_status"], "PERSISTENT_FRAMESHIFT")
        self.assertEqual(consequence["effect"], "frameshift")
        self.assertIn("fs", consequence["combined_aa_change"])


if __name__ == "__main__":
    unittest.main()
