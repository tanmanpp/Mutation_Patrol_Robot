# -*- coding: utf-8 -*-

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from modules.bam_utils import prepare_sorted_bam


class Logger:
    def info(self, _message):
        pass

    def error(self, _message):
        pass


class BamUtilsTests(unittest.TestCase):
    def test_uploaded_bam_is_sorted_into_canonical_location(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source = root / "uploaded sample.bam"
            destination = root / "result" / "work" / "bam" / "merged.sorted.bam"
            source.write_bytes(b"source")
            calls = []

            def fake_run(command, **_kwargs):
                calls.append(command)
                if command[:2] == ["samtools", "sort"]:
                    output = Path(command[command.index("-o") + 1])
                    output.write_bytes(b"sorted")
                elif command[:2] == ["samtools", "index"]:
                    Path(str(command[2]) + ".bai").write_bytes(b"index")
                return 0

            with patch("modules.bam_utils.run_cmd", side_effect=fake_run):
                result = prepare_sorted_bam(
                    source,
                    destination,
                    threads=4,
                    force=True,
                    logger=Logger(),
                )

            self.assertEqual(result, destination)
            self.assertEqual(destination.read_bytes(), b"sorted")
            self.assertTrue(Path(str(destination) + ".bai").exists())
            self.assertEqual(calls[0][:2], ["samtools", "sort"])
            self.assertEqual(calls[1][:2], ["samtools", "index"])


if __name__ == "__main__":
    unittest.main()
