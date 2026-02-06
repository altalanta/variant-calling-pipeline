#!/usr/bin/env python3

"""Unit tests for entropy-aware variant helpers."""

import unittest
import sys
import types
from pathlib import Path

SCRIPTS_DIR = Path(__file__).parent.parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

# Provide a lightweight pysam stub so helper-only tests can import the module
if "pysam" not in sys.modules:
    pysam_stub = types.SimpleNamespace(
        AlignmentFile=object,
        FastaFile=object,
        VariantFile=object,
        tabix_index=lambda *args, **kwargs: None,
    )
    sys.modules["pysam"] = pysam_stub

from entropy_variant_pipeline import (  # noqa: E402
    DeBruijnGraph,
    IntervalIndex,
    VariantCandidate,
    shannon_entropy,
    variant_type,
)


class TestEntropyHelpers(unittest.TestCase):
    def test_shannon_entropy_balanced(self):
        """Balanced distribution should have 1 bit of entropy."""
        value = shannon_entropy([10, 10])
        self.assertAlmostEqual(value, 1.0, places=3)

    def test_shannon_entropy_empty(self):
        self.assertEqual(shannon_entropy([]), 0.0)

    def test_variant_type_detection(self):
        snv = VariantCandidate("chr1", 10, "A", "G", "baseline")
        mnv = VariantCandidate("chr1", 20, "AC", "GT", "assembly")
        ins = VariantCandidate("chr1", 30, "A", "AT", "assembly")
        micro_ins = VariantCandidate("chr1", 40, "A", "A" + "T" * 25, "assembly")
        del_var = VariantCandidate("chr1", 50, "AT", "A", "assembly")
        micro_del = VariantCandidate("chr1", 60, "A" + "C" * 25, "A", "assembly")
        self.assertEqual(variant_type(snv), "SNV")
        self.assertEqual(variant_type(mnv), "MNV")
        self.assertEqual(variant_type(ins), "insertion")
        self.assertEqual(variant_type(micro_ins), "micro_insertion")
        self.assertEqual(variant_type(del_var), "deletion")
        self.assertEqual(variant_type(micro_del), "micro_deletion")

    def test_debruijn_graph_simple_contig(self):
        graph = DeBruijnGraph(k=4)
        graph.add_sequence("AAGGTCC")
        graph.add_sequence("AAGGTGC")
        contigs = graph.assemble()
        self.assertTrue(any(seq.startswith("AAGGT") for seq, _ in contigs))

    def test_interval_index(self):
        idx = IntervalIndex()
        idx.add("chr1", 0, 100)
        idx.add("chr1", 200, 300)
        idx.finalize()
        self.assertTrue(idx.contains("chr1", 50))
        self.assertFalse(idx.contains("chr1", 150))
        self.assertTrue(idx.contains("chr1", 250))
        self.assertFalse(idx.contains("chr2", 10))


if __name__ == "__main__":
    unittest.main()
