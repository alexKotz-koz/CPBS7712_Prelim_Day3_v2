import os
import sys

sys.path.insert(0, "../src/components")
import unittest
from collections import deque
from createContigs import CreateContigs


class TestCreateContigs(unittest.TestCase):
    def setUp(self):
        self.graph = {
            "ACTG": ["CTGG"],
            "CTGG": ["TGGA"],
            "TGGA": ["GGAT"],
            "TTTA": ["TTAG"],
            "TTAG": ["TAGC"],
            "TAGC": ["AGCT"],
            "AGCT": ["GCTA"],
            "GCTA": ["CTAG"],
        }
        self.edgesCount = {
            "ACTG": [0, 1],
            "CTGG": [1, 1],
            "TGGA": [1, 1],
            "GGAT": [1, 0],
            "TTTA": [0, 1],
            "TTAG": [1, 1],
            "TAGC": [1, 1],
            "AGCT": [1, 1],
            "GCTA": [1, 1],
            "CTAG": [1, 0],
        }
        self.startNodes = ["ACTG", "TTTA"]
        self.createContigs = CreateContigs(self.graph)

    def test_findStartNodes(self):
        edgesCount, startNodes = self.createContigs.findStartNodes()
        self.assertEqual(edgesCount, self.edgesCount)
        self.assertEqual(startNodes, self.startNodes)

    def test_createContigs(self):
        # os management for contig file creation and testing
        original_cwd = os.getcwd()
        os.chdir(os.path.join(original_cwd, "../src"))

        contigs = self.createContigs.createContigs()
        self.assertEqual(contigs, ["ACTGGAT", "TTTAGCTAG"])
        os.chdir(original_cwd)


if __name__ == "__main__":
    unittest.main()
