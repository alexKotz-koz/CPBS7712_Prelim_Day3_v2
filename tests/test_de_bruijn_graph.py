import unittest
import sys

sys.path.insert(0, "../src/components")
from deBruijnGraph import DeBruijnGraph


class TestDeBruijnGraph(unittest.TestCase):
    def setUp(self):
        self.kmerPool = {
            "ACTGG": {"Read1": [{0: 5}]},
            "CTGGA": {"Read1": [{1: 6}]},
            "TGGAT": {"Read1": [{2: 7}]},
            "GGATC": {"Read1": [{3: 8}]},
            "GATCT": {"Read1": [{4: 9}]},
            "ATCTT": {"Read1": [{5: 10}]},
            "TCTTC": {"Read1": [{6: 11}]},
            "CTTCA": {"Read1": [{7: 12}]},
            "TTCAG": {"Read1": [{8: 13}]},
            "CTAGC": {"Read2": [{0: 5}]},
            "TAGCC": {"Read2": [{1: 6}]},
            "AGCCT": {"Read2": [{2: 7}], "Read3": [{0: 5}]},
            "GCCTT": {"Read2": [{3: 8}], "Read3": [{1: 6}]},
            "CCTTA": {"Read2": [{4: 9}]},
            "CTTAT": {"Read2": [{5: 10}]},
            "TTATC": {"Read2": [{6: 11}]},
            "CCTTC": {"Read3": [{2: 7}]},
            "CTTCG": {"Read3": [{3: 8}]},
            "TTTAG": {"Read4": [{0: 5}]},
            "TTAGC": {"Read4": [{1: 6}]},
            "TAGCT": {"Read4": [{2: 7}]},
            "AGCTA": {"Read4": [{3: 8}]},
            "GCTAG": {"Read4": [{4: 9}]},
        }
        self.k = 5
        self.dbg = DeBruijnGraph(self.kmerPool, self.k)

    def testGetPrefixSuffix(self):
        kmer = "ACTGG"
        prefix, suffix = self.dbg.getPrefixSuffix(kmer)
        self.assertEqual(prefix, "ACTG")
        self.assertEqual(suffix, "CTGG")
        self.assertEqual(self.k - 1, len(prefix))
        self.assertEqual(self.k - 1, len(suffix))

    def testPresenceInGraph(self):
        nodes, edges = self.dbg.constructGraph()
        for kmer in self.kmerPool:
            prefix, suffix = self.dbg.getPrefixSuffix(kmer)
            assert (
                prefix in edges and suffix in edges[prefix]
            ), "The key-value pair is not in the dictionary."


if __name__ == "__main__":
    unittest.main()
