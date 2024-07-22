import json
import unittest
import sys
from unittest.mock import patch, mock_open

sys.path.insert(0, "../src/components")
from searchForViruses import SearchString


class TestSearchString(unittest.TestCase):
    def setUp(self):
        self.viruses = {"virus1": {"sequence": "ATCG"}}
        self.readsKmerPoolFile = "test_data/testkmerpool.json"
        self.contigs = ["ATCG", "CGTA"]
        self.k = 2
        self.searchString = SearchString(
            self.viruses, self.readsKmerPoolFile, self.contigs, self.k
        )
        with patch(
            "builtins.open",
            mock_open(
                read_data=json.dumps(
                    {
                        "GGCCAAGCTCGATGCCGAAATCAAGGCCAGGGCCGTAGACATCAA": {
                            "@SRR24581287.201": [{"0": 45}]
                        },
                        "GCCAAGCTCGATGCCGAAATCAAGGCCAGGGCCGTAGACATCAAC": {
                            "@SRR24581287.201": [{"1": 46}]
                        },
                    }
                )
            ),
        ):
            self.searchString = SearchString(
                self.viruses, self.readsKmerPoolFile, self.contigs, self.k
            )

    def test_virusToKmers(self):
        sequence = "ATCG"
        expected = {"AT": [0], "TC": [1], "CG": [2]}
        result = self.searchString.virusToKmers(sequence)
        self.assertEqual(result, expected)

    def test_hammingDistance(self):
        virus = "ATCG"
        contig = "ATGC"
        expected = 2
        result = self.searchString.hammingDistance(virus, contig)
        self.assertEqual(result, expected)


if __name__ == "__main__":
    unittest.main()
