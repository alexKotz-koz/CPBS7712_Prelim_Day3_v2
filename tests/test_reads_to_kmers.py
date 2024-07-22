import unittest
import pandas as pd
import sys

# os management for importing readsToKmers
sys.path.insert(0, "../src/components")
from readsToKmers import ReadsToKmers


class TestReadsToKmers(unittest.TestCase):
    def setUp(self):
        readsData = []

        with open("test_data/testReadsToKmers.fastq", "r") as inputReadsFile:
            reads = inputReadsFile.readlines()
            for index, line in enumerate(reads):
                if "@" == line[0]:
                    readString = reads[index + 1]
                    plus = reads[index + 2]
                    quality = reads[index + 3]
                    readsData.append(
                        {
                            "id": line.lstrip(">").rstrip("\n"),
                            "sequence": readString.rstrip("\n"),
                        }
                    )
        # convert readsData to dataframe and add a length column that is the length of each read string
        self.dfReadsData = pd.DataFrame(readsData)
        self.dfReadsData["length"] = self.dfReadsData["sequence"].str.len()

    def testExtractKmers(self):
        rtk = ReadsToKmers(self.dfReadsData, 5)
        kmerPool = rtk.extractKmers()
        self.checkKmerPool(kmerPool=kmerPool)

    def checkKmerPool(self, kmerPool):
        self.assertDictEqual(
            kmerPool,
            {
                "ACTGG": {"@Read1": [{0: 5}]},
                "CTGGA": {"@Read1": [{1: 6}]},
                "TGGAT": {"@Read1": [{2: 7}]},
                "GGATC": {"@Read1": [{3: 8}]},
                "GATCT": {"@Read1": [{4: 9}]},
                "ATCTT": {"@Read1": [{5: 10}]},
                "TCTTC": {"@Read1": [{6: 11}]},
                "CTTCA": {"@Read1": [{7: 12}]},
                "TTCAG": {"@Read1": [{8: 13}]},
                "CTAGC": {"@Read2": [{0: 5}]},
                "TAGCC": {"@Read2": [{1: 6}]},
                "AGCCT": {"@Read2": [{2: 7}], "@Read3": [{0: 5}]},
                "GCCTT": {"@Read2": [{3: 8}], "@Read3": [{1: 6}]},
                "CCTTA": {"@Read2": [{4: 9}]},
                "CTTAT": {"@Read2": [{5: 10}]},
                "TTATC": {"@Read2": [{6: 11}]},
                "CCTTC": {"@Read3": [{2: 7}]},
                "CTTCG": {"@Read3": [{3: 8}]},
                "TTTAG": {"@Read4": [{0: 5}]},
                "TTAGC": {"@Read4": [{1: 6}]},
                "TAGCT": {"@Read4": [{2: 7}]},
                "AGCTA": {"@Read4": [{3: 8}]},
                "GCTAG": {"@Read4": [{4: 9}]},
            },
        )


if __name__ == "__main__":
    unittest.main()
