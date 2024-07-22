import os
import sys

sys.path.insert(0, "../src/components")
import unittest
from importVirus import ImportVirus


class TestImportVirus(unittest.TestCase):
    def setUp(self):
        self.importVirus = ImportVirus()
        self.virusFile = ["test_data/testviruses.fasta"]

    def test_importVirusData(self):
        result = self.importVirus.importVirusData(self.virusFile)
        self.assertIsInstance(result, dict)
        if result:
            firstKey = list(result.keys())[0]
            self.assertIn("name", result[firstKey])
            self.assertIn("sequence", result[firstKey])


if __name__ == "__main__":
    unittest.main()
