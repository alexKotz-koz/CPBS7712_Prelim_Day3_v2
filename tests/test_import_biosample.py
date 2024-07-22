import os
import sys

sys.path.insert(0, "../src/components")
import unittest
import pandas as pd
from importBioSample import ImportBioSample


class TestImportBioSample(unittest.TestCase):
    def setUp(self):
        self.file = "testbiosample.fastq"
        self.import_biosample = ImportBioSample(self.file)

    def test_importBioSample(self):
        resultDict = self.import_biosample.importBioSample()
        fileExtension = os.path.splitext(self.file)[-1]
        self.assertIsInstance(resultDict, dict)
        self.assertTrue(fileExtension == ".fastq")
        if resultDict:
            firstKey = list(resultDict.keys())[0]
            self.assertIn("sequence", resultDict[firstKey])
            self.assertIn("quality", resultDict[firstKey])


if __name__ == "__main__":
    unittest.main()
