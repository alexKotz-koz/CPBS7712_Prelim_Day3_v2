import os
import sys


sys.path.insert(0, "../src/components")
import json
import unittest
import pandas as pd

from qc import QualityControl


class TestQualityControl(unittest.TestCase):
    def setUp(self):

        with open("test_data/biosample_ex.json") as file:
            self.test_data = json.load(file)

    def test_qualityControl(self):
        qc = QualityControl(self.test_data)

        returnedChar = qc.asciiToPhred("A")
        self.assertEqual(returnedChar, 32)

        biosampleDf, minSeqLength, qcReport, qcMetaData = qc.qualityControl()

        # Assert the output data
        self.assertIsInstance(biosampleDf, pd.DataFrame)
        self.assertIsInstance(minSeqLength, int)
        self.assertIsInstance(qcReport, dict)
        # Add more assertions as needed
        if qcReport:
            firstKey = list(qcReport.keys())[0]
            self.assertIn("meanQ", qcReport[firstKey])
            self.assertIn("medianQ", qcReport[firstKey])
            self.assertIn("length", qcReport[firstKey])


if __name__ == "__main__":
    unittest.main()
