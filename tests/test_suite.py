import unittest
from test_reads_to_kmers import TestReadsToKmers
from test_de_bruijn_graph import TestDeBruijnGraph
from test_create_contigs import TestCreateContigs
from test_search_string import TestSearchString
from test_qc import TestQualityControl
from test_import_biosample import TestImportBioSample
from test_import_virus import TestImportVirus


# runs all tests
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestImportBioSample))
    suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestImportVirus))
    suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestQualityControl))
    suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestReadsToKmers))
    suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestDeBruijnGraph))
    suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestCreateContigs))
    suite.addTest(unittest.makeSuite(TestSearchString))

    return suite


if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    runner.run(suite())
