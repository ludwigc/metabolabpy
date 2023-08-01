import unittest # pragma: no cover

def suite(): # pragma: no cover
    test_loader = unittest.TestLoader() # pragma: no cover
    test_suite  = test_loader.discover('.', pattern='test_nmrDataSet.py') # pragma: no cover
    return test_suite # pragma: no cover
