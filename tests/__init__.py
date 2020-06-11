import unittest

def suite():
    test_loader = unittest.TestLoader()
    test_suite  = test_loader.discover('.', pattern='test_nmrDataSet.py')
    return test_suite
