#!/usr/bin/env python


"""

test_nmrDataSet

author(s): C. Ludwig
origin: 26-07-2019


"""


import unittest
import metabolabpy.nmr.acqPars as acqPars
import os


class procParsTestCase(unittest.TestCase):

    def test_readAcqPars(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData", "1")  # directory of spectrum
        ap = acqPars.AcqPars()  # create nmrDataSet object
        ap.read(pName)  # check if Bruker acqus data can be read


if __name__ == "__main__":
    unittest.main()