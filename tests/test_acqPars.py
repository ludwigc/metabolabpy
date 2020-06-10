#!/usr/bin/env python


"""

test_acqPars

author(s): C. Ludwig
origin: 05-05-2020


"""


import unittest
import metabolabpy.nmr.acqPars as acqPars
import os


class acqParsTestCase(unittest.TestCase):

    def test_readAcqPars(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData", "1")  # directory of spectrum
        ap = acqPars.AcqPars()  # create nmrDataSet object
        ap.read(pName)  # check if Bruker acqus data can be read


if __name__ == "__main__":
    unittest.main()