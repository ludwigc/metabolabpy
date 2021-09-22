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

    def test_read_acq_pars(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData", "1")  # directory of spectrum
        ap = acqPars.AcqPars()  # create nmrDataSet object
        ap.read(pName)  # check if Bruker acqus data can be read

    def test_read_v_acq_pars(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData", "varianData.fid")  # directory of spectrum
        ap = acqPars.AcqPars()  # create nmrDataSet object
        ap.read(pName)  # check if Varian procpar data can be read


if __name__ == "__main__":
    unittest.main()