#!/usr/bin/env python


"""

test_procPars

author(s): C. Ludwig
origin: 26-07-2019


"""


import unittest
import metabolabpy.nmr.procPars as procPars
import os


class procParsTestCase(unittest.TestCase):
    
    def test_readProcPars(self):
        pName = os.path.join(os.path.dirname(__file__),"data","nmrData","1")  # directory of spectrum
        pp = procPars.ProcPars()                                              # create nmrDataSet object
        pp.read(pName)                                                        # check if Bruker procs data can be read


    def test_readVProcPars(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData", "varianData.fid")  # directory of spectrum
        pp = procPars.ProcPars()  # create nmrDataSet object
        pp.read(pName)  # check if Varian procpar data can be read


if __name__ == "__main__":
    unittest.main()