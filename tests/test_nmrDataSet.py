#!/usr/bin/env python


"""

test_nmrDataSet

author(s): C. Ludwig
origin: 26-07-2019


"""


import unittest
import metabolabpy.nmr.nmrDataSet as nmrDataSet
import os


class nmrDataSetTestCase(unittest.TestCase):
    
    def test_readSpc(self):
        pName = os.path.join(os.path.dirname(__file__),"data","nmrData") # directory of test data set
        eName = "1"                                                      # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()                                 # create nmrDataSet object
        nd.readSpc(pName,eName)                                      # check if Bruker data can be read
        self.assertEqual(len(nd.nmrdat[0][0].fid[0]), 32768)         # check number of data points in fid
        self.assertEqual(len(nd.nmrdat[0][0].spc[0]), 65536)         # check number of data points in 1r
    
if __name__ == "__main__":
    unittest.main()