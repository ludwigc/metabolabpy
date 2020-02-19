#!/usr/bin/env python


"""

test_nmrDataSet

author(s): C. Ludwig
origin: 26-07-2019


"""


import unittest
import metabolabpy.nmr.nmrDataSet as nmrDataSet
import os
import numpy as np


class nmrDataSetTestCase(unittest.TestCase):
    
    def test_readSpc(self):
        pName = os.path.join(os.path.dirname(__file__),"data","nmrData") # directory of test data set
        eName = "1"                                                      # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()                                 # create nmrDataSet object
        nd.readSpc(pName,eName)                                      # check if Bruker data can be read
        self.assertEqual(len(nd.nmrdat[0][0].fid[0]), 32768)         # check number of data points in fid
        self.assertEqual(len(nd.nmrdat[0][0].spc[0]), 65536)         # check number of data points in 1r

    def test_ft(self):
        pName = os.path.join(os.path.dirname(__file__),"data","nmrData") # directory of test data set
        eName = "1"                                                      # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()                                 # create nmrDataSet object
        nd.readSpc(pName,eName)                                      # check if Bruker data can be read
        nd.ft()
        self.assertEqual(len(nd.nmrdat[0][0].spc[0]), 65536)         # check number of data points in Fourier transformed spectrum

    def test_setZeroFill(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.setZeroFill([131072])
        nd.ft()
        self.assertEqual(len(nd.nmrdat[0][0].spc[0]),131072)  # check number of data points in Fourier transformed spectrum

    def test_autobaseline1d(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.autobaseline1d()

    def test_autophase1d(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.autophase1d()

    def test_autoref(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.autoref()

    def test_clear(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.clear()
        self.assertEqual(len(nd.nmrdat),1)
        self.assertEqual(len(nd.nmrdat[0]),0)

    def test_excludeRegion(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.pp.excludeStart = [8.9]
        nd.pp.excludeEnd   = [11.0]
        nd.pp.flagExcludeRegion = True
        nd.dataPreProcessing()
        pts = len(nd.nmrdat[0][0].spc[0]) - nd.nmrdat[0][0].ppm2points([nd.pp.excludeStart, nd.pp.excludeEnd],0)
        pts2 = range(int(min(pts)), int(max(pts)))
        ssum = np.sum(nd.nmrdat[0][0].spc[0][pts2].real)
        self.assertEqual(ssum,0.0)

    def test_noiseFiltering(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.pp.noiseStart     = 10.0
        nd.pp.noiseEnd       = 10.5
        nd.pp.noiseThreshold = 8.0
        nd.pp.thLineWidth    = 2.0
        nd.pp.flagNoiseFiltering = True
        nd.dataPreProcessing()
        pts = len(nd.nmrdat[0][0].spc[0]) - nd.nmrdat[0][0].ppm2points([9.5, 11.5],0)
        pts2 = range(int(min(pts)), int(max(pts)))
        ssum = np.sum(nd.nmrdat[0][0].spc[0][pts2].real)
        self.assertEqual(ssum,0.0)

    def test_bucketSpectra(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.pp.bucketPoints = 64
        nd.pp.flagBucketSpectra = True
        nd.dataPreProcessing()
        pts = len(nd.nmrdat[0][0].spc[0])
        self.assertEqual(pts,1024)

    def test_setGb(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.setGb([15.0])
        self.assertEqual(nd.nmrdat[0][0].proc.gb[0],15.0)

    def test_setLb(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.setLb([15.0])
        self.assertEqual(nd.nmrdat[0][0].proc.lb[0],15.0)

    def test_setPhFromExp(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.setPhFromExp()
        self.assertEqual(nd.nmrdat[0][0].proc.ph0[0],22.0)
        self.assertEqual(nd.nmrdat[0][0].proc.ph1[0],22.0)

    def test_ph0(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.setPh0([15.0])
        self.assertEqual(nd.nmrdat[0][0].proc.ph0[0],15.0)

    def test_ph1(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.setPh1([15.0])
        self.assertEqual(nd.nmrdat[0][0].proc.ph1[0],15.0)

    def test_setSsb(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.setSsb([1.0])
        self.assertEqual(nd.nmrdat[0][0].proc.ssb[0],1.0)

    def test_setWindowType(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.setWindowType([2])
        self.assertEqual(nd.nmrdat[0][0].proc.windowType[0],2)

    def test_shiftRef(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.shiftRef()
        self.assertEqual(nd.nmrdat[0][0].refShift[0],0.0)




if __name__ == "__main__":
    unittest.main()