#!/usr/bin/env python


"""

test_nmrDataSet

author(s): C. Ludwig
origin: 26-07-2019


"""


import unittest
import metabolabpy.nmr.nmrData as nmrData
import os
import numpy as np


class nmrDataTestCase(unittest.TestCase):
    
    def test_apodise(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        for k in range(6):
            f = nd.apodise(nd.fid, 0, 0.5, 0.0, 0.0, nd.acq.groupDelay, nd.acq.sw_h[0])


    def test_autobaseline1d(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.autobaseline1d()

    def test_autophase1d(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.autophase1d()

    def test_autoRef(self, tmsp=True):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.autoRef(True)
        self.assertEqual(nd.refShift[0], 0.0)
        self.assertEqual(nd.refPoint[0], 10342)

    #def test_baseline1d(self):

    def test_calcPPM(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.calcPPM()
        self.assertAlmostEqual(nd.ppm1[0], 13.996205596)
        self.assertAlmostEqual(nd.ppm1[len(nd.ppm1)-1], 0.0)

    def test_conv(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        f = nd.conv(nd.fid[0])

    def test_fidOffsetCorrection(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        f = nd.fidOffsetCorrection(nd.fid[0])


    def test_gibbs(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        f = nd.gibbs(nd.fid[0])


    def test_hilbert(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        m = nd.hilbert(nd.spc.real, 0)
        self.assertTrue(np.iscomplex(m[0][0]))

    def test_hilbert1(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        m = nd.hilbert1(nd.spc[0].real, 0)
        self.assertTrue(np.iscomplex(m[0]))

    def test_phase(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        s1 = np.copy(nd.spc[0].real)
        nd.phase(360.0, 0.0, len(nd.spc[0]))
        s2 = np.copy(nd.spc[0].real)
        self.assertAlmostEqual(np.sum(s1-s2), 0.0, 4)

    def test_phase2(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        s1 = np.copy(nd.spc[0].real)
        s2 = nd.phase2(nd.spc[0], 360.0, 0.0)
        s2 = s2.real
        self.assertAlmostEqual(np.sum(s1-s2), 0.0, 4)

    def test_phase2a(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        s1 = np.copy(nd.spc[0].real)
        nd.phase2a(360.0, 0.0, 0)
        s2 = np.copy(nd.spc[0].real)
        self.assertAlmostEqual(np.sum(s1-s2), 0.0, 4)

    #def test_phase2d(self, ph0, ph1, dim):

    def test_phase3(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        s1 = np.copy(nd.spc[0].real)
        s2 = nd.phase3(nd.spc[0], 360.0, 0.0)
        s2 = s2.real
        self.assertAlmostEqual(np.sum(s1-s2), 0.0, 4)

    def test_points2Hz(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        hz1 = nd.points2Hz(0, 0)
        hz2 = nd.points2Hz(1, 0)
        self.assertAlmostEqual(hz2 - hz1, 0.128, 3)

    def test_points2ppm(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        ppm1 = nd.points2ppm(0, 0)
        ppm2 = nd.points2ppm(1, 0)
        self.assertAlmostEqual(ppm2 - ppm1, 0.0002135684)

    def test_ppm2points(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        pts1 = nd.ppm2points(0, 0)
        pts2 = nd.ppm2points(1, 0)
        self.assertAlmostEqual(pts2 - pts1, 4682)

    #def test_ppm2points2d(self, ppm):

    def test_procSpc(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.procSpc()
        self.assertEqual(len(nd.spc[0]), 65536)         # check number of data points in Fourier transformed spectrum

    def test_procSpc1D(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.procSpc1D()
        self.assertEqual(len(nd.spc[0]), 65536)         # check number of data points in Fourier transformed spectrum


    #def test_procSpc2D(self):
    #def test_quad2D(self, fid):
    #def test_readPipe2D(self, pName, fName):

    def test_readSpc(self):
        pName = os.path.join(os.path.dirname(__file__),"data","nmrData") # directory of test data set
        eName = "1"                                                      # 1D NMR data in exp 1
        nd               = nmrData.NmrData()
        nd.dataSetName   = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        self.assertEqual(len(nd.fid[0]), 32768)         # check number of data points in fid
        self.assertEqual(len(nd.spc[0]), 65536)         # check number of data points in 1r

    def test_setRef(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.setRef([-10.0], [2048])
        self.assertEqual(nd.refShift[0], -10.0)
        self.assertEqual(nd.refPoint[0], 2048)

    def test_setWindowFunction(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.setWindowFunction(0, 'none')
        self.assertEqual(nd.proc.windowType[0],0)
        nd.setWindowFunction(0, 'exponential')
        self.assertEqual(nd.proc.windowType[0],1)
        nd.setWindowFunction(0, 'gaussian')
        self.assertEqual(nd.proc.windowType[0],2)
        nd.setWindowFunction(0, 'sine')
        self.assertEqual(nd.proc.windowType[0],3)
        nd.setWindowFunction(0, 'qsine')
        self.assertEqual(nd.proc.windowType[0],4)
        nd.setWindowFunction(0, 'sem')
        self.assertEqual(nd.proc.windowType[0],5)

    def test_smo(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        f = nd.smo(nd.fid[0])

    #def test_symjres(self):
    #def test_tiltJRes(self):

    def test_waterSupp(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.proc.waterSuppression = 1
        f = nd.waterSupp(nd.spc[0])
        nd.proc.waterSuppression = 2
        f = nd.waterSupp(nd.spc[0])

    def test_zeroFill(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.proc.nPoints[0] = 131072
        f = nd.zeroFill(nd.fid[0], 0)
        npts2 = len(f)
        self.assertEqual(npts2, 131072)


if __name__ == "__main__":
    unittest.main()