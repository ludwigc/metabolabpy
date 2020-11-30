#!/usr/bin/env python


"""

test_nmrData

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
        fpMultiplier = [1.0, 0.9999405017700899, 0.9716680181001582, 0.9999931873299882, 0.9999863747063888, 0.9954327753052931]
        fidNo = int(nd.acq.groupDelay+1)
        for k in range(6):
            nd.proc.windowType[0] = k
            f = np.copy(nd.fid[0])
            v1 = f[fidNo].real
            f = nd.apodise(f, 0, 0.5, 1.0, 90.0, nd.acq.groupDelay, nd.acq.sw_h[0])
            v2 = f[fidNo].real
            self.assertAlmostEqual(v2/v1, fpMultiplier[k], 8)

    def test_autobaseline1d(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.autobaseline1d()

    def test_apcbcGetHist(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.autobaseline1d()
        h1 = nd.apc.getHist(nd.spc[0], np.linspace(-nd.apc.nMax, nd.apc.nMax, nd.apc.npts))
        self.assertEqual(len(h1), 5001)

    def test_apcbcSetVars(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.apc.setVars(nd.spc[0])
        self.assertEqual(nd.apc.npts, 65536)

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
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "2"  # 2D Jres NMR data in exp 2
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.procSpc2D()
        nd.autoRef(True)
        self.assertEqual(nd.refShift[0], 0.0)
        self.assertEqual(nd.refPoint[0], 817)
        self.assertAlmostEqual(nd.refShift[1], 4.7, 5)
        self.assertEqual(nd.refPoint[1], 64)

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
        self.assertAlmostEqual(nd.ppm1[len(nd.ppm1) - 1], 0.0)
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "2"  # 2D Jres NMR data in exp 2
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.procSpc2D()
        nd.calcPPM()
        self.assertAlmostEqual(nd.ppm1[0], 12.01657, 4)
        self.assertAlmostEqual(nd.ppm1[-1], 0.0, )
        self.assertAlmostEqual(nd.ppm2[0], 0.08333, 4)
        self.assertAlmostEqual(nd.ppm2[-1], 0.0, 4)

    def test_conv(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        fpMultiplier = [0.8537643573214553, 0.881140405836974]
        for k in range(2):
            nd.proc.convWindowType[0] = k
            f = nd.conv(nd.fid[0])
            fidNo = int(nd.acq.groupDelay + 1)
            v1 = nd.fid[0][fidNo].real
            v2 = f[fidNo].real
            self.assertAlmostEqual(v2/v1, fpMultiplier[k], 8)

    def test_fidOffsetCorrection(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.fidOffsetCorr = 128
        f = nd.fidOffsetCorrection(nd.fid[0])
        self.assertEqual(np.mean(f[:nd.fidOffsetCorr].real), 0.0)

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
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "2"  # 2D Jres NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.procSpc2D()
        m = nd.hilbert(nd.spc, 1)
        self.assertTrue(np.iscomplex(m[0][0]))

    def test_hilbert1(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        m = nd.hilbert1(nd.spc[0].real, 1)
        self.assertTrue(np.iscomplex(m[0]))
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "2"  # 2D Jres NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.procSpc2D()
        m = nd.hilbert1(nd.spc, 1)
        self.assertTrue(np.iscomplex(m[0][0]))

    def test_multiply(self):
        pName = os.path.join(os.path.dirname(__file__),"data","nmrData") # directory of test data set
        eName = "1"                                                      # 1D NMR data in exp 1
        nd               = nmrData.NmrData()
        nd.dataSetName   = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        int1 = nd.spc[0][nd.refPoint[0]].real
        nd.multiply(10.0)
        int2 = nd.spc[0][nd.refPoint[0]].real
        self.assertEqual(int2/int1, 10.0)

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
        self.assertAlmostEqual(np.sum(s1 - s2), 0.0, 4)
        nd.spc[0] = nd.spc[0].real
        s3 = np.copy(nd.spc[0].real)
        nd.phase2a(360.0, 0.0, 0)
        s4 = np.copy(nd.spc[0].real)
        self.assertAlmostEqual(np.sum(s3 - s4), 0.0, 4)
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "3"  # 2D HSQC NMR data in exp 3
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.procSpc2D()
        spc1 = np.copy(nd.spc)
        nd.phase2a(360.0, 0.0, 1)
        spc2 = np.copy(nd.spc)
        self.assertLess(np.sum(np.sum(spc1.real - spc2.real)), 1e-5)

    def test_phase2d(self):
        pName = os.path.join(os.path.dirname(__file__), "data",
                             "nmrData")  # directory of test data set
        eName = "3"  # 2D HSQC NMR data in exp 3
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.procSpc2D()
        ppm = [[]]
        ppm[0].append(1.64164)
        ppm[0].append(20.7576)
        pts = nd.ppm2points2d(ppm)
        int1 = nd.spc[int(pts[0][1])][int(pts[0][0])].real
        nd.phase2d(180.0, 0.0, 1)
        int2 = nd.spc[int(pts[0][1])][int(pts[0][0])].real
        self.assertAlmostEqual(int1, -int2)

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
        self.assertEqual(pts2 - pts1, 4682)
        nd = [[]]
        pName = os.path.join(os.path.dirname(__file__), "data",
                             "nmrData")  # directory of test data set
        eName = "3"  # 2D HSQC NMR data in exp 3
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.procSpc2D()
        pts1 = nd.ppm2points(0, 0)
        pts2 = nd.ppm2points(1, 0)
        pts3 = nd.ppm2points(0, 1)
        pts4 = nd.ppm2points(10, 1)
        self.assertEqual(pts2 - pts1, 79)
        self.assertEqual(pts4 - pts3, 223)

    def test_ppm2points2d(self):
        pName = os.path.join(os.path.dirname(__file__), "data",
                             "nmrData")  # directory of test data set
        eName = "3"  # 2D HSQC NMR data in exp 3
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.procSpc2D()
        ppm = [[]]
        ppm[0].append(1.64164)
        ppm[0].append(20.7576)
        pts = nd.ppm2points2d(ppm)
        self.assertEqual(int(pts[0][0]), 129)
        self.assertEqual(int(pts[0][1]), 462)

    def test_procSpc(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.procSpc()
        self.assertEqual(len(nd.spc[0]), 65536)         # check number of data points in Fourier transformed spectrum
        nd = [[]]
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "2"  # 2D NMR data in exp 2
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.procSpc()
        self.assertEqual(len(nd.spc), 128)  # check number of data points in Fourier transformed spectrum  (indirect dimension)
        self.assertEqual(len(nd.spc[0]), 8192)  # check number of data points in Fourier transformed spectrum (direct dimension)

    def test_procSpc1D(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.procSpc1D()
        self.assertEqual(len(nd.spc[0]), 65536)         # check number of data points in Fourier transformed spectrum

    def test_procSpc2D(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "2"  # 2D Jres NMR data in exp 2
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.procSpc2D()
        self.assertEqual(len(nd.spc), 128)  # check number of data points in Fourier transformed spectrum  (indirect dimension)
        self.assertEqual(len(nd.spc[0]), 8192)  # check number of data points in Fourier transformed spectrum (direct dimension)
        nd = [[]]
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "3"  # 2D HSQC NMR data in exp 3
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.proc.nPoints[0] = 1024
        nd.proc.nPoints[1] = 4096
        nd.procSpc2D()
        self.assertEqual(len(nd.spc), 4096)  # check number of data points in Fourier transformed spectrum  (indirect dimension)
        self.assertEqual(len(nd.spc[0]), 1024)  # check number of data points in Fourier transformed spectrum (direct dimension)
        nd = [[]]
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "4"  # 2D HSQC NMR data in exp 3
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.proc.nPoints[0] = 1024
        nd.proc.nPoints[1] = 1024
        nd.procSpc2D()
        self.assertEqual(len(nd.spc), 1024)  # check number of data points in Fourier transformed spectrum  (indirect dimension)
        self.assertEqual(len(nd.spc[0]), 1024)  # check number of data points in Fourier transformed spectrum (direct dimension)

    def test_quad2D(self):
        # test QF data (FnMODE = 1)
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "2"  # 2D Jres NMR data in exp 2
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.procSpc2D(True)
        f = np.copy(nd.quad2D(nd.fid))
        self.assertEqual(f[0][100].real, -309905)
        f = [[]]
        nd = [[]]
        # test E/A data (FnMODE = 6)
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "3"  # 2D Jres NMR data in exp 2
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.procSpc2D(True)
        f = np.copy(nd.quad2D(nd.fid))
        self.assertEqual(f[0][100].real, 1206.0)
        f = [[]]
        nd = [[]]
        # test States & States-TPPI data data (FnMODE = 5)
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "4"  # 2D Jres NMR data in exp 2
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.procSpc2D(True)
        f = np.copy(nd.quad2D(nd.fid))
        self.assertEqual(f[0][100].real, -9689.0)

    def test_readPipe2D(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "5"  # 2D HSQC NMR data rocessed with NMRPipe in exp 15
        dName = "test.dat"
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.readPipe2D(pName + os.sep + eName + '.proc', dName)
        self.assertEqual(len(nd.spc), 16384)
        self.assertEqual(len(nd.spc[0]), 922)

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

    def test_symjres(self):
        pName = os.path.join(os.path.dirname(__file__), "data",
                             "nmrData")  # directory of test data set
        eName = "2"  # 2D HSQC NMR data in exp 3
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.proc.tilt = True
        nd.proc.symj = True
        nd.procSpc2D()
        pt = np.where(np.transpose(nd.spc)[1406].real == np.amax(np.transpose(nd.spc)[1406].real))[0][0]
        self.assertEqual(pt, 23)

    def test_tiltJRes(self):
        pName = os.path.join(os.path.dirname(__file__), "data",
                             "nmrData")  # directory of test data set
        eName = "2"  # 2D HSQC NMR data in exp 3
        nd = nmrData.NmrData()
        nd.dataSetName = pName
        nd.dataSetNumber = eName
        nd.readSpc()
        nd.proc.tilt = True
        nd.proc.symj = False
        nd.procSpc2D()
        pt = np.where(np.transpose(nd.spc)[1406].real == np.amax(np.transpose(nd.spc)[1406].real))[0][0]
        self.assertEqual(pt, 78)

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
