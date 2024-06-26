#!/usr/bin/env python


"""

test_nmrData

author(s): C. Ludwig
origin: 26-07-2019


"""

import unittest
import os
import numpy as np
from metabolabpy.nmr import nmrData
from metabolabpy.nmr import nmrDataSet
from metabolabpy.nmr.phase3 import *
import pandas as pd
import matplotlib.pyplot as pl

class nmrDataTestCase(unittest.TestCase):

    def test_apodise(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        fp_multiplier = [1.0, 0.9999405017700899, 0.9716680181001582, 0.9999931873299882, 0.9999863747063888,
                         0.9954327753052931]
        fid_no = int(nd.acq.group_delay + 1)
        for k in range(6):
            nd.proc.window_type[0] = k
            f = np.copy(nd.fid[0])
            v1 = f[fid_no].real
            f = nd.apodise(f, 0, 0.5, 1.0, 90.0, nd.acq.group_delay, nd.acq.sw_h[0])
            v2 = f[fid_no].real
            self.assertAlmostEqual(v2 / v1, fp_multiplier[k], 3)
    
    def test_autobaseline1d(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.autobaseline1d()
    
    def test_apcbc_get_hist(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.autobaseline1d()
        h1 = nd.apc.get_hist(nd.spc[0], np.linspace(-nd.apc.n_max, nd.apc.n_max, nd.apc.npts))
        self.assertEqual(len(h1), 5001)
    
    def test_apcbc_set_vars(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.apc.set_vars(nd.spc[0])
        self.assertEqual(nd.apc.npts, 65536)
    
    def test_autophase1d(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.autophase1d()

    def test_auto_ref(self, tmsp=True):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.auto_ref(True)
        self.assertEqual(nd.ref_shift[0], 0.0)
        self.assertEqual(nd.ref_point[0], 10342)
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "2"  # 2D Jres NMR data in exp 2
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.proc_spc2d()
        nd.auto_ref(True)
        self.assertEqual(nd.ref_shift[0], 0.0)
        self.assertEqual(nd.ref_point[0], 817)
        self.assertAlmostEqual(nd.ref_shift[1], 4.7, 5)
        self.assertEqual(nd.ref_point[1], 64)

    # def test_baseline1d(self):

    def test_calc_ppm(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.calc_ppm()
        self.assertAlmostEqual(nd.ppm1[0], 13.996205596)
        self.assertAlmostEqual(nd.ppm1[len(nd.ppm1) - 1], 0.0)
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "2"  # 2D Jres NMR data in exp 2
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.proc_spc2d()
        nd.calc_ppm()
        self.assertAlmostEqual(nd.ppm1[0], 12.01657, 4)
        self.assertAlmostEqual(nd.ppm1[-1], 0.0, )
        self.assertAlmostEqual(nd.ppm2[0], 0.08333, 4)
        self.assertAlmostEqual(nd.ppm2[-1], 0.0, 4)
    
    def test_conv(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        fp_multiplier = [0.8537643573214553, 0.881140405836974]
        for k in range(2):
            nd.proc.conv_window_type[0] = k
            f = nd.conv(nd.fid[0])
            fid_no = int(nd.acq.group_delay + 1)
            v1 = nd.fid[0][fid_no].real
            v2 = f[fid_no].real
            self.assertAlmostEqual(v2 / v1, fp_multiplier[k], 8)
    
    def test_fid_offset_correction(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.fid_offset_corr = 128
        f = nd.fid_offset_correction(nd.fid[0])
        self.assertEqual(np.mean(f[:nd.fid_offset_corr].real), 0.0)
    
    def test_gibbs(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        f = nd.gibbs(nd.fid[0])
    
    def test_hilbert(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        m = nd.hilbert(nd.spc.real, 0)
        self.assertTrue(np.iscomplex(m[0][0]))
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "2"  # 2D Jres NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.proc_spc2d()
        m = nd.hilbert(nd.spc, 1)
        self.assertTrue(np.iscomplex(m[0][0]))
    
    def test_hilbert1(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        m = nd.hilbert1(nd.spc[0].real, 1)
        self.assertTrue(np.iscomplex(m[0]))
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "2"  # 2D Jres NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.proc_spc2d()
        m = nd.hilbert1(nd.spc, 1)
        self.assertTrue(np.iscomplex(m[0][0]))
    
    def test_multiply(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        int1 = nd.spc[0][nd.ref_point[0]].real
        nd.multiply(10.0)
        int2 = nd.spc[0][nd.ref_point[0]].real
        self.assertEqual(int2 / int1, 10.0)
    
    def test_phase(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        s1 = np.copy(nd.spc[0].real)
        nd.phase(360.0, 0.0, len(nd.spc[0]))
        s2 = np.copy(nd.spc[0].real)
        self.assertAlmostEqual(np.sum(s1 - s2), 0.0, 4)
    
    def test_phase2(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        s1 = np.copy(nd.spc[0].real)
        s2 = nd.phase2(nd.spc[0], 360.0, 0.0)
        s2 = s2.real
        self.assertAlmostEqual(np.sum(s1 - s2), 0.0, 4)
    
    def test_phase2a(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        s1 = np.copy(nd.spc[0].real)
        nd.phase2a(360.0, 0.0, 0)
        s2 = np.copy(nd.spc[0].real)
        self.assertAlmostEqual(np.sum(s1 - s2), 0.0, 4)
        nd.spc[0] = nd.spc[0].real
        s3 = np.copy(nd.spc[0].real)
        nd.phase2a(360.0, 0.0, 0)
        s4 = np.copy(nd.spc[0].real)
        self.assertAlmostEqual(np.sum(s3 - s4), 0.0, 4)
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "3"  # 2D HSQC NMR data in exp 3
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.proc_spc2d()
        spc1 = np.copy(nd.spc)
        nd.phase2a(360.0, 0.0, 1)
        spc2 = np.copy(nd.spc)
        self.assertLess(np.sum(np.sum(spc1.real - spc2.real)), 1e-5)
    
    def test_phase2d(self):
        p_name = os.path.join(os.path.dirname(__file__), "data",
                              "nmrData")  # directory of test data set
        e_name = "3"  # 2D HSQC NMR data in exp 3
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.proc_spc2d()
        ppm = [[]]
        ppm[0].append(1.64164)
        ppm[0].append(20.7576)
        pts = nd.ppm2points2d(ppm)
        int1 = nd.spc[int(pts[0][1])][int(pts[0][0])].real
        nd.phase2d(180.0, 0.0, 1)
        int2 = nd.spc[int(pts[0][1])][int(pts[0][0])].real
        self.assertAlmostEqual(int1, -int2)
    
    def test_phase3(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        s1 = np.copy(nd.spc[0].real)
        s2 = phase3(nd.spc[0], 360.0, 0.0, len(nd.spc[0]))
        s2 = s2.real
        s2 -= (s2[0] - s1[0])
        self.assertAlmostEqual(np.sum(s1 - s2), 0.0, 4)
    
    def test_points2hz(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        hz1 = nd.points2hz(0, 0)
        hz2 = nd.points2hz(1, 0)
        self.assertAlmostEqual(hz2 - hz1, 0.128, 3)
    
    def test_points2ppm(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        ppm1 = nd.points2ppm(0, 0)
        ppm2 = nd.points2ppm(1, 0)
        self.assertAlmostEqual(ppm2 - ppm1, 0.0002135684)
    
    def test_ppm2points(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        pts1 = nd.ppm2points(0, 0)
        pts2 = nd.ppm2points(1, 0)
        self.assertEqual(pts2 - pts1, 4682)
        nd = [[]]
        p_name = os.path.join(os.path.dirname(__file__), "data",
                              "nmrData")  # directory of test data set
        e_name = "3"  # 2D HSQC NMR data in exp 3
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.proc_spc2d()
        pts1 = nd.ppm2points(0, 0)
        pts2 = nd.ppm2points(1, 0)
        pts3 = nd.ppm2points(0, 1)
        pts4 = nd.ppm2points(10, 1)
        self.assertEqual(pts2 - pts1, 157)
        self.assertEqual(pts4 - pts3, 254)
    
    def test_ppm2points2d(self):
        p_name = os.path.join(os.path.dirname(__file__), "data",
                              "nmrData")  # directory of test data set
        e_name = "3"  # 2D HSQC NMR data in exp 3
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.proc_spc2d()
        ppm = [[]]
        ppm[0].append(1.64164)
        ppm[0].append(20.7576)
        pts = nd.ppm2points2d(ppm)
        self.assertEqual(int(pts[0][0]), 258)
        self.assertEqual(int(pts[0][1]), 528)
    
    def test_proc_spc(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.proc_spc()
        self.assertEqual(len(nd.spc[0]), 65536)  # check number of data points in Fourier transformed spectrum
        nd = [[]]
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "2"  # 2D NMR data in exp 2
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.proc_spc()
        self.assertEqual(len(nd.spc),
                         128)  # check number of data points in Fourier transformed spectrum  (indirect dimension)
        self.assertEqual(len(nd.spc[0]),
                         8192)  # check number of data points in Fourier transformed spectrum (direct dimension)
    
    def test_proc_spc1d(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.proc_spc1d()
        self.assertEqual(len(nd.spc[0]), 65536)  # check number of data points in Fourier transformed spectrum
    
    def test_proc_spc2d(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "2"  # 2D Jres NMR data in exp 2
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.proc_spc2d()
        self.assertEqual(len(nd.spc),
                         128)  # check number of data points in Fourier transformed spectrum  (indirect dimension)
        self.assertEqual(len(nd.spc[0]),
                         8192)  # check number of data points in Fourier transformed spectrum (direct dimension)
        nd = [[]]
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "3"  # 2D HSQC NMR data in exp 3
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.proc.n_points[0] = 1024
        nd.proc.n_points[1] = 4096
        nd.proc_spc2d()
        self.assertEqual(len(nd.spc),
                         4096)  # check number of data points in Fourier transformed spectrum  (indirect dimension)
        self.assertEqual(len(nd.spc[0]),
                         2048)  # check number of data points in Fourier transformed spectrum (direct dimension)
        nd = [[]]
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "4"  # 2D HSQC NMR data in exp 3
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.proc.n_points[0] = 1024
        nd.proc.n_points[1] = 1024
        nd.proc_spc2d()
        self.assertEqual(len(nd.spc),
                         1024)  # check number of data points in Fourier transformed spectrum  (indirect dimension)
        self.assertEqual(len(nd.spc[0]),
                         1024)  # check number of data points in Fourier transformed spectrum (direct dimension)
    
    def test_quad2d(self):
        # test QF data (FnMODE = 1)
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "2"  # 2D Jres NMR data in exp 2
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.proc_spc2d(True)
        f = np.copy(nd.quad_2d(nd.fid))
        self.assertEqual(f[0][100].real, -309905)
        f = [[]]
        nd = [[]]
        # test E/A data (FnMODE = 6)
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "3"  # 2D Jres NMR data in exp 2
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.proc_spc2d(True)
        f = np.copy(nd.quad_2d(nd.fid))
        self.assertEqual(f[0][100].real, 1206.0)
        f = [[]]
        nd = [[]]
        # test States & States-TPPI data data (FnMODE = 5)
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "4"  # 2D Jres NMR data in exp 2
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.proc_spc2d(True)
        f = np.copy(nd.quad_2d(nd.fid))
        self.assertEqual(f[0][100].real, -9689.0)
    
    def test_read_pipe_2d(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "5"  # 2D HSQC NMR data processed with NMRPipe in exp 15
        d_name = "test.dat"
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.read_pipe_2d(p_name + os.sep + e_name + '.proc', d_name)
        self.assertEqual(len(nd.spc), 16384)
        self.assertEqual(len(nd.spc[0]), 922)
    
    def test_read_spc(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        self.assertEqual(len(nd.fid[0]), 32768)  # check number of data points in fid
        self.assertEqual(len(nd.spc[0]), 65536)  # check number of data points in 1r
        nd = [[]]
        nd = nmrData.NmrData()
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "varianData.fid"
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.proc_spc1d()
        self.assertEqual(len(nd.fid[0]), 16384)  # check number of data points in fid
        self.assertEqual(len(nd.spc[0]), 32768)  # check number of data points in 1r
        # end test_read_spc
    
    def test_set_ref(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.set_ref([-10.0], [2048])
        self.assertEqual(nd.ref_shift[0], -10.0)
        self.assertEqual(nd.ref_point[0], 2048)
    
    def test_set_window_function(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.set_window_function(0, 'none')
        self.assertEqual(nd.proc.window_type[0], 0)
        nd.set_window_function(0, 'exponential')
        self.assertEqual(nd.proc.window_type[0], 1)
        nd.set_window_function(0, 'gaussian')
        self.assertEqual(nd.proc.window_type[0], 2)
        nd.set_window_function(0, 'sine')
        self.assertEqual(nd.proc.window_type[0], 3)
        nd.set_window_function(0, 'qsine')
        self.assertEqual(nd.proc.window_type[0], 4)
        nd.set_window_function(0, 'sem')
        self.assertEqual(nd.proc.window_type[0], 5)
    
    def test_smo(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        f = nd.smo(nd.fid[0])
    
    def test_sym_jres(self):
        p_name = os.path.join(os.path.dirname(__file__), "data",
                              "nmrData")  # directory of test data set
        e_name = "2"  # 2D HSQC NMR data in exp 3
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.proc.tilt = True
        nd.proc.symj = True
        nd.proc_spc2d()
        pt = np.where(np.transpose(nd.spc)[1406].real == np.amax(np.transpose(nd.spc)[1406].real))[0][0]
        self.assertEqual(pt, 23)
    
    def test_tilt_jres(self):
        p_name = os.path.join(os.path.dirname(__file__), "data",
                              "nmrData")  # directory of test data set
        e_name = "2"  # 2D HSQC NMR data in exp 3
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.proc.tilt = True
        nd.proc.symj = False
        nd.proc_spc2d()
        pt = np.where(np.transpose(nd.spc)[1406].real == np.amax(np.transpose(nd.spc)[1406].real))[0][0]
        self.assertEqual(pt, 78)
    
    def test_water_supp(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.proc.water_suppression = 1
        f = nd.water_supp(nd.spc[0])
        nd.proc.water_suppression = 2
        f = nd.water_supp(nd.spc[0])
    
    def test_zero_fill(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.proc.n_points[0] = 131072
        f = nd.zero_fill(nd.fid[0], 0)
        npts2 = len(f)
        self.assertEqual(npts2, 131072)
    
    def test_acq_pars(self):
        f_name = os.path.join(os.path.dirname(__file__), "data", "loadData.mlpy")  # directory of test data set
        nd = nmrDataSet.NmrDataSet()
        nd.load(f_name)
        self.assertEqual(nd.nmrdat[0][0].acq.acqus_text.index('Z814601_0087'), 5005)
        text = nd.nmrdat[0][0].acq_pars()
        self.assertEqual(text.index('Z814601_0087'), 446)

        # end

    def test_add_baseline_points(self):
        f_name = os.path.join(os.path.dirname(__file__), "data", "loadData.mlpy")  # directory of test data set
        nd = nmrDataSet.NmrDataSet()
        nd.load(f_name)
        nd.nmrdat[0][0].spline_baseline.baseline_points = [11.787481124553866, -2.2087244721275536]
        nd.nmrdat[0][0].add_baseline_points()
        self.assertEqual(nd.nmrdat[0][0].spline_baseline.baseline_points[0], 11.7768)
        self.assertEqual(nd.nmrdat[0][0].spline_baseline.baseline_points[1], -2.1978    )
    # end
    
    def test_add_peak(self):
        f_name = os.path.join(os.path.dirname(__file__), "data", "loadData.mlpy")  # directory of test data set
        nd = nmrDataSet.NmrDataSet()
        nd.load(f_name)
        nd.nmrdat[0][0].add_peak(np.array([0.01, -0.01]))
        self.assertAlmostEqual(nd.nmrdat[0][0].peak_max_ppm, 0.0)
    # end
    
    def test_set_peak(self):
        f_name = os.path.join(os.path.dirname(__file__), "data", "loadData.mlpy")  # directory of test data set
        nd = nmrDataSet.NmrDataSet()
        nd.load(f_name)
        nd.nmrdat[0][0].set_peak(np.array([0.01]), np.array([-0.01]), np.array(['TMSP']))
        self.assertAlmostEqual(nd.nmrdat[0][0].peak_max_ppm[0], 0.0)
    # end
    
    def test_clear_peak(self):
        f_name = os.path.join(os.path.dirname(__file__), "data", "loadData.mlpy")  # directory of test data set
        nd = nmrDataSet.NmrDataSet()
        nd.load(f_name)
        nd.nmrdat[0][0].set_peak(np.array([0.01]), np.array([-0.01]), np.array(['TMSP']))
        nd.nmrdat[0][0].clear_peak()
        self.assertEqual(len(nd.nmrdat[0][0].peak_max_ppm), 0)
    # end
    
    def test_add_tmsp(self):
        f_name = os.path.join(os.path.dirname(__file__), "data", "loadData.mlpy")  # directory of test data set
        nd = nmrDataSet.NmrDataSet()
        nd.load(f_name)
        ppm_values = [0, 0.01, -0.01]
        points = len(nd.nmrdat[0][0].spc[0]) - nd.nmrdat[0][0].ppm2points(ppm_values)
        nd.nmrdat[0][0].spc[0][points[1]:points[2]] = np.zeros(points[2] - points[1])
        nd.nmrdat[0][0].add_tmsp(1e5, 0.0001)
        nd.nmrdat[0][0].add_peak(np.array([0.01, -0.01]))
        self.assertAlmostEqual(nd.nmrdat[0][0].peak_max_ppm[0], 0, places=3)
    # end
    
    def test_autobaseline1d(self):
        f_name = os.path.join(os.path.dirname(__file__), "data", "loadData.mlpy")  # directory of test data set
        nd = nmrDataSet.NmrDataSet()
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='notThere')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='airpls')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='arpls')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='asls')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='aspls')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='derpsalsa')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='drpls')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='iarpls')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='iasls')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='psalsa')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='mpls')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='mor')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='imor')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='mormol')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='amormol')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='rolling_ball')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='mwmv')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='tophat')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='mpspline')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='jbcd')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='noise_median')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='snip')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='swima')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='ipsa')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='ria')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='dietrich')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='std_distribution')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='fastchrom')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='cwt_br')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='fabc')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='beads')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='poly')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='modpoly')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='imodpoly')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='penalized_poly')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='quant_reg')
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d(alg='goldindec')
    
    def test_autobaseline1d_old(self):
        f_name = os.path.join(os.path.dirname(__file__), "data", "loadData.mlpy")  # directory of test data set
        nd = nmrDataSet.NmrDataSet()
        nd.load(f_name)
        nd.nmrdat[0][0].autobaseline1d_old()
    
    def test_autophase1d1(self):
        f_name = os.path.join(os.path.dirname(__file__), "data", "loadData.mlpy")  # directory of test data set
        nd = nmrDataSet.NmrDataSet()
        nd.load(f_name)
        nd.nmrdat[0][0].exclude_water = True
        nd.nmrdat[0][0].autophase1d1()
        nd.nmrdat[0][0].exclude_water = False
        nd.nmrdat[0][0].autophase1d1()
        nd.nmrdat[0][0].set_peak(np.array([0.01]), np.array([-0.01]), np.array(['TMSP']))
        self.assertAlmostEqual(nd.nmrdat[0][0].peak_max_ppm[0], 0.0, places=3)
    
    def test_create_title(self):
        excel_name = os.path.join(os.path.dirname(__file__), "data", "sampleTitleSpreadSheet.xlsx")
        xls = pd.read_excel(excel_name).fillna('')
        f_name = os.path.join(os.path.dirname(__file__), "data", "loadData.mlpy")  # directory of test data set
        nd = nmrDataSet.NmrDataSet()
        nd.load(f_name)
        rack_label = 'Rack'
        pos_label = 'Position'
        dataset_label = 'dataset'
        replace_title = True
        c_dict = {}
        for k in range(len(xls[pos_label])):
            if str(xls[pos_label][k]) != 'nan':
                if not xls[dataset_label][k] in c_dict.keys():
                    c_dict[xls[dataset_label][k]] = {}

                c_dict[xls[dataset_label][k]][str(xls[rack_label][k]) + " " + str(xls[pos_label][k])] = k

        nd.nmrdat[0][0].create_title(xls=xls, dataset_label=dataset_label, pos_label=pos_label, rack_label=rack_label, replace_title=replace_title, c_dict=c_dict, excel_name=excel_name)
        nd.nmrdat[0][0].title.index('sample : control')
        # end test_create_title
    

    def test_calc_spline_baseline(self):
        f_name = os.path.join(os.path.dirname(__file__), "data", "loadData.mlpy")  # directory of test data set
        nd = nmrDataSet.NmrDataSet()
        nd.load(f_name)
        nd.nmrdat[0][0].spline_baseline.baseline_points = [11.787481124553866, -2.2087244721275536]
        nd.nmrdat[0][0].add_baseline_points()
        baseline = nd.nmrdat[0][0].calc_spline_baseline()
        self.assertAlmostEqual(baseline[51], 6854.771940828485, places=1)
        self.assertAlmostEqual(baseline[65485], 4821.232145140779, places=1)
        # end test_calc_spline_baseline
    
    def test_corr_spline_baseline(self):
        f_name = os.path.join(os.path.dirname(__file__), "data", "loadData.mlpy")  # directory of test data set
        nd = nmrDataSet.NmrDataSet()
        nd.load(f_name)
        nd.nmrdat[0][0].spline_baseline.baseline_points = [11.787481124553866, -2.2087244721275536]
        nd.nmrdat[0][0].add_baseline_points()
        nd.nmrdat[0][0].corr_spline_baseline()
        self.assertAlmostEqual(nd.nmrdat[0][0].spc[0][51].real / np.max(nd.nmrdat[0][0].spc[0].real), 0.0, places=3)
        self.assertAlmostEqual(nd.nmrdat[0][0].spc[0][65485].real / np.max(nd.nmrdat[0][0].spc[0].real), 0.0, places=3)
        # end test_corr_spline_baseline
    
    def test_add_hsqc_peak(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")
        e_name = "5"
        nd = nmrDataSet.NmrDataSet()
        nd.read_nmrpipe_spc(p_name, e_name, "test.dat")
        nd.nmrdat[0][0].proc.ph0 = [153.05550133, -1.55127223, 0.0]
        nd.nmrdat[0][0].proc.ph1 = [165.12883873, 6.06083814, 0.0]
        nd.nmrdat[0][0].ref_shift = [1.314, 22.8972, 0.0]
        nd.nmrdat[0][0].ref_point = [143, 2355, 0]
        nd.nmrdat[0][0].calc_ppm()
        nd.nmrdat[0][0].phase2d(nd.nmrdat[0][0].proc.ph0[0], nd.nmrdat[0][0].proc.ph1[0], 0)
        nd.nmrdat[0][0].phase2d(nd.nmrdat[0][0].proc.ph0[1], nd.nmrdat[0][0].proc.ph1[1], 1)
        nd.nmrdat[0][0].autobaseline2d()
        metabolite_name = 'L-AsparticAcid'
        cur_peak = 2
        nd.nmrdat[0][0].hsqc.read_metabolite_information(metabolite_name)
        nd.nmrdat[0][0].hsqc.set_metabolite_information(metabolite_name, nd.nmrdat[0][0].hsqc.metabolite_information)
        nd.nmrdat[0][0].hsqc.cur_metabolite = metabolite_name
        nd.nmrdat[0][0].hsqc.cur_peak = cur_peak
        nd.nmrdat[0][0].hsqc.set_peak_information()
        x_data = [3.90869]
        y_data = [54.932118]
        nd.nmrdat[0][0].add_hsqc_peak(x_data, y_data)
        self.assertEqual(nd.nmrdat[0][0].hsqc.hsqc_data[metabolite_name].spin_systems[cur_peak - 1]['c13_shifts'][0][0], 54.932118)
        # end
    
    def test_autobaseline2d(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")
        e_name = "5"
        nd = nmrDataSet.NmrDataSet()
        nd.read_nmrpipe_spc(p_name, e_name, "test.dat")
        nd.nmrdat[0][0].proc.ph0 = [153.05550133, -1.55127223, 0.0]
        nd.nmrdat[0][0].proc.ph1 = [165.12883873, 6.06083814, 0.0]
        nd.nmrdat[0][0].ref_shift = [1.314, 22.8972, 0.0]
        nd.nmrdat[0][0].ref_point = [143, 2355, 0]
        nd.nmrdat[0][0].calc_ppm()
        nd.nmrdat[0][0].phase2d(nd.nmrdat[0][0].proc.ph0[0], nd.nmrdat[0][0].proc.ph1[0], 0)
        nd.nmrdat[0][0].phase2d(nd.nmrdat[0][0].proc.ph0[1], nd.nmrdat[0][0].proc.ph1[1], 1)
        nd.nmrdat[0][0].autobaseline2d()
        # end test_autobaseline2d
    
    def test_autopick_hsqc(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")
        e_name = "5"
        nd = nmrDataSet.NmrDataSet()
        nd.read_nmrpipe_spc(p_name, e_name, "test.dat")
        nd.nmrdat[0][0].proc.ph0 = [153.05550133, -1.55127223, 0.0]
        nd.nmrdat[0][0].proc.ph1 = [165.12883873, 6.06083814, 0.0]
        nd.nmrdat[0][0].ref_shift = [1.314, 22.8972, 0.0]
        nd.nmrdat[0][0].ref_point = [143, 2355, 0]
        nd.nmrdat[0][0].calc_ppm()
        nd.nmrdat[0][0].phase2d(nd.nmrdat[0][0].proc.ph0[0], nd.nmrdat[0][0].proc.ph1[0], 0)
        nd.nmrdat[0][0].phase2d(nd.nmrdat[0][0].proc.ph0[1], nd.nmrdat[0][0].proc.ph1[1], 1)
        nd.nmrdat[0][0].autobaseline2d()
        metabolite_name = 'L-AsparticAcid'
        cur_peak = 2
        nd.nmrdat[0][0].hsqc.read_metabolite_information(metabolite_name)
        nd.nmrdat[0][0].hsqc.set_metabolite_information(metabolite_name, nd.nmrdat[0][0].hsqc.metabolite_information)
        nd.nmrdat[0][0].hsqc.cur_metabolite = metabolite_name
        nd.nmrdat[0][0].hsqc.cur_peak = cur_peak
        nd.nmrdat[0][0].hsqc.set_peak_information()
        nd.nmrdat[0][0].autopick_hsqc([metabolite_name])
        self.assertEqual(nd.nmrdat[0][0].hsqc.hsqc_data[metabolite_name].spin_systems[0]['c13_shifts'][0][0], 55.086)
        # end test_autopick_hsqc
    
    def test_autofit_hsqc(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")
        e_name = "5"
        nd = nmrDataSet.NmrDataSet()
        nd.read_nmrpipe_spc(p_name, e_name, "test.dat")
        nd.nmrdat[0][0].proc.ph0 = [153.05550133, -1.55127223, 0.0]
        nd.nmrdat[0][0].proc.ph1 = [165.12883873, 6.06083814, 0.0]
        nd.nmrdat[0][0].ref_shift = [1.314, 22.8972, 0.0]
        nd.nmrdat[0][0].ref_point = [143, 2355, 0]
        nd.nmrdat[0][0].calc_ppm()
        nd.nmrdat[0][0].phase2d(nd.nmrdat[0][0].proc.ph0[0], nd.nmrdat[0][0].proc.ph1[0], 0)
        nd.nmrdat[0][0].phase2d(nd.nmrdat[0][0].proc.ph0[1], nd.nmrdat[0][0].proc.ph1[1], 1)
        nd.nmrdat[0][0].autobaseline2d()
        metabolite_name = 'L-AsparticAcid'
        cur_peak = 2
        nd.nmrdat[0][0].hsqc.read_metabolite_information(metabolite_name)
        nd.nmrdat[0][0].hsqc.set_metabolite_information(metabolite_name, nd.nmrdat[0][0].hsqc.metabolite_information)
        nd.nmrdat[0][0].hsqc.cur_metabolite = metabolite_name
        nd.nmrdat[0][0].hsqc.cur_peak = cur_peak
        nd.nmrdat[0][0].hsqc.set_peak_information()
        nd.nmrdat[0][0].autopick_hsqc([metabolite_name])
        nd.nmrdat[0][0].autofit_hsqc([metabolite_name])
        self.assertAlmostEqual(nd.nmrdat[0][0].hsqc.hsqc_data[metabolite_name].spin_systems[0]['c13_shifts'][0][0]/100.0, 0.5568829020857215, places=1)
        # end test_autofit_hsqc
    def test_sim_hsqc_1d(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")
        #import os
        #import numpy as np
        #from metabolabpy.nmr import nmrDataSet
        #p_name = "nmrData"
        e_name = "5"
        nd = nmrDataSet.NmrDataSet()
        nd.read_nmrpipe_spc(p_name, e_name, "test.dat")
        nd.nmrdat[0][0].proc.ph0 = [153.05550133, -1.55127223, 0.0]
        nd.nmrdat[0][0].proc.ph1 = [165.12883873, 6.06083814, 0.0]
        nd.nmrdat[0][0].ref_shift = [1.314, 22.8972, 0.0]
        nd.nmrdat[0][0].ref_point = [143, 2355, 0]
        nd.nmrdat[0][0].calc_ppm()
        nd.nmrdat[0][0].phase2d(nd.nmrdat[0][0].proc.ph0[0], nd.nmrdat[0][0].proc.ph1[0], 0)
        nd.nmrdat[0][0].phase2d(nd.nmrdat[0][0].proc.ph0[1], nd.nmrdat[0][0].proc.ph1[1], 1)
        nd.nmrdat[0][0].autobaseline2d()
        metabolite_name = 'L-AsparticAcid'
        cur_peak = 2
        nd.nmrdat[0][0].hsqc.read_metabolite_information(metabolite_name)
        nd.nmrdat[0][0].hsqc.set_metabolite_information(metabolite_name, nd.nmrdat[0][0].hsqc.metabolite_information)
        nd.nmrdat[0][0].hsqc.cur_metabolite = metabolite_name
        nd.nmrdat[0][0].hsqc.cur_peak = cur_peak
        nd.nmrdat[0][0].hsqc.set_peak_information()
        nd.nmrdat[0][0].autopick_hsqc([metabolite_name])
        nd.nmrdat[0][0].hsqc.cur_peak = cur_peak
        nd.nmrdat[0][0].sim_hsqc_1d()
        self.assertEqual(len(nd.nmrdat[0][0].hsqc.hsqc_data[metabolite_name].sim_spc[cur_peak - 1]), 16384)
        # end test_sim_hsqc_1d
    
    def test_fit_hsqc_1d(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")
        #import os
        #import numpy as np
        #from metabolabpy.nmr import nmrDataSet
        #p_name = "nmrData"
        e_name = "5"
        nd = nmrDataSet.NmrDataSet()
        nd.read_nmrpipe_spc(p_name, e_name, "test.dat")
        nd.nmrdat[0][0].proc.ph0 = [153.05550133, -1.55127223, 0.0]
        nd.nmrdat[0][0].proc.ph1 = [165.12883873, 6.06083814, 0.0]
        nd.nmrdat[0][0].ref_shift = [1.314, 22.8972, 0.0]
        nd.nmrdat[0][0].ref_point = [143, 2355, 0]
        nd.nmrdat[0][0].calc_ppm()
        nd.nmrdat[0][0].phase2d(nd.nmrdat[0][0].proc.ph0[0], nd.nmrdat[0][0].proc.ph1[0], 0)
        nd.nmrdat[0][0].phase2d(nd.nmrdat[0][0].proc.ph0[1], nd.nmrdat[0][0].proc.ph1[1], 1)
        nd.nmrdat[0][0].autobaseline2d()
        metabolite_name = 'L-AsparticAcid'
        cur_peak = 2
        nd.nmrdat[0][0].hsqc.read_metabolite_information(metabolite_name)
        nd.nmrdat[0][0].hsqc.set_metabolite_information(metabolite_name, nd.nmrdat[0][0].hsqc.metabolite_information)
        nd.nmrdat[0][0].hsqc.cur_metabolite = metabolite_name
        nd.nmrdat[0][0].hsqc.cur_peak = cur_peak
        nd.nmrdat[0][0].hsqc.set_peak_information()
        nd.nmrdat[0][0].autopick_hsqc([metabolite_name])
        nd.nmrdat[0][0].hsqc.cur_peak = cur_peak
        nd.nmrdat[0][0].fit_hsqc_1d()
        self.assertAlmostEqual(nd.nmrdat[0][0].hsqc.hsqc_data[metabolite_name].spin_systems[cur_peak - 1]['contribution'][0], 28.973, places=1)
        # end test_fit_hsqc_1d

if __name__ == "__main__":
    unittest.main()
