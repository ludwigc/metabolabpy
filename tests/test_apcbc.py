#!/usr/bin/env python


"""

test_acqPars

author(s): C. Ludwig
origin: 02-08-2023


"""


import unittest
from metabolabpy.nmr import nmrData
from metabolabpy.nmr import apcbc
import os
import numpy as np


class acqParsTestCase(unittest.TestCase):

    def test_baseline_fit_func(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        ab = apcbc.Apcbc()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.auto_ref()
        npts = len(nd.spc[0])
        l = ab.calc_line(nd.spc[0], nd.ppm1)
        n_coeff_real = np.polynomial.chebyshev.chebfit(nd.ppm1, nd.spc[0].real - l.real, ab.n_order)
        n_coeff_imag = np.polynomial.chebyshev.chebfit(nd.ppm1, nd.spc[0].imag - l.imag, ab.n_order)
        pars = np.array([])
        pars = np.append(pars, n_coeff_real)
        pars = np.append(pars, n_coeff_imag)
        baseline_func = ab.baseline_fit_func(pars, nd.spc[0], nd.ppm1)
        self.assertEqual(baseline_func, 0.0)

    def test_baseline_fit_func_eval(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        ab = apcbc.Apcbc()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.auto_ref()
        npts = len(nd.spc[0])
        l = ab.calc_line(nd.spc[0], nd.ppm1)
        n_coeff_real = np.polynomial.chebyshev.chebfit(nd.ppm1, nd.spc[0].real - l.real, ab.n_order)
        n_coeff_imag = np.polynomial.chebyshev.chebfit(nd.ppm1, nd.spc[0].imag - l.imag, ab.n_order)
        pars = np.array([])
        pars = np.append(pars, n_coeff_real)
        pars = np.append(pars, n_coeff_imag)
        baseline_func_eval = ab.baseline_fit_func_eval(pars, nd.spc[0], nd.ppm1)
        self.assertAlmostEqual(baseline_func_eval[0].real, -371640.24939898413)

    def test_calc_line(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        ab = apcbc.Apcbc()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.auto_ref()
        line = ab.calc_line(nd.spc[0], nd.ppm1)
        self.assertAlmostEqual(line[0].real, -5281.861632170241)

    def test_fit_baseline(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        ab = apcbc.Apcbc()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.auto_ref()
        plsq = ab.fit_baseline(nd.spc[0], nd.ppm1)
        self.assertEqual(round(plsq[0]), 1659853)

    def test_fit_phase(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        ab = apcbc.Apcbc()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.auto_ref()
        plsq = ab.fit_phase(nd.spc[0], nd.ppm1)
        self.assertAlmostEqual(plsq[0], -5.948091860368079e-11)

    def test_get_hist(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        ab = apcbc.Apcbc()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.auto_ref()
        hist = ab.get_hist(nd.spc[0], nd.ppm1)
        self.assertAlmostEqual(hist[0].real, -1.499368547343163)

    def test_hist(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        ab = apcbc.Apcbc()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.auto_ref()
        hist = ab.hist(nd.spc[0].real)
        self.assertEqual(hist[0][0], 1)

    def test_phase2(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        ab = apcbc.Apcbc()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.auto_ref()
        mat = ab.phase2(nd.spc[0], 0.0, 0.0)
        self.assertEqual(mat[0].real, nd.spc[0][0].real)

    def test_phase_fit_func(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        ab = apcbc.Apcbc()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.auto_ref()
        l = ab.calc_line(nd.spc[0], nd.ppm1)
        n_coeff_real = np.polynomial.chebyshev.chebfit(nd.ppm1, nd.spc[0].real - l.real, ab.n_order)
        n_coeff_imag = np.polynomial.chebyshev.chebfit(nd.ppm1, nd.spc[0].imag - l.imag, ab.n_order)
        ph0 = 0.0
        ph1 = 0.0
        pars = np.array([])
        pars = np.append(pars, ph0)
        pars = np.append(pars, ph1)
        pars = np.append(pars, n_coeff_real)
        pars = np.append(pars, n_coeff_imag)
        phase_fit_func = ab.phase_fit_func(pars, nd.spc[0], nd.ppm1)
        self.assertEqual(phase_fit_func, 0.0)

    def test_phase_fit_func_eval(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrData.NmrData()
        ab = apcbc.Apcbc()
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()
        nd.auto_ref()
        l = ab.calc_line(nd.spc[0], nd.ppm1)
        n_coeff_real = np.polynomial.chebyshev.chebfit(nd.ppm1, nd.spc[0].real - l.real, ab.n_order)
        n_coeff_imag = np.polynomial.chebyshev.chebfit(nd.ppm1, nd.spc[0].imag - l.imag, ab.n_order)
        ph0 = 0.0
        ph1 = 0.0
        pars = np.array([])
        pars = np.append(pars, ph0)
        pars = np.append(pars, ph1)
        pars = np.append(pars, n_coeff_real)
        pars = np.append(pars, n_coeff_imag)
        phase_fit_func_eval = ab.phase_fit_func_eval(pars, nd.spc[0], nd.ppm1)
        self.assertAlmostEqual(phase_fit_func_eval[0].real, -371640.24939898413)


if __name__ == "__main__":
    unittest.main()