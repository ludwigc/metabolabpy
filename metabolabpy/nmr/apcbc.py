"""
Auomatic phase and baseline correction for 1D NMR spectra

Using the histogram method
"""

from metabolabpy.nmr import nmrData
import numpy as np
# import matplotlib.pyplot as pl
from scipy.optimize import fmin
import math
from scipy.fftpack import fft, ifft, fftshift


class Apcbc:

    def __init__(self):
        self.correct_baseline = 0
        self.npts = 0
        self.n_max = 10.00
        self.n_order = 1
        self.nbins = 5001
        self.m_fact0 = 20 * 10 ** (self.n_order + 1)
        self.m_fact1 = 5 * 10 ** (self.n_order + 1)
        self.start_pts = 3500
        self.end_pts = 4000
        self.pars = np.array([])
        self.r_spc = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self.i_spc = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    def baseline_fit_func(self, par, spc, xaxis):
        l = self.calc_line(spc, xaxis)
        spc = spc - l
        nc = int((len(par)) / 2)
        cr = par[:nc]
        ci = par[nc:]
        blr = np.polynomial.chebyshev.chebval(xaxis, cr)
        bli = np.polynomial.chebyshev.chebval(xaxis, ci)
        spc = spc - blr - 1j * bli
        hr = self.hist(spc.real)
        hi = self.hist(spc.imag)
        n1 = -self.qr(hr)
        n2 = self.qi(hi)
        # print("%4.2f & %4.2f" % (n1, n2))
        return 1.0 * n1 + 100.0 * n2
        # end baseline_fit_func

    def baseline_fit_func_eval(self, par, spc, xaxis):  # , plot):
        l = self.calc_line(spc, xaxis)
        spc = spc - l
        nc = int((len(par)) / 2)
        cr = par[:nc]
        ci = par[nc:]
        blr = np.polynomial.chebyshev.chebval(xaxis, cr)
        bli = np.polynomial.chebyshev.chebval(xaxis, ci)
        spc = spc - blr - 1j * bli
        return spc
        # end baseline_fit_func_eval

    def calc_line(self, spc, xaxis):
        x1 = np.mean(xaxis[self.start_pts:self.end_pts])
        x2 = np.mean(xaxis[-(self.end_pts + 1):-(self.start_pts + 1)])
        y1r = np.mean(spc.real[self.start_pts:self.end_pts])
        y2r = np.mean(spc.real[-(self.end_pts + 1):-(self.start_pts + 1)])
        y1i = np.mean(spc.imag[self.start_pts:self.end_pts])
        y2i = np.mean(spc.imag[-(self.end_pts + 1):-(self.start_pts + 1)])
        slope_r = (y2r - y1r) / (x2 - x1)
        slope_i = (y2i - y1i) / (x2 - x1)
        line_r = slope_r * (xaxis - x1) + y1r
        line_i = slope_i * (xaxis - x1) + y1i
        line = line_r + 1j * line_i
        return line
        # end calc_line

    def fit_baseline(self, spc, xaxis):
        npts = len(spc)
        l = self.calc_line(spc, xaxis)
        n_coeff_real = np.polynomial.chebyshev.chebfit(xaxis, spc.real - l.real, self.n_order)
        n_coeff_imag = np.polynomial.chebyshev.chebfit(xaxis, spc.imag - l.imag, self.n_order)
        pars = np.array([])
        pars = np.append(pars, n_coeff_real)
        pars = np.append(pars, n_coeff_imag)
        print("...Fitting phase and baseline...")
        plsq = fmin(self.baseline_fit_func, pars, args=(spc, xaxis), ftol=1e-3, xtol=1e-3, maxfun=1000)
        print("...finished")
        r_spc_par = plsq[:int(len(plsq) / 2)]
        i_spc_par = plsq[int(len(plsq) / 2):]
        self.r_spc[:len(r_spc_par)] = r_spc_par
        self.i_spc[:len(i_spc_par)] = i_spc_par
        return plsq
        # end fit_baseline

    def fit_phase(self, spc, xaxis):
        npts = len(spc)
        ph0 = 0.0
        ph1 = 0.0
        l = self.calc_line(spc, xaxis)
        n_coeff_real = np.polynomial.chebyshev.chebfit(xaxis, spc.real - l.real, self.n_order)
        n_coeff_imag = np.polynomial.chebyshev.chebfit(xaxis, spc.imag - l.imag, self.n_order)
        pars = np.array([])
        pars = np.append(pars, ph0)
        pars = np.append(pars, ph1)
        pars = np.append(pars, n_coeff_real)
        pars = np.append(pars, n_coeff_imag)
        print("...Fitting phase and baseline...")
        plsq = fmin(self.phase_fit_func, pars, args=(spc, xaxis), ftol=1e-3, xtol=1e-3, maxfun=1000)
        print("...finished")
        r_spc_par = plsq[2:int(len(plsq) / 2 + 2)]
        i_spc_par = plsq[int(len(plsq) / 2 + 2):]
        self.r_spc[:len(r_spc_par)] = r_spc_par
        self.i_spc[:len(i_spc_par)] = i_spc_par
        return plsq
        # end fit_phase

    #def get_hist(self, spc, xaxis, plot=True, ph0=0.0, ph1=0.0):
    def get_hist(self, spc, xaxis, ph0=0.0, ph1=0.0):
        hr = self.hist(spc.real)
        hi = self.hist(spc.imag)
        xr = hr[1][1:self.nbins + 1]
        xi = hi[1][1:self.nbins + 1]
        w = self.hist_r_weight(xr)
        hr2 = hr[0] * w
        hi2 = hi[0]
        pp = False
        if (ph0 != 0.0) or (ph1 != 0):
            pp = True

        if pp == True:
            spc2a = self.phase2(spc, ph0, ph1)
            hra = self.hist(spc2a.real)
            hia = self.hist(spc2a.imag)
            xra = hra[1][1:self.nbins + 1]
            xia = hia[1][1:self.nbins + 1]
            wa = self.hist_r_weight(xra)
            hr2a = hra[0] * wa
            hi2a = hia[0]
            max_pos = np.where(hi2a == np.max(hi2a))
            if (np.abs(max_pos[0][0] - len(hi2a) / 2) > 2):
                hi2a = hi2a * np.sign(xia)

        #if (plot == True):
        #    pl.subplot(2, 1, 1)
        #    if (pp == True):
        #        pl.plot(xra, hra[0])
        #        pl.plot(xra, hr2a)
        #
        #    pl.plot(xr, hr[0])
        #    pl.plot(xr, hr2)
        #    pl.subplot(2, 1, 2)
        #    if (pp == True):
        #        pl.plot(xia, hi2a)
        #
        #    pl.plot(xi, hi2)
        #    pl.show()

        return hr2 + 1j * hi2
        # end get_hist

    def hist_r_weight(self, x):
        n1 = x[np.where(x < 0.0)]
        n2 = x[np.where(x < 0.05)]
        n2 = n2[np.where(n2 >= 0.0)]
        n3 = x[np.where(x < 0.5)]
        n3 = n3[np.where(n3 >= 0.05)]
        y4 = x[np.where(x >= 0.5)]
        y1 = (n1 + 0.05) * 1.5 / (1 - 0.05)
        y2 = 1 - np.abs(n2) / 0.05
        y3 = -(n3 - 0.05) * 0.5 / (0.5 - 0.05)
        y4 = -0.5 * np.ones(len(y4))
        weight = np.append(y1, y2)
        weight = np.append(weight, y3)
        weight = np.append(weight, y4)
        return weight
        # end hist_r_weight

    def hist_i_weight(self, x):
        weight = -x * x
        return weight
        # end hist_i_weight

    def hist(self, spc):
        h = np.histogram(spc, bins=self.nbins, range=(-1.0, 1.0))
        return h

    def phase2(self, mat, ph0, ph1):
        npts = len(mat)
        ph0 = -ph0 * math.pi / 180.0
        ph1 = -ph1 * math.pi / 180.0
        t = complex()
        frac = np.linspace(0, 1, npts)
        ph = ph0 + frac * ph1
        mat = np.cos(ph) * mat.real + np.sin(ph) * mat.imag + 1j * (-np.sin(ph) * mat.real + np.cos(ph) * mat.imag)
        return mat
        # end phase2

    def phase_fit_func(self, par, spc, xaxis):
        ph0 = self.m_fact0 * par[0]
        ph1 = self.m_fact1 * par[1]
        spc = self.phase2(spc, ph0, ph1)
        l = self.calc_line(spc, xaxis)
        spc = spc - l
        nc = int((len(par) - 2) / 2)
        cr = par[2:2 + nc]
        ci = par[2 + nc:2 + 2 * nc]
        blr = np.polynomial.chebyshev.chebval(xaxis, cr)
        bli = np.polynomial.chebyshev.chebval(xaxis, ci)
        spc = spc - blr - 1j * bli
        hr = self.hist(spc.real)
        hi = self.hist(spc.imag)
        n1 = -self.qr(hr)
        n2 = self.qi(hi)
        return 1.0 * n1 + 100.0 * n2
        # end phase_fit_func

    def phase_fit_func_eval(self, par, spc, xaxis): #, plot):
        ph0 = self.m_fact0 * par[0]
        ph1 = self.m_fact1 * par[1]
        spc = self.phase2(spc, ph0, ph1)
        l = self.calc_line(spc, xaxis)
        spc = spc - l
        nc = int((len(par) - 2) / 2)
        cr = par[2:2 + nc]
        ci = par[2 + nc:2 + 2 * nc]
        blr = np.polynomial.chebyshev.chebval(xaxis, cr)
        bli = np.polynomial.chebyshev.chebval(xaxis, ci)
        spc = spc - blr - 1j * bli
        #if (plot == True):
        #    hr = self.hist(spc.real)
        #    hi = self.hist(spc.imag)
        #    xr = hr[1][1:nbins + 1]
        #    xi = hi[1][1:nbins + 1]
        #    pl.subplot(2, 1, 1)
        #    pl.plot(xr, hr[0])
        #    pl.subplot(2, 1, 2)
        #    pl.plot(xi, hi[0])
        #    pl.show()

        return spc
        # end phase_fit_func_eval

    def qr(self, hist_real):
        x = hist_real[1][1:self.nbins + 1]
        w = self.hist_r_weight(x)
        h = hist_real[0] * w
        return np.sum(h)
        # end qr

    def qi(self, hist_imag):
        x = hist_imag[1][1:self.nbins + 1]
        w = self.hist_i_weight(x)
        h = hist_imag[0] * w
        return np.sum(h)
        # end qi

    def set_vars(self, spc):
        self.npts = len(spc)
        # end set_vars
