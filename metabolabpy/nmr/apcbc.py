'''
Auomatic phase and baseline correction for 1D NMR spectra

Using the histogram method
'''

from metabolabpy.nmr import nmrData
import numpy as np
#import matplotlib.pyplot as pl
from scipy.optimize import fmin
import math
from scipy.fftpack import fft, ifft, fftshift


class Apcbc:

    def __init__(self):
        self.correctBaseline = 0
        self.npts = 0
        self.nMax = 10.00
        self.nOrder = 1
        self.nbins = 5001
        self.mFact0 = 20 * 10 ** (self.nOrder + 1)
        self.mFact1 = 5 * 10 ** (self.nOrder + 1)
        self.startPts = 3500
        self.endPts = 4000
        self.pars = np.array([])
        self.rSpc = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self.iSpc = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    def baselineFitFunc(self, par, spc, xaxis):
        l = self.calcLine(spc, xaxis)
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
        # end baselineFitFunc

    def baselineFitFuncEval(self, par, spc, xaxis):  # , plot):
        l = self.calcLine(spc, xaxis)
        spc = spc - l
        nc = int((len(par)) / 2)
        cr = par[:nc]
        ci = par[nc:]
        blr = np.polynomial.chebyshev.chebval(xaxis, cr)
        bli = np.polynomial.chebyshev.chebval(xaxis, ci)
        spc = spc - blr - 1j * bli
        # if(plot == True):
        #    hr  = self.hist(spc.real)
        #    hi  = self.hist(spc.imag)
        #    xr  = hr[1][1:nbins+1]
        #    xi  = hi[1][1:nbins+1]
        #    pl.subplot(2,1,1)
        #    pl.plot(xr, hr[0])
        #    pl.subplot(2,1,2)
        #    pl.plot(xi, hi[0])
        #    pl.show()

        return spc
        # end baselineFitFuncEval

    def calcLine(self, spc, xaxis):
        x1 = np.mean(xaxis[self.startPts:self.endPts])
        x2 = np.mean(xaxis[-self.endPts:-self.startPts])
        y1r = np.mean(spc.real[self.startPts:self.endPts])
        y2r = np.mean(spc.real[-self.endPts:-self.startPts])
        y1i = np.mean(spc.imag[self.startPts:self.endPts])
        y2i = np.mean(spc.imag[-self.endPts:-self.startPts])
        slopeR = (y2r - y1r) / (x2 - x1)
        slopeI = (y2i - y1i) / (x2 - x1)
        lineR = slopeR * (xaxis - x1) + y1r
        lineI = slopeI * (xaxis - x1) + y1i
        line = lineR + 1j * lineI
        return line
        # end calcLine

    def fitBaseline(self, spc, xaxis):
        npts = len(spc)
        l = self.calcLine(spc, xaxis)
        nCoeffReal = np.polynomial.chebyshev.chebfit(xaxis, spc.real - l.real, self.nOrder)
        nCoeffImag = np.polynomial.chebyshev.chebfit(xaxis, spc.imag - l.imag, self.nOrder)
        pars = np.array([])
        pars = np.append(pars, nCoeffReal)
        pars = np.append(pars, nCoeffImag)
        print("...Fitting phase and baseline...")
        plsq = fmin(self.baselineFitFunc, pars, args=(spc, xaxis), ftol=1e-3, xtol=1e-3, maxfun=1000)
        print("...finished")
        rSpcPar = plsq[:int(len(plsq) / 2)]
        iSpcPar = plsq[int(len(plsq) / 2):]
        self.rSpc[:len(rSpcPar)] = rSpcPar
        self.iSpc[:len(iSpcPar)] = iSpcPar
        return plsq
        # end fitBaseline

    def fitPhase(self, spc, xaxis):
        npts = len(spc)
        ph0 = 0.0
        ph1 = 0.0
        l = self.calcLine(spc, xaxis)
        nCoeffReal = np.polynomial.chebyshev.chebfit(xaxis, spc.real - l.real, self.nOrder)
        nCoeffImag = np.polynomial.chebyshev.chebfit(xaxis, spc.imag - l.imag, self.nOrder)
        pars = np.array([])
        pars = np.append(pars, ph0)
        pars = np.append(pars, ph1)
        pars = np.append(pars, nCoeffReal)
        pars = np.append(pars, nCoeffImag)
        print("...Fitting phase and baseline...")
        plsq = fmin(self.phaseFitFunc, pars, args=(spc, xaxis), ftol=1e-3, xtol=1e-3, maxfun=1000)
        print("...finished")
        rSpcPar = plsq[2:int(len(plsq) / 2 + 2)]
        iSpcPar = plsq[int(len(plsq) / 2 + 2):]
        self.rSpc[:len(rSpcPar)] = rSpcPar
        self.iSpc[:len(iSpcPar)] = iSpcPar
        return plsq
        # end fitPhase

    #def getHist(self, spc, xaxis, plot=True, ph0=0.0, ph1=0.0):
    def getHist(self, spc, xaxis, ph0=0.0, ph1=0.0):
        hr = self.hist(spc.real)
        hi = self.hist(spc.imag)
        xr = hr[1][1:self.nbins + 1]
        xi = hi[1][1:self.nbins + 1]
        w = self.histRWeight(xr)
        hr2 = hr[0] * w
        hi2 = hi[0]
        pp = False
        if (ph0 != 0.0) or (ph1 != 0):
            pp = True

        if pp == True:
            spc2a = phase2(spc, ph0, ph1)
            hra = self.hist(spc2a.real)
            hia = self.hist(spc2a.imag)
            xra = hra[1][1:self.nbins + 1]
            xia = hia[1][1:self.nbins + 1]
            wa = self.histRWeight(xra)
            hr2a = hra[0] * wa
            hi2a = hia[0]
            maxPos = np.where(hi2a == np.max(hi2a))
            if (np.abs(maxPos[0][0] - len(hi2a) / 2) > 2):
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
        # end getHist

    def histRWeight(self, x):
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
        # end histRWeight

    def histIWeight(self, x):
        weight = -x * x
        return weight
        # end histIWeight

    def hist(self, spc):
        h = np.histogram(spc, bins=self.nbins, range=(-1.0, 1.0))
        return h

    def phase2(elf, mat, ph0, ph1):
        npts = len(mat)
        ph0 = -ph0 * math.pi / 180.0
        ph1 = -ph1 * math.pi / 180.0
        t = complex()
        frac = np.linspace(0, 1, npts)
        ph = ph0 + frac * ph1
        mat = np.cos(ph) * mat.real + np.sin(ph) * mat.imag + 1j * (-np.sin(ph) * mat.real + np.cos(ph) * mat.imag)
        return mat
        # end phase2

    def phaseFitFunc(self, par, spc, xaxis):
        ph0 = self.mFact0 * par[0]
        ph1 = self.mFact1 * par[1]
        spc = self.phase2(spc, ph0, ph1)
        l = self.calcLine(spc, xaxis)
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
        # print("%4.2f & %4.2f | %4.8f & %4.8f" % (n1, n2, ph0, ph1))
        return 1.0 * n1 + 100.0 * n2
        # end phaseFitFunc

    def phaseFitFuncEval(self, par, spc, xaxis): #, plot):
        ph0 = self.mFact0 * par[0]
        ph1 = self.mFact1 * par[1]
        spc = self.phase2(spc, ph0, ph1)
        l = self.calcLine(spc, xaxis)
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
        # end phaseFitFuncEval

    def qr(self, histReal):
        x = histReal[1][1:self.nbins + 1]
        w = self.histRWeight(x)
        h = histReal[0] * w
        return np.sum(h)
        # end qr

    def qi(self, histImag):
        x = histImag[1][1:self.nbins + 1]
        w = self.histIWeight(x)
        h = histImag[0] * w
        return np.sum(h)
        # end qi

    def setVars(self, spc):
        self.npts = len(spc)
        # end setVars
