import numpy as np
# from metabolabpy.nmr import wdwf
from metabolabpy.nmr import acqPars
from metabolabpy.nmr import procPars
from metabolabpy.nmr import dispPars
import os
from scipy.fftpack import fft, ifft, fftshift
import math
import matplotlib.pyplot as pl
from scipy import signal
import time
from metabolabpy.nmr import nmrpipeData
from scipy.optimize import leastsq
from metabolabpy.nmr import apcbc


class NmrData:

    def __init__(self):
        self.fid = np.array([[]], dtype='complex')
        self.spc = np.array([[]], dtype='complex')
        self.ppm1 = np.array([], dtype='float64')
        self.ppm2 = np.array([], dtype='float64')
        self.ppm3 = np.array([], dtype='float64')
        self.dim = 0
        self.title = np.array([], dtype='str')
        self.origDataSet = str('')
        self.phCorrMode = 0
        self.acq = acqPars.AcqPars()
        self.proc = procPars.ProcPars()
        self.display = dispPars.DispPars()
        self.fidOffsetCorr = 0
        self.dataSetName = ''
        self.dataSetNumber = ''
        self.title = str('Empty NMR data set')
        self.pulseProgram = str('')
        self.windowFunction = {'none': 0, 'exponential': 1, 'gaussian': 2, 'sine': 3, 'qsine': 4, 'sem': 5}
        self.refShift = np.array([0, 0, 0], dtype='float64')
        self.refPoint = np.array([0, 0, 0], dtype='int')
        self.refsw = np.array([0, 0, 0], dtype='float64')
        self.refTmspRange = 0.3  # [ppm]
        self.apc = apcbc.Apcbc()
        self.projectedJres = False
        self.origJresSet = -1
        self.origJresExp = -1
        self.pjresMode = ''
        self.jRes = False
        self.acqu = str('')
        self.acqu2 = str('')
        self.acqu3 = str('')
        self.audita_txt = str('')
        self.format_temp = str('')
        self.fq1list = str('')
        self.scon2 = str('')
        self.shimvalues = str('')
        self.uxnmr_par = str('')
        self.proc1 = str('')
        self.proc2 = str('')
        self.proc3 = str('')
        self.outd = str('')
        # end __init__

    def __str__(self):  # pragma: no cover
        rString = 'MetaboLabPy NMR Data (v. 0.1)\n'
        rString += '__________________________________________________________________\n'
        rString += self.title
        return rString
        # end __str__

    def apodise(self, fid, dim, lb, gb, ssb, groupDelay, sw_h):
        fid = np.copy(fid)
        if (self.proc.windowType[dim] == 0):  # no window
            wdwf = np.ones(len(fid))

        if (self.proc.windowType[dim] == 1):  # exponential window
            t = (np.linspace(0.0, len(fid) - 1, len(fid)) - groupDelay) / sw_h
            wdwf = np.exp(-lb * t)

        if (self.proc.windowType[dim] == 2):  # gaussian window
            t = (np.linspace(0.0, len(fid) - 1 - groupDelay, len(fid))) / sw_h
            wdwf = np.exp(-lb * 2 * math.pi * t - 2 * math.pi * (t ** 2) / ((2 * math.pi * gb * len(fid) / sw_h)))

        if (self.proc.windowType[dim] == 3):  # sine window
            if self.acq.fnMode == 1 or self.acq.fnMode == 2:
                npts = int(min(self.acq.nDataPoints[dim], len(fid)))
            else:
                npts = int(min(self.acq.nDataPoints[dim] / 2, len(fid)))

            t = (np.linspace(0.0, npts - 1, npts)) / (npts - 1)
            wdwf = np.zeros(len(fid))
            ssb2 = ssb * math.pi / 180.0
            wdwf[:npts] = np.sin(ssb2 + (math.pi - ssb2) * t)

        if (self.proc.windowType[dim] == 4):  # qsine window
            if self.acq.fnMode == 1 or self.acq.fnMode == 2:
                npts = int(min(self.acq.nDataPoints[dim], len(fid)))
            else:
                npts = int(min(self.acq.nDataPoints[dim] / 2, len(fid)))

            t = (np.linspace(0.0, npts - 1, npts)) / (npts - 1)
            wdwf = np.zeros(len(fid))
            ssb2 = ssb * math.pi / 180.0
            wdwf[:npts] = np.sin(ssb2 + (math.pi - ssb2) * t) ** 2

        if (self.proc.windowType[dim] == 5):  # sem window
            if self.acq.fnMode == 1 or self.acq.fnMode == 2:
                npts = int(min(self.acq.nDataPoints[dim], len(fid)))
            else:
                npts = int(min(self.acq.nDataPoints[dim] / 2, len(fid)))

            t1 = (np.linspace(0.0, npts - 1 - groupDelay, npts)) / sw_h
            t2 = (np.linspace(0.0, npts - 1, npts)) / (npts - 1)
            wdwf = np.zeros(len(fid))
            ssb2 = ssb * math.pi / 180.0
            wdwf[:npts] = np.exp(-lb * t1) * np.sin(ssb2 + (math.pi - ssb2) * t2)

        fid = np.copy(wdwf * fid)
        return fid
        # end apodise

    def autobaseline1d(self):
        spc = self.spc[0]
        self.apc.npts = len(spc)
        scaleFact = np.max(np.abs(spc))
        spc /= scaleFact
        xaxis = np.linspace(-self.apc.nMax, self.apc.nMax, self.apc.npts)
        parEval = self.apc.fitBaseline(spc, xaxis)
        spc2 = self.apc.baselineFitFuncEval(parEval, spc, xaxis)  # , False)
        self.apc.pars = parEval
        self.apc.rSpc = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self.apc.iSpc = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self.apc.rSpc[:int(len(parEval) / 2)] = parEval[:int(len(parEval) / 2)]
        self.apc.iSpc[:int(len(parEval) / 2)] = parEval[int(len(parEval) / 2):]
        self.spc[0] = spc2 * scaleFact
        self.apc.correctBaseline = 1
        self.procSpc1D()
        self.baseline1d()
        # end autobaseline1d

    def autophase1d(self):
        spc = self.spc[0]
        self.apc.npts = len(spc)
        scaleFact = np.max(np.abs(spc))
        spc /= scaleFact
        xaxis = np.linspace(-self.apc.nMax, self.apc.nMax, self.apc.npts)
        parEval = self.apc.fitPhase(spc, xaxis)
        spc2 = self.apc.phaseFitFuncEval(parEval, spc, xaxis)  # , False)
        if (np.min(spc2.real) == -np.max(np.abs(spc2.real))):
            parEval[0] = ((parEval[0] * self.apc.mFact0 + 180.0) % 360) / self.apc.mFact0

        parEval[0] = (((parEval[0] * self.apc.mFact0 + 180.0) % 360) - 180.0) / self.apc.mFact0
        spc2 = self.apc.phaseFitFuncEval(parEval, spc, xaxis)  # , False)
        self.apc.pars = parEval
        self.apc.rSpc = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self.apc.iSpc = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self.apc.rSpc[:int((len(parEval) - 2) / 2)] = parEval[2:2 + int((len(parEval) - 2) / 2)]
        self.apc.iSpc[:int((len(parEval) - 2) / 2)] = parEval[2 + int((len(parEval) - 2) / 2):]
        self.spc[0] = spc2 * scaleFact
        self.proc.ph0[0] += parEval[0] * self.apc.mFact0
        self.proc.ph1[0] += parEval[1] * self.apc.mFact1
        self.apc.correctBaseline = 1
        self.autoRef()
        self.procSpc1D()
        self.baseline1d()
        # end autophase1d

    def autoRef(self, tmsp=True):
        if self.acq.o1 == 0:
            self.refShift[0] = 4.76
        else:
            self.refShift[0] = self.acq.o1 / self.acq.bf1

        self.refPoint[0] = int(len(self.spc[0]) / 2)
        if (self.dim == 2):
            self.refShift[1] = (self.acq.spcFrequency[1] + self.acq.spcOffset[1]) / self.acq.spcFrequency[1] - 1.0
            self.refPoint[1] = int(len(self.spc) / 2)
            if (tmsp == True):
                self.refPoint[0] = self.ppm2points(0.0, 0)
                self.refShift[0] = 0.0
                pts = self.ppm2points(np.array([-self.refTmspRange, self.refTmspRange]), 0)
                npts = len(self.spc[0])
                r = np.arange(npts - max(pts), npts - min(pts))
                spc = np.sum(self.spc, 0)
                spc = spc[r].real
                refP = np.where(spc == np.amax(spc))
                self.refPoint[0] -= refP[0][0] - int((max(pts) - min(pts)) / 2) + 1

            self.proc.refPoint[0] = self.refPoint[0]
            self.proc.refPoint[1] = self.refPoint[1]

        if (self.dim == 1):
            if (tmsp == True):
                self.refPoint[0] = self.ppm2points(0.0, 0)
                self.refShift[0] = 0.0
                pts = self.ppm2points(np.array([-self.refTmspRange, self.refTmspRange]), 0)
                npts = len(self.spc[0])
                r = np.arange(npts - max(pts), npts - min(pts))
                r = np.copy(r[np.where(r < len(self.spc[0]))])
                spc = self.spc[0][r].real
                refP = np.where(spc == np.amax(spc))
                self.refPoint[0] -= refP[0][0] - int((max(pts) - min(pts)) / 2) + 1

        self.calcPPM()

        # end autoRef

    def baseline1d(self):
        spc = self.spc[0]
        self.apc.npts = len(spc)
        # scaleFact         = np.max(np.abs(spc))
        # spc              /= scaleFact
        xaxis = np.linspace(-self.apc.nMax, self.apc.nMax, self.apc.npts)
        parEval = np.copy(self.apc.rSpc)
        parEval = np.append(parEval, self.apc.iSpc)
        # print(parEval)
        spc2 = self.apc.baselineFitFuncEval(parEval, spc, xaxis)  # , False)
        self.spc[0] = spc2  # *scaleFact
        # end baseline1d

    def calcPPM(self):
        if self.proc.stripStart == 0:
            stsr = 1
        else:
            stsr = self.proc.stripStart

        if self.proc.stripEnd > stsr:
            stsi = self.proc.stripEnd
        else:
            stsi = int(len(self.spc[0]))

        if (self.display.axisType1 == 'ppm'):
            self.ppm1 = self.points2ppm(np.linspace(self.proc.nPoints[0] - 1, 0, self.proc.nPoints[0]), 0)
        else:
            self.ppm1 = self.points2Hz(np.linspace(self.proc.nPoints[0] - 1, 0, self.proc.nPoints[0]), 0)

        self.ppm1 = self.ppm1[stsr-1:stsi]

        if (self.dim > 1):
            npts = int(len(self.spc))
            if (self.display.axisType2 == 'ppm'):
                self.ppm2 = self.points2ppm(np.linspace(npts - 1, 0, npts), 1)
            else:
                self.ppm2 = self.points2Hz(np.linspace(npts - 1, 0, npts), 1)

        # end calcPPM

    def conv(self, fid):
        f1 = np.copy(fid)
        x = np.linspace(0, len(fid) - 1, len(fid))
        f0 = np.copy(fid)
        filtFid = fid[list(range(round(self.acq.groupDelay)))]
        fid = fid[list(range(round(self.acq.groupDelay), len(fid)))]
        ws2 = int(self.proc.convWindowSize[0] / 2)
        es = self.proc.convExtrapolationSize[0]
        wt = int(self.proc.convWindowType[0])
        win = np.zeros(int(ws2 * 2 + 1))
        sWin = 0
        for k in range(int(ws2 * 2 + 1)):
            if (wt == 0):  # Gaussian window
                win[k] = math.exp(-4.0 * ((k - ws2) ** 2) / (ws2 ** 2))
                sWin += win[k]
            else:  # sine bell window
                win[k] = math.cos((math.pi * (k - ws2)) / (2 * ws2 + 2))
                sWin += win[k]

        fid2 = np.convolve(win, fid) / sWin
        # Extrapolation of first sw2 data points
        idx = np.linspace(2 * ws2 - 1, 0, 2 * ws2, dtype='int')
        idx2 = np.linspace(0, 2 * ws2 - 1, 2 * ws2, dtype='int')
        fid2[idx2].real = fid2[int(2 * ws2)].real + np.mean(np.diff(fid2[idx + int(2 * ws2) - 1].real)) * idx
        fid2[idx2].imag = fid2[int(2 * ws2)].imag + np.mean(np.diff(fid2[idx + int(2 * ws2) - 1].imag)) * idx
        # Extrapolation of last sw2 data points
        idx = np.linspace(0, 2 * ws2 - 1, 2 * ws2, dtype='int')
        idx2 = np.linspace(len(fid) - int(2 * ws2) - 1, len(fid) - 1, int(2 * ws2), dtype='int')
        fid2[idx2].real = fid2[len(fid) - int(2 * ws2) - 1].real - np.mean(
            np.diff(fid2[len(fid) - int(2 * ws2) - 1 - idx].real)) * idx
        fid2[idx2].imag = fid2[len(fid) - int(2 * ws2) - 1].imag - np.mean(
            np.diff(fid2[len(fid) - int(2 * ws2) - 1 - idx].imag)) * idx
        fid2 = np.delete(fid2, np.linspace(0, 2 * ws2 - 1, 2 * ws2, dtype='int'))
        fid -= fid2  # (fidRe + 1j*fidIm)
        fid = np.concatenate([filtFid, fid])
        return fid
        # end conv

    def exportBruker1d(self, pathName, expName, scaleFactor = -1):
        if self.acq.manufacturer != 'Bruker':
            return

        if os.path.isdir(pathName + os.sep + expName) is False:
            os.makedirs(pathName + os.sep + expName)

        if os.path.isdir(pathName + os.sep + expName + os.sep + 'pdata' + os.sep + '1') is False:
            os.makedirs(pathName + os.sep + expName + os.sep + 'pdata' + os.sep + '1')

        fidFile1 = pathName + os.sep + expName + os.sep + 'fid'
        spcFile1r = pathName + os.sep + expName + os.sep + 'pdata' + os.sep + '1' + os.sep + '1r'
        # write 1D FID file
        fid = np.zeros(2*len(self.fid[0]))
        fid[0::2] = self.fid[0].real
        fid[1::2] = self.fid[0].imag
        if self.acq.byteOrder == 1:
            fid.astype('>i4').tofile(fidFile1)
        else:
            if self.acq.dataType == 0:
                fid.astype('<i4').tofile(fidFile1)
            else:
                if self.acq.dataType == 1:
                    fid.astype(np.float64).tofile(fidFile1)
                else:
                    fid.astype(np.float32).tofile(fidFile1)


        # write 1r file
        spc = self.spc[0].real
        if scaleFactor == -1:
            scaleFactor = int(2*np.max(spc)/2147483647)

        spc /= scaleFactor
        spc.real.astype(np.int32).tofile(spcFile1r)
        # write parameter files
        # acqu file
        fName = pathName + os.sep + expName + os.sep + 'acqu'
        f = open(fName, 'w')
        f.write(self.acq.acqusText)
        f.close()
        # acqus file
        fName = pathName + os.sep + expName + os.sep + 'acqus'
        f = open(fName, 'w')
        f.write(self.acq.acqusText)
        f.close()
        # pulseprogram file
        fName = pathName + os.sep + expName + os.sep + 'pulseprogram'
        f = open(fName, 'w')
        f.write(self.pulseProgram)
        f.close()
        # audita.txt file
        fName = pathName + os.sep + expName + os.sep + 'audita.txt'
        f = open(fName, 'w')
        f.write(self.audita_txt)
        f.close()
        # format.temp file
        fName = pathName + os.sep + expName + os.sep + 'format.temp'
        f = open(fName, 'w')
        f.write(self.format_temp)
        f.close()
        # fq1list file
        fName = pathName + os.sep + expName + os.sep + 'fq1list'
        f = open(fName, 'w')
        f.write(self.fq1list)
        f.close()
        # scon2 file
        fName = pathName + os.sep + expName + os.sep + 'scon2'
        f = open(fName, 'w')
        f.write(self.scon2)
        f.close()
        # shimvalues file
        fName = pathName + os.sep + expName + os.sep + 'shimvalues'
        f = open(fName, 'w')
        f.write(self.shimvalues)
        f.close()
        # uxnmr.par file
        fName = pathName + os.sep + expName + os.sep + 'uxnmr.par'
        f = open(fName, 'w')
        f.write(self.uxnmr_par)
        f.close()
        # procs file
        origSI = self.proc.regEx.si.findall(self.proc.procsText)[0]
        ftMod  = int(self.proc.regEx.ftMod.findall(self.proc.procsText)[0])
        ncProc = int(self.proc.regEx.ncProc.findall(self.proc.procsText)[0])
        procsText = self.proc.procsText.replace(origSI, str(len(self.spc[0])))
        procsText = procsText.replace('FT_mod= ' + str(ftMod), 'FT_mod= 6')
        procsText = procsText.replace('NC_proc= '+ str(ncProc), 'NC_proc= 0')
        fName = pathName + os.sep + expName + os.sep + 'pdata' + os.sep + '1' + os.sep + 'procs'
        f = open(fName, 'w')
        f.write(procsText)
        f.close()
        # proc file
        fName = pathName + os.sep + expName + os.sep + 'pdata' + os.sep + '1' + os.sep + 'proc'
        f = open(fName, 'w')
        f.write(procsText)
        f.close()
        # outd file
        fName = pathName + os.sep + expName + os.sep + 'pdata' + os.sep + '1' + os.sep + 'outd'
        f = open(fName, 'w')
        f.write(self.outd)
        f.close()
        # title file
        fName = pathName + os.sep + expName + os.sep + 'pdata' + os.sep + '1' + os.sep + 'title'
        f = open(fName, 'w')
        f.write(self.title)
        f.close()
        # end exportBruker1d

    def fidOffsetCorrection(self, fid):
        fid = np.copy(fid)
        if (self.fidOffsetCorr > 0):
            mMean = np.mean(fid[:self.fidOffsetCorr])
            fid -= mMean

        return fid
        # end fidOffsetCorrection

    def gibbs(self, fid):
        if (self.proc.gibbs[0] == True):
            fid[0] /= 2.0

        return fid
        # end gibbs

    def hilbert(self, mat, dim=0):
        if (dim == 1):
            mat = np.ndarray.transpose(mat)

        npts = len(mat[0])
        npts1 = len(mat)
        v1 = np.ones(npts1)
        mat1 = np.array([[]], dtype='complex')
        mat1 = np.resize(mat1, (npts1, npts))
        bMat = np.zeros(int(2 * npts), dtype='complex')
        bMat[:(npts + 1)] = np.ones(npts + 1)
        bMat[1:npts] += bMat[1:npts]
        zMat = np.zeros(int(2 * npts), dtype='complex')
        bMat = np.outer(v1, bMat)
        zMat = np.outer(v1, zMat)
        zMat[:, :npts] = mat
        zMat = np.fft.ifft(bMat * np.fft.fft(zMat))
        mat = zMat[:, :npts]
        if (dim == 1):
            mat = np.ndarray.transpose(mat)

        return mat
        # end hilbert

    def hilbert1(self, mat, dim):
        if (self.dim == 1):
            npts = len(mat)
            bMat = np.zeros(int(2 * npts), dtype='complex')
            zMat = np.zeros(int(2 * npts), dtype='complex')
            zMat[:int(len(mat))] = mat
            bMat[:(npts + 1)] = np.ones(npts + 1)
            bMat[1:npts] += bMat[1:npts]
            mat2 = np.zeros(2 * npts, dtype='complex')
            mat2[:len(mat)] = mat
            mat2 = ifft(bMat * fft(mat2))
            mat1 = np.copy(mat2[:len(mat)])

        if (self.dim == 2):
            if (dim == 1):
                mat = np.ndarray.transpose(mat)

            npts = len(mat[0])
            npts1 = len(mat)
            mat1 = np.array([[]], dtype='complex')
            mat1 = np.resize(mat1, (npts1, npts))
            for k in range(len(mat)):
                bMat = np.zeros(int(2 * npts), dtype='complex')
                zMat = np.zeros(int(2 * npts), dtype='complex')
                zMat[:int(len(mat[k]))] = mat[k]
                bMat[:(npts + 1)] = np.ones(npts + 1)
                bMat[1:npts] += bMat[1:npts]
                mat2 = np.zeros(2 * npts, dtype='complex')
                mat2[:len(mat[k])] = mat[k]
                mat2 = ifft(bMat * fft(mat2))
                mat1[k] = mat2[:len(mat[k])]

            if (dim == 1):
                mat1 = np.ndarray.transpose(mat1)

        return mat1
        # end hilbert1

    def multiply(self, factor=1.0):
        self.spc = factor * self.spc

    def phase(self, ph0, ph1, npts):
        ph0 = -ph0 * math.pi / 180.0
        ph1 = -ph1 * math.pi / 180.0
        t = complex()
        frac = np.linspace(0, 1, npts)
        ph = ph0 + frac * ph1
        self.spc[0] = np.cos(ph) * self.spc[0].real + np.sin(ph) * self.spc[0].imag + 1j * (-np.sin(ph) * self.spc[0].real + np.cos(ph) * self.spc[0].imag)
        # end phase

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

    def phase2a(self, ph0, ph1, dim=0):
        mat = np.copy(self.spc.real)
        mat = self.hilbert(mat, dim)
        if (dim == 1):
            mat = np.ndarray.transpose(mat)

        npts = len(mat[0])
        npts1 = len(mat)
        v1 = np.ones(npts1)
        ph0 = -ph0 * math.pi / 180.0
        ph1 = -ph1 * math.pi / 180.0
        frac = np.outer(v1, np.linspace(0, 1, npts))
        ph = ph0 + frac * ph1
        mat = np.cos(ph) * mat.real + np.sin(ph) * mat.imag + 1j * (-np.sin(ph) * mat.real + np.cos(ph) * mat.imag)
        if (dim == 1):
            mat = np.ndarray.transpose(mat)

        self.spc = np.copy(mat.real)
        # end phase2a

    def phase2d(self, ph0, ph1, dim):
        mat = self.hilbert1(self.spc, dim)
        if (dim == 1):
            mat = np.ndarray.transpose(mat)

        npts = len(mat)
        for k in range(npts):
            mat[k] = self.phase2(mat[k], ph0, ph1)

        if (dim == 1):
            mat = np.ndarray.transpose(mat)

        self.spc = np.copy(mat.real)
        # end phase2d

    def phase3(self, mat, ph0, ph1):
        npts = len(mat)
        ph0 = -ph0 * math.pi / 180.0
        ph1 = -ph1 * math.pi / 180.0
        t = complex()
        for k in range(int(npts)):
            frac = float(k) / float(npts)
            ph = ph0 + frac * ph1
            t = complex(math.cos(ph) * mat[k].real + math.sin(ph) * mat[k].imag,
                        -math.sin(ph) * mat[k].real + math.cos(ph) * mat[k].imag)
            mat[k] = t

        return mat
        # end phase3

    def points2Hz(self, points, dim):
        sw = self.acq.sw_h[dim]
        if (dim == 0):
            npts = int(len(self.spc[0]))

        if (dim == 1):
            npts = int(len(self.spc))

        hz = sw * (points - npts / 2) / npts
        return hz
        # end points2Hz

    def points2ppm(self, points, dim):
        sw = self.acq.sw_h[dim]
        if (dim == 0):
            sfo = self.acq.sfo1
            npts = self.proc.nPoints[0] #int(len(self.spc[0]))

        if (dim == 1):
            if self.acq.manufacturer == 'Bruker':
                sfo = self.proc.sf[1]
                if sfo == 0.0:
                    sfo = self.acq.sfo2

            else:
                sfo = self.acq.sfo2

            npts = int(len(self.spc))

        ppm = (sw / sfo) * (points / (npts - 1) - self.refPoint[dim] / (npts - 1)) + self.refShift[dim]
        return ppm
        # end points2ppm

    def ppm2points(self, ppm, dim):
        sw = self.acq.sw_h[dim]
        if (dim == 0):
            sfo = self.acq.sfo1
            npts = int(len(self.spc[0]))

        if (dim == 1):
            if self.acq.manufacturer == 'Bruker':
                sfo = self.proc.sf[1]
            else:
                sfo = self.acq.sfo2

            npts = int(len(self.spc))

        points = np.round(((ppm - self.refShift[dim]) * (sfo / sw) + self.refPoint[dim] / (npts - 1)) * (npts - 1))
        return points.astype(int)
        # end ppm2points

    def ppm2points2d(self, ppm):
        points = np.copy(ppm)
        for k in range(len(ppm)):
            points[k][0] = self.ppm2points(ppm[k][0], 0)
            points[k][1] = self.ppm2points(ppm[k][1], 1)

        return points
        # end ppm2points2d

    def procSpc(self):
        if (self.dim == 1):
            self.procSpc1D()

        if (self.dim == 2):
            self.procSpc2D()

        #self.autoRef()
        self.calcPPM()
        # end procSpc

    def procSpc1D(self):
        fid = np.copy(self.fid[0])
        fid = self.waterSupp(fid)
        fid = self.fidOffsetCorrection(fid)
        fid = self.gibbs(fid)
        fid = self.apodise(fid, 0, self.proc.lb[0], self.proc.gb[0], self.proc.ssb[0], self.acq.groupDelay,
                           self.acq.sw_h[0])
        fid = self.zeroFill(fid)
        spc = fftshift(fft(fid))
        self.spc = np.resize(self.spc, (1, len(spc)))
        self.spc[0] = np.copy(spc)
        self.phase(0, 360.0 * self.acq.groupDelay, self.proc.nPoints[0])
        if self.proc.invertMatrix[0] == True:
            self.spc[0] = np.copy(np.flip(self.spc[0]))

        self.phase(self.proc.ph0[0], self.proc.ph1[0], self.proc.nPoints[0])
        # end procSpc1D

    def procSpc2D(self, testQuad2d=False, noAbs=False):
        fid = np.copy(self.fid)
        if self.proc.multFactor[0] == 0:
            self.proc.multFactor[0] = self.proc.nPoints[0] / len(self.fid[0])

        if self.proc.multFactor[1] == 0:
            self.proc.multFactor[1] = self.proc.nPoints[1] / len(self.fid)

        npts2 = len(self.spc)
        npts1 = len(self.spc[0])
        if npts1 > 0:
            if self.proc.stripStart < 1:
                stsr = 1
            else:
                stsr = self.proc.stripStart

            if self.proc.nPoints[0] != npts1:
                self.refPoint[0] = int((self.proc.refPoint[0])*self.proc.nPoints[0] / (self.proc.multFactor[0] * len(self.fid[0])))  #- stsr + 1

        if npts2 > 1:
            if self.proc.nPoints[1] != npts2:
                self.refPoint[1] = int(self.proc.refPoint[1] * self.proc.nPoints[1] / (self.proc.multFactor[1] * len(self.fid)))

        self.spc = np.copy(np.array([[]], dtype='complex'))
        self.spc = np.copy(np.resize(self.spc, (self.proc.nPoints[0], self.proc.nPoints[1])))
        self.spc *= 0
        if (self.proc.nPoints[0] > len(fid[0])):
            fid = np.resize(fid, (self.proc.nPoints[1], self.proc.nPoints[0]))
            fid *= 0
            for k in range(len(self.fid)):
                fid[k][:len(self.fid[k])] = np.copy(self.fid[k][:])

        if (self.proc.nPoints[0] < len(fid[0])):
            fid = np.resize(fid, (len(fid), self.proc.nPoints[0]))
            for k in range(len(fid)):
                fid[k] = np.copy(self.fid[k][:self.proc.nPoints[0]])

        for k in range(len(self.fid)):
            fid2 = np.copy(fid[k])
            fid2 = self.waterSupp(fid2)
            fid2 = self.fidOffsetCorrection(fid2)
            fid2 = self.gibbs(fid2)
            fid2 = self.apodise(fid2, 0, self.proc.lb[0], self.proc.gb[0], self.proc.ssb[0], self.acq.groupDelay,
                                self.acq.sw_h[0])
            fid2 = self.zeroFill(fid2, 0)
            fid2 = fftshift(fft(fid2))
            if self.acq.fnMode != 1:
                fid2 = self.phase2(fid2, 0, 360.0 * self.acq.groupDelay)

            if self.proc.invertMatrix[0] == True:
                fid2 = np.flip(fid2)

            if self.acq.fnMode != 1:
                fid2 = self.phase2(fid2, self.proc.ph0[0], self.proc.ph1[0])

            fid[k] = np.copy(fid2)

        fid = np.copy(np.ndarray.transpose(np.conj(fid)))
        if testQuad2d == False:
            fid = self.quad2D(fid)
            for k in range(len(fid)):
                fid2 = np.copy(fid[k])
                fid2 = self.gibbs(fid2)
                fid2 = self.apodise(fid2, 1, self.proc.lb[1], self.proc.gb[1], self.proc.ssb[1], 0.0, self.acq.sw_h[1])
                fid2 = self.zeroFill(fid2, 1)
                fid2 = fftshift(fft(fid2))
                if self.proc.invertMatrix[1] == True:
                    fid2 = np.flip(fid2)

                if (self.acq.fnMode != 1):
                    fid2 = self.phase2(fid2, self.proc.ph0[1], self.proc.ph1[1])
                    self.spc[k] = fid2.real
                else:
                    self.spc[k] = np.copy(fid2)

            self.spc = np.copy(np.ndarray.transpose(self.spc))
            if ((self.acq.fnMode == 1) and (self.proc.tilt == True)):
                #self.spc = np.copy(self.hilbert(self.spc, 0))
                self.tiltJRes()

            if self.acq.fnMode == 1 and noAbs is False:
                for k in range(len(self.spc)):
                    self.spc[k] = np.copy(np.abs(self.spc[k]))

                if (self.proc.symj == True):
                    self.symjres()

        else:
            self.spc = np.copy(self.fid)

        if self.proc.stripStart == 0:
            stsr = 1
        else:
            stsr = self.proc.stripStart

        if self.proc.stripEnd > stsr:
            stsi = self.proc.stripEnd
            self.spc = np.ndarray.transpose(self.spc)
            self.spc = self.spc[stsr-1:stsi]
            self.ppm1 = self.ppm1[stsr-1:stsi]
            self.spc = np.ndarray.transpose(self.spc)

        self.spc = np.copy(self.spc.real)
        # end procSpc2D

    def quad2D(self, fid):
        fid = np.copy(fid)
        inc = self.acq.fnMode
        # MetLab: 6,     0,    3,      2,            1,   4
        # FnMode: 1,     2,    3,      4,            5,   6
        #         jres, QF, TPPI, States, States-TPPI , E/A
        rFid = np.copy(fid)
        rFid = np.resize(rFid, (len(fid), int(len(fid[0]) / 2)))
        iFid = np.copy(fid)
        iFid = np.resize(iFid, (len(fid), int(len(fid[0]) / 2)))
        if inc == 1 or inc == 2:
            fid = fid  # rFid + 1j * iFid

        if inc == 5:
            for k in range(len(fid)):
                npts = int(len(fid[k]) / 2)
                rFid[k] = fid[k][0::2].real * (-2 * (np.linspace(0, npts - 1, npts, dtype='int') % 2) + 1)
                iFid[k] = fid[k][1::2].real * (2 * (np.linspace(0, npts - 1, npts, dtype='int') % 2) - 1)

            fid = np.copy(rFid + 1j * iFid)

        if inc == 6:
            for k in range(len(fid)):
                rFid[k] = np.conj(fid[k][0::2] + fid[k][1::2])
                iFid[k] = np.conj(fid[k][0::2] - fid[k][1::2])
                iFid[k] = iFid[k].imag - 1j * iFid[k].real
                rFid[k] = rFid[k].imag + 1j * iFid[k].imag

            fid = np.copy(rFid)

        return fid
        # end quad2D

    def readPipe2D(self, pName, fName):
        npd = nmrpipeData.NmrPipeData()
        npd.readPipe(pName, fName)
        self.spc = npd.spc
        self.proc.sw_h[0] = npd.fdf2sw
        self.refShift[0] = npd.fdf2orig / self.acq.sfo1
        self.refPoint[0] = 0
        self.proc.sw_h[1] = npd.fdf1sw
        self.refShift[1] = npd.fdf1orig / self.acq.sfo2
        self.refPoint[1] = 0
        self.proc.phaseInversion = False
        self.proc.multFactor = [1, 1]
        self.proc.nPoints[0] = len(self.spc[0])
        self.proc.nPoints[1] = len(self.spc)
        self.calcPPM()
        # end readPipe2D

    def readSpc(self):
        self.acq.read(self.dataSetName + os.sep + self.dataSetNumber)
        self.proc.read(self.dataSetName + os.sep + self.dataSetNumber)
        # self.proc.sw_h = self.acq.sw_h
        self.origDataSet = self.dataSetName + os.sep + self.dataSetNumber
        if self.acq.manufacturer == 'Bruker':
            titleFile = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'pdata' + os.sep + '1' + os.sep + 'title'
            if (os.path.isfile(titleFile)):
                fid = open(titleFile, "r")
                self.title = fid.read()
                fid.close()

            acquFile = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'acqu'
            if (os.path.isfile(acquFile)):
                fid = open(acquFile, "r")
                self.acqu = fid.read()
                fid.close()

            acqu2File = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'acqu2'
            if (os.path.isfile(acqu2File)):
                fid = open(acqu2File, "r")
                self.acqu2 = fid.read()
                fid.close()

            acqu3File = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'acqu3'
            if (os.path.isfile(acqu3File)):
                fid = open(acqu3File, "r")
                self.acqu3 = fid.read()
                fid.close()

            auditaFile = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'audita.txt'
            if (os.path.isfile(auditaFile)):
                fid = open(auditaFile, "r")
                self.audita_txt = fid.read()
                fid.close()

            formatFile = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'format.temp'
            if (os.path.isfile(formatFile)):
                fid = open(formatFile, "r")
                self.format_temp = fid.read()
                fid.close()

            fq1listFile = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'fq1list'
            if (os.path.isfile(fq1listFile)):
                fid = open(fq1listFile, "r")
                self.fq1list = fid.read()
                fid.close()

            scon2File = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'scon2'
            if (os.path.isfile(scon2File)):
                fid = open(scon2File, "r")
                self.scon2 = fid.read()
                fid.close()

            shimvaluesFile = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'shimvalues'
            if (os.path.isfile(shimvaluesFile)):
                fid = open(shimvaluesFile, "r")
                self.shimvalues = fid.read()
                fid.close()

            uxnmrparFile = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'uxnmr.par'
            if (os.path.isfile(uxnmrparFile)):
                fid = open(uxnmrparFile, "r")
                self.uxnmr_par = fid.read()
                fid.close()

            procFile = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'pdata' + os.sep + '1' + os.sep + 'proc'
            if (os.path.isfile(procFile)):
                fid = open(procFile, "r")
                self.proc1 = fid.read()
                fid.close()

            proc2File = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'pdata' + os.sep + '1' + os.sep + 'proc2'
            if (os.path.isfile(proc2File)):
                fid = open(proc2File, "r")
                self.proc2 = fid.read()
                fid.close()

            proc3File = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'pdata' + os.sep + '1' + os.sep + 'proc3'
            if (os.path.isfile(proc3File)):
                fid = open(proc3File, "r")
                self.proc3 = fid.read()
                fid.close()

            outdFile = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'pdata' + os.sep + '1' + os.sep + 'outd'
            if (os.path.isfile(outdFile)):
                fid = open(outdFile, "r")
                self.outd = fid.read()
                fid.close()

            self.acq.sfo2 = self.proc.sf[1]
            self.display.yLabel = self.proc.axisNucleus[1]

        else:
            titleFile1 = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'text'
            titleFile2 = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'TEXT'
            if (os.path.isfile(titleFile1)):
                fid = open(titleFile1, "r")
                self.title = fid.read()
                fid.close()

            if (os.path.isfile(titleFile2)):
                fid = open(titleFile2, "r")
                self.title = fid.read()
                fid.close()

        if self.acq.manufacturer == 'Bruker':
            pulProgFile = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'pulseprogram'
        else:
            pulProgFile = self.dataSetName + os.sep + self.dataSetNumber + os.sep + self.acq.pulProgName + ".c"

        if (os.path.isfile(pulProgFile)):
            fid = open(pulProgFile, "r")
            self.pulseProgram = fid.read()
            fid.close()

        if self.acq.manufacturer == 'Bruker':
            fidFile1 = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'fid'
            spcFile1r = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'pdata' + os.sep + '1' + os.sep + '1r'
            spcFile1i = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'pdata' + os.sep + '1' + os.sep + '1i'
            fidFile2 = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'ser'
            if (os.path.isfile(fidFile1)):
                # read 1D FID file
                self.fid = np.resize(self.fid, (1, int(self.acq.nDataPoints[0] / 2)))
                f = open(fidFile1, 'rb')
                if self.acq.byteOrder == 1:
                    fid = np.fromfile(f, dtype='>i4')
                else:
                    if self.acq.dataType == 0:
                        fid = np.fromfile(f, dtype='<i4')
                    else:
                        if self.acq.dataType == 1:
                            fid = np.fromfile(f, dtype=np.float64)
                        else:
                            if self.acq.dataType == 2:
                                fid = np.fromfile(f, dtype=np.double)
                            else:
                                fid = np.fromfile(f, dtype=np.float32)


                f.close()
                self.fid[0].real = fid[0::2]
                self.fid[0].imag = -fid[1::2]
                self.dim = 1

            if (os.path.isfile(spcFile1r)):
                # read 1D spectrum file (real part)
                f = open(spcFile1r, 'rb')
                fid = np.fromfile(f, dtype=np.int32)
                self.spc = np.resize(self.spc, (1, int(len(fid))))
                self.spc[0].real = fid;
                f.close()
                if (os.path.isfile(spcFile1i)):
                    # read 1D spectrum file (imaginary part)
                    f = open(spcFile1i, 'rb')
                    fid = np.fromfile(f, dtype=np.int32)
                    f.close()
                    self.spc[0].imag = fid


            elif (os.path.isfile(fidFile2)):
                # read 2D spectrum
                np1 = int(self.acq.nDataPoints[0])
                np2 = int(self.acq.nDataPoints[1])
                self.fid = np.resize(self.fid, (int(np2), int(np1 / 2)))
                f = open(fidFile2, 'rb')
                fid = np.fromfile(f, dtype=np.int32)
                f.close()
                fid = fid.reshape(int(np2), int(np1))
                for x in range(np2):
                    self.fid[x].real = fid[x][0::2]
                    self.fid[x].imag = -fid[x][1::2]

                self.dim = 2
                if self.acq.pulProgName.find("jres"):
                    self.jRes = True
                    self.proc.tilt = True
                    self.proc.symj = True

            self.proc.sw_h = np.copy(self.acq.sw_h)
        elif self.acq.manufacturer == 'Varian':
            fidFile1 = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'fid'
            fidFile2 = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'FID'
            fidFile = ''
            if os.path.isfile(fidFile1):
                fidFile = fidFile1
            elif os.path.isfile(fidFile2):
                fidFile = fidFile2

            if len(fidFile) > 0:
                self.dim = self.acq.np + self.acq.np2 - 1
                if self.dim == 1 and self.acq.ni > 1:
                    self.dim = 2
                    self.jRes = True
                    self.proc.tilt = True
                    self.proc.symj = True
                    self.acq.fnMode = 1

                f = open(fidFile, 'rb')
                dh1 = np.fromfile(f, dtype='>i4', count=6)
                dh2 = np.fromfile(f, dtype='>i2', count=2)
                dh3 = np.fromfile(f, dtype='>i4', count=1)
                status = np.binary_repr(dh2[1], 16)
                dataFormat = '>i2'
                if status[12] == '1':
                    dataFormat = '>f'
                elif status[13] == '1':
                    dataFormat = '>i4'

                nBlocks = dh1[0]
                nTraces = dh1[1]
                nPoints = dh1[2]
                if self.dim == 1:
                    self.fid = np.resize(self.fid, (1, int(self.acq.nDataPoints[0] / 2)))
                    fid = np.array([])
                    for k in range(nBlocks):
                        bh1 = np.fromfile(f, dtype='>i2', count=4)
                        bh2 = np.fromfile(f, dtype='>i4', count=1)
                        bh3 = np.fromfile(f, dtype='>f', count=4)
                        if status[10] == '1':
                            bh4 = np.fromfile(f, dtype='>i2', count=4)
                            bh5 = np.fromfile(f, dtype='>i4', count=1)
                            bh6 = np.fromfile(f, dtype='>f', count=4)

                        data = np.fromfile(f, dtype=dataFormat, count=nTraces * nPoints)
                        fid = np.append(fid, data[0::2] + 1j * data[1::2])

                    self.fid[0] = fid

                else:
                    ni = int(self.acq.nDataPoints[1] / 2)
                    td = int(self.acq.nDataPoints[0] / 2)
                    self.fid = np.zeros((ni, td), dtype='complex')
                    for k in range(ni):
                        bh1 = np.fromfile(f, dtype='>i2', count=4)
                        bh2 = np.fromfile(f, dtype='>i4', count=1)
                        bh3 = np.fromfile(f, dtype='>f', count=4)
                        status = np.binary_repr(bh1[1], 8)
                        dataFormat = '>i2'
                        if status[4] == '1':
                            dataFormat = '>f'
                        elif status[5] == '1':
                            dataFormat = '>i4'

                        data = np.fromfile(f, dtype=dataFormat, count=nTraces * nPoints)
                        self.fid[k] = data[0::2] + 1j * data[1::2]

                f.close()

        # end readSpc

    def setRef(self, refShift, refPoint):
        for k in range(len(refShift)):
            self.refShift[k] = refShift[k]

        for k in range(len(refPoint)):
            self.refPoint[k] = refPoint[k]

        self.calcPPM()
        # end setRef    

    def setWindowFunction(self, dim, wdwf):
        try:
            self.proc.windowType[dim] = self.windowFunction[wdwf]

        except:
            self.proc.windowType[dim] = self.windowFunction['none']
            print('Unknown windowType, setting windowFunction to "none"')

        # end setWindowFunction

    def smo(self, fid):
        x = np.linspace(0, len(fid) - 1, len(fid))
        f0 = np.copy(fid)
        fid = np.roll(fid, math.floor(-self.acq.groupDelay))
        pp = np.polyfit(x, fid, self.proc.polyOrder)
        fid = fid - np.polyval(pp, x)
        fid = np.roll(fid, math.ceil(self.acq.groupDelay))
        return fid
        # end smo

    def symjres(self):
        for k in range(int(len(self.spc) / 2)):
            tmp = np.minimum(self.spc[k], self.spc[len(self.spc) - k - 1])
            self.spc[k] = tmp
            self.spc[len(self.spc) - k - 1] = tmp

        # end symj

    def tiltJRes(self):
        hzPerPoint = self.acq.sw_h[0] / len(self.spc[0])
        sw2 = self.acq.sw_h[1]
        npts2 = len(self.spc)
        hzVect = np.linspace( sw2 / 2.0, -sw2 / 2.0, npts2)
        #hzVect = hzVect - (hzVect[1] - hzVect[0]) / 2.0
        for k in range(npts2):
            npts = hzVect[k] / hzPerPoint
            fid1 = np.copy(ifft(self.spc[k]))
            fid1 = self.phase2(fid1, 0, -npts * 360.0)

            self.spc[k] = np.copy(fft(fid1))

        # end tiltJRes

    def waterSupp(self, fid):
        if (self.proc.waterSuppression == 1):
            fid = self.conv(fid)

        if (self.proc.waterSuppression == 2):
            fid = self.smo(fid)

        return fid
        # end waterSupp

    def zeroFill(self, fid, dim=0):
        fid1 = np.zeros(self.proc.nPoints[dim], dtype='complex')
        fid1[:int(len(fid))] = fid
        return fid1
        # end zeroFill
