import numpy as np
from metabolabpy.nmr import wdwf
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
        self.fid            = np.array([[]], dtype = 'complex')
        self.spc            = np.array([[]], dtype = 'complex')
        self.ppm1           = np.array([], dtype = 'float64')
        self.ppm2           = np.array([], dtype = 'float64')
        self.ppm3           = np.array([], dtype = 'float64')
        self.dim            = 0
        self.title          = np.array([], dtype = 'str')
        self.origDataSet    = str('')
        self.phCorrMode     = 0
        self.acq            = acqPars.AcqPars()
        self.proc           = procPars.ProcPars()
        self.disp           = dispPars.DispPars()
        self.fidOffsetCorr  = 0
        self.dataSetName    = ''
        self.dataSetNumber  = ''
        self.title          = str('Empty NMR data set')
        self.pulseProgram   = str('')
        self.windowFunction = {'none':0, 'exponential':1, 'gaussian':2, 'sine':3, 'qsine':4, 'sem':5}
        self.refShift       = np.array([0, 0, 0], dtype = 'float64')
        self.refPoint       = np.array([0, 0, 0], dtype = 'int')
        self.refsw          = np.array([0, 0, 0], dtype = 'float64')
        self.refTmspRange   = 0.3 # [ppm]
        self.apc            = apcbc.Apcbc()
        # end __init__

    def __str__(self):
        rString  = 'pyMetaboLab Data Set (v. 0.1)\n'
        rString += '__________________________________________________________________\n'
        rString += self.title
        return rString
        # end __str__
    
    def calcPPM(self):
        if(self.dim==1):
            self.ppm.resize(1,self.proc.nPoints[0])
            x           = np.linspace(0.0, self.proc.nPoints[0] - 1.0, self.proc.nPoints[0])
            self.ppm[0] = self.proc.refShift[0] + (self.proc.refPoint[0] - x)*self.acq.sw[0]/self.proc.nPoints[0]
        
        # end calcPPM
        
    def apodise(self, fid, dim, lb, gb, ssb, groupDelay, sw_h):
        if(self.proc.windowType[dim]==0):            # no window
            wdwf = np.ones(len(fid))
            
        if(self.proc.windowType[dim]==1):            # exponential window
            t    = (np.linspace(0.0, len(fid) - 1, len(fid)) - groupDelay)/sw_h
            wdwf = np.exp(-lb*t)
            
        if(self.proc.windowType[dim]==2):            # gaussian window
            t    = (np.linspace(0.0, len(fid) - 1 - groupDelay, len(fid)))/sw_h
            wdwf = np.exp(-lb*2*math.pi*t - 2*math.pi*(t**2)/((2*math.pi*gb*len(fid)/sw_h)))
            
        if(self.proc.windowType[dim]==3):            # sine window
            t    = (np.linspace(0.0, len(fid) - 1, len(fid)))/len(fid)
            ssb2 = ssb*math.pi/180.0
            wdwf = np.sin(ssb + (math.pi - ssb)*t)
            
        if(self.proc.windowType[dim]==4):            # qsine window
            t    = (np.linspace(0.0, len(fid) - 1, len(fid)))/len(fid)
            ssb2 = ssb*math.pi/180.0
            wdwf = np.sin(ssb + (math.pi - ssb)*t)**2
            
        if(self.proc.windowType[dim]==5):            # sem window
            t1   = (np.linspace(0.0, len(fid) - 1 - groupDelay, len(fid)))/sw_h
            t2   = (np.linspace(0.0, len(fid) - 1, len(fid)))/len(fid)
            ssb2 = ssb*math.pi/180.0
            wdwf = np.exp(-lb*t1)*np.sin(ssb + (math.pi - ssb)*t2)
            
        fid *= wdwf
        return fid
        # end apodise
    
    def autobaseline1d(self):
        spc                                 = self.spc[0]
        self.apc.npts                       = len(spc)
        scaleFact                           = np.max(np.abs(spc))
        spc                                /= scaleFact
        xaxis                               = np.linspace(-self.apc.nMax, self.apc.nMax, self.apc.npts)
        parEval                             = self.apc.fitBaseline(spc, xaxis)
        spc2                                = self.apc.baselineFitFuncEval(parEval, spc, xaxis, False)
        self.apc.pars                       = parEval
        self.apc.rSpc                       = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self.apc.iSpc                       = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self.apc.rSpc[:int(len(parEval)/2)] = parEval[:int(len(parEval)/2)] 
        self.apc.iSpc[:int(len(parEval)/2)] = parEval[int(len(parEval)/2):] 
        self.spc[0]                         = spc2*scaleFact
        self.apc.correctBaseline            = 1
        self.procSpc1D()
        self.baseline1d()
        # end autobaseline1d

    def autophase1d(self):
        spc               = self.spc[0]
        self.apc.npts     = len(spc)
        scaleFact         = np.max(np.abs(spc))
        spc              /= scaleFact
        xaxis             = np.linspace(-self.apc.nMax, self.apc.nMax, self.apc.npts)
        parEval           = self.apc.fitPhase(spc, xaxis)
        spc2              = self.apc.phaseFitFuncEval(parEval, spc, xaxis, False)
        if(np.min(spc2.real) == -np.max(np.abs(spc2.real))):
            parEval[0]    = ((parEval[0]*self.apc.mFact0 + 180.0) % 360)/self.apc.mFact0
    
        parEval[0]                              = (((parEval[0]*self.apc.mFact0 + 180.0) % 360) - 180.0)/self.apc.mFact0
        spc2                                    = self.apc.phaseFitFuncEval(parEval, spc, xaxis, False)
        self.apc.pars                           = parEval
        self.apc.rSpc                           = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self.apc.iSpc                           = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self.apc.rSpc[:int((len(parEval)-2)/2)] = parEval[2:2+int((len(parEval)-2)/2)] 
        self.apc.iSpc[:int((len(parEval)-2)/2)] = parEval[2+int((len(parEval)-2)/2):] 
        self.spc[0]                             = spc2*scaleFact
        self.proc.ph0[0]                       += parEval[0]*self.apc.mFact0
        self.proc.ph1[0]                       += parEval[1]*self.apc.mFact1
        self.apc.correctBaseline                = 1
        self.autoRef()
        self.procSpc1D()
        self.baseline1d()
        # end autophase1d

    def autoRef(self,tmsp=True):
        self.refShift[0] = self.acq.o1/self.acq.bf1
        self.refPoint[0] = int(len(self.spc[0])/2)
        if(self.dim==2):
            self.refShift[1] = self.acq.o2/self.acq.bf2
            self.refPoint[1] = int(len(self.spc)/2)
            
        if(self.dim==1):
            if(tmsp==True):
                self.refPoint[0]  = self.ppm2points(0.0, 0)
                self.refShift[0]  = 0.0
                pts               = self.ppm2points(np.array([-self.refTmspRange, self.refTmspRange]),0)
                npts              = len(self.spc[0])
                r                 = np.arange(npts-max(pts),npts-min(pts))
                spc               = self.spc[0][r].real
                refP              = np.where(spc == np.amax(spc))
                self.refPoint[0] -= refP[0][0] - int((max(pts)-min(pts))/2) + 1
                self.calcPPM()
                
                
        # end autoRef

    def baseline1d(self):
        spc               = self.spc[0]
        self.apc.npts     = len(spc)
        #scaleFact         = np.max(np.abs(spc))
        #spc              /= scaleFact
        xaxis             = np.linspace(-self.apc.nMax, self.apc.nMax, self.apc.npts)
        parEval           = np.copy(self.apc.rSpc)
        parEval           = np.append(parEval, self.apc.iSpc)
        print(parEval)
        spc2              = self.apc.baselineFitFuncEval(parEval, spc, xaxis, False)
        self.spc[0]       = spc2 #*scaleFact
        # end baseline1d
        
    def calcPPM(self):
        npts      = int(len(self.spc[0]))
        self.ppm1 = self.points2ppm(np.linspace(npts-1, 0, npts), 0)
        if(self.dim>1):
            npts      = int(len(self.spc))
            self.ppm2 = self.points2ppm(np.linspace(npts-1, 0, npts), 1)
        
        # end calcPPM
    
    def conv(self, fid):
        f1 = fid
        x       = np.linspace(0, len(fid)-1,len(fid))
        f0      = np.copy(fid)
        filtFid = fid[list(range(round(self.acq.groupDelay)))]
        fid     = fid[list(range(round(self.acq.groupDelay),len(fid)))]
        ws2  = self.proc.convWindowSize[0]/2
        es   = self.proc.convExtrapolationSize[0]
        wt   = int(self.proc.convWindowType[0])
        win  = np.zeros(int(ws2*2 + 1))
        sWin = 0
        for k in range(int(ws2*2 + 1)):
            if(wt==0):  # Gaussian window
                win[k] = math.exp(-4.0*((k-ws2)**2)/(ws2**2))
                sWin  += win[k]
            else:       # sine bell window
                win[k] = math.cos((math.pi*(k-ws2))/(2*ws2+2))
                sWin  += win[k]
                
        fid2            = np.convolve(win, fid)/sWin
        # Extrapolation of first sw2 data points
        idx             = np.linspace(2*ws2-1,0,2*ws2, dtype = 'int')
        idx2            = np.linspace(0, 2*ws2-1, 2*ws2, dtype = 'int')
        fid2[idx2].real = fid2[int(2*ws2)].real + np.mean(np.diff(fid2[idx + int(2*ws2) - 1].real))*idx
        fid2[idx2].imag = fid2[int(2*ws2)].imag + np.mean(np.diff(fid2[idx + int(2*ws2) - 1].imag))*idx
        # Extrapolation of last sw2 data points
        idx             = np.linspace(0, 2*ws2 - 1, 2*ws2, dtype = 'int')
        idx2            = np.linspace(len(fid) - int(2*ws2) - 1, len(fid) - 1, int(2*ws2), dtype = 'int')
        fid2[idx2].real = fid2[len(fid) - int(2*ws2) - 1].real - np.mean(np.diff(fid2[len(fid) - int(2*ws2) - 1 - idx].real))*idx
        fid2[idx2].imag = fid2[len(fid) - int(2*ws2) - 1].imag - np.mean(np.diff(fid2[len(fid) - int(2*ws2) - 1 - idx].imag))*idx
        fid2            = np.delete(fid2,  np.linspace(0, 2*ws2-1, 2*ws2, dtype = 'int'))
        fid            -= fid2 #(fidRe + 1j*fidIm)
        fid             = np.concatenate([filtFid, fid])
        return fid
        # end conv

    def fidOffsetCorrection(self, fid):
        if(self.fidOffsetCorr>0):
            meanRange = np.linspace(len(fid) - self.fidOffsetCorr - 1, len(fid) - 1, self.fidOffsetCorr, dtype = 'int')
            rMean     = np.mean(fid[meanRange].real)
            iMean     = np.mean(fid[meanRange].imag)
            fid      -= (rMean + 1j*iMean)
            
        return fid
        # end fidOffsetCorrection
    
    def gibbs(self, fid):
        if(self.proc.gibbs[0]==True):
            fid[0] /= 2.0
            
        return fid
        # end gibbs
    
    def hilbert(self, mat, dim=0):
        if(dim==1):
            mat = np.ndarray.transpose(mat)
                
        npts                        = len(mat[0])
        npts1                       = len(mat)
        v1                          = np.ones(npts1)
        mat1                        = np.array([[]], dtype = 'complex')
        mat1                        = np.resize(mat1,(npts1,npts))
        bMat                    = np.zeros(int(2*npts), dtype = 'complex')
        bMat[:(npts+1)]         = np.ones(npts+1)
        bMat[1:npts]           += bMat[1:npts]
        zMat                    = np.zeros(int(2*npts), dtype = 'complex')
        bMat = np.outer(v1,bMat)
        zMat = np.outer(v1,zMat)
        zMat[:,:npts] = mat
        zMat = np.fft.ifft(bMat*np.fft.fft(zMat))
        mat = zMat[:,:npts]            
        if(dim==1):
            mat = np.ndarray.transpose(mat)
        
        return mat
        # end hilbert

    def hilbert1(self, mat, dim):
        if(self.dim==1):
            npts                 = len(mat)
            bMat                 = np.zeros(int(2*npts), dtype = 'complex')
            zMat                 = np.zeros(int(2*npts), dtype = 'complex')
            zMat[:int(len(mat))] = mat
            bMat[:(npts+1)]      = np.ones(npts+1)
            bMat[1:npts]        += bMat[1:npts]
            mat2                 = np.zeros(2*npts, dtype = 'complex')
            mat2[:len(mat)]      = mat
            mat2                 = ifft(bMat * fft(mat2))
            mat1                 = np.copy(mat2[:len(mat)])
            
        if(self.dim==2):
            if(dim==1):
                mat = np.ndarray.transpose(mat)
                
            npts                        = len(mat[0])
            npts1                       = len(mat)
            mat1                        = np.array([[]], dtype = 'complex')
            mat1                        = np.resize(mat1,(npts1,npts))
            for k in range(len(mat)):
                bMat                    = np.zeros(int(2*npts), dtype = 'complex')
                zMat                    = np.zeros(int(2*npts), dtype = 'complex')
                zMat[:int(len(mat[k]))] = mat[k]
                bMat[:(npts+1)]         = np.ones(npts+1)
                bMat[1:npts]           += bMat[1:npts]
                mat2                    = np.zeros(2*npts, dtype = 'complex')
                mat2[:len(mat[k])]      = mat[k]
                mat2                    = ifft(bMat * fft(mat2))
                mat1[k]                 = mat2[:len(mat[k])]

            if(dim==1):
                mat1 = np.ndarray.transpose(mat1)
                
            
        return mat1
        # end hilbert1

    def phase(self, ph0, ph1, npts):
        ph0         = -ph0*math.pi/180.0
        ph1         = -ph1*math.pi/180.0
        t           = complex()
        frac        = np.linspace(0, 1, npts)
        ph          = ph0 + frac*ph1
        self.spc[0] = np.cos(ph)*self.spc[0].real + np.sin(ph)*self.spc[0].imag + 1j*(-np.sin(ph)*self.spc[0].real + np.cos(ph)*self.spc[0].imag)
        # end phase

    def phase2(self, mat, ph0, ph1):
        npts  = len(mat)
        ph0   = -ph0*math.pi/180.0
        ph1   = -ph1*math.pi/180.0
        t     = complex()
        frac  = np.linspace(0, 1, npts)
        ph             = ph0 + frac*ph1
        mat   = np.cos(ph)*mat.real + np.sin(ph)*mat.imag + 1j*(-np.sin(ph)*mat.real + np.cos(ph)*mat.imag)
        return mat
        # end phase2

    def phase2a(self, ph0, ph1, dim=0):
        mat = np.copy(self.spc)
        if(mat.imag.sum()==0):
            mat = self.hilbert(mat,dim)
            
        if(dim==1):
            mat = np.ndarray.transpose(mat)
            
        npts  = len(mat[0])
        npts1 = len(mat)
        v1    = np.ones(npts1)
        ph0   = -ph0*math.pi/180.0
        ph1   = -ph1*math.pi/180.0
        frac  = np.outer(v1,np.linspace(0, 1, npts))
        ph    = ph0 + frac*ph1
        mat   = np.cos(ph)*mat.real + np.sin(ph)*mat.imag + 1j*(-np.sin(ph)*mat.real + np.cos(ph)*mat.imag)
        if(dim==1):
            mat = np.ndarray.transpose(mat)
            
        self.spc = np.copy(mat)
        # end phase2a

    def phase2d(self, ph0, ph1, dim):
        mat = self.hilbert1(self.spc, dim)
        if(dim==1):
            mat = np.ndarray.transpose(mat)
            
        npts = len(mat)
        for k in range(npts):
            mat[k] = self.phase2(mat[k], ph0, ph1)
        
        if(dim==1):
            mat = np.ndarray.transpose(mat)
        
        self.spc = np.copy(mat.real)
        # end phase2d
        
    def phase3(self, mat, ph0, ph1):
        npts = len(mat)
        ph0  = -ph0*math.pi/180.0
        ph1  = -ph1*math.pi/180.0
        t    = complex()
        for k in range(int(npts)):
            frac           = float(k)/float(npts)
            ph             = ph0 + frac*ph1
            t = complex(math.cos(ph)*mat[k].real + math.sin(ph)*mat[k].imag, -math.sin(ph)*mat[k].real + math.cos(ph)*mat[k].imag)
            mat[k] = t
        
        return mat
        # end phase3

    def points2ppm(self, points, dim):
        sw = self.proc.sw_h[dim]
        if(dim==0):
            sfo  = self.acq.sfo1
            npts = int(len(self.spc[0]))
            
        if(dim==1):
            sfo  = self.acq.sfo2
            npts = int(len(self.spc))
        
        ppm = (sw/sfo)*(points/(npts-1) - self.refPoint[dim]/(npts-1)) + self.refShift[dim]
        return ppm
        # end points2ppm
    
    def ppm2points(self, ppm, dim):
        sw = self.proc.sw_h[dim]
        if(dim==0):
            sfo  = self.acq.sfo1
            npts = int(len(self.spc[0]))
            
        if(dim==1):
            sfo  = self.acq.sfo2
            npts = int(len(self.spc))
        
        points = np.round(((ppm - self.refShift[dim])*(sfo/sw) + self.refPoint[dim]/(npts-1))*(npts-1))
        return points.astype(int)
        # end ppm2points

    def ppm2points2d(self, ppm):
        points    = np.copy(ppm)
        for k in range(len(ppm)):
            points[k][0] = self.ppm2points(ppm[k][0], 0)
            points[k][1] = self.ppm2points(ppm[k][1], 1)
            
        return points
        # end ppm2points2d
    
    def procSpc(self):
        if(self.dim==1):
            self.procSpc1D()
            
        if(self.dim==2):
            self.procSpc2D()
            
        self.autoRef()
        self.calcPPM()
        # end procSpc

    def procSpc1D(self):
        fid         = np.copy(self.fid[0])
        fid         = self.waterSupp(fid)
        fid         = self.fidOffsetCorrection(fid)
        fid         = self.gibbs(fid)
        fid         = self.apodise(fid, 0, self.proc.lb[0], self.proc.gb[0], self.proc.ssb[0], self.acq.groupDelay, self.acq.sw_h[0])
        fid         = self.zeroFill(fid)
        spc         = fftshift(fft(fid))
        self.spc    = np.resize(self.spc,(1,len(spc)))
        self.spc[0] = np.copy(spc)
        self.phase(self.proc.ph0[0], 360.0*self.acq.groupDelay + self.proc.ph1[0], self.proc.nPoints[0])
        # end procSpc1D
    
    def procSpc2D(self):
        fid         = np.copy(self.fid)
        self.spc    = np.resize(self.spc,(self.proc.nPoints[0],self.proc.nPoints[1]))
        if(self.proc.nPoints[0]<len(fid[0])):
            fid = np.resize(fid,(len(fid), self.proc.nPoints[0]))
            for k in range(len(fid)):
                fid[k] = self.fid[k][:self.proc.nPoints[0]]
        
        for k in range(len(self.fid)):
            fid2   = fid[k]
            fid2   = self.waterSupp(fid2)
            fid2   = self.fidOffsetCorrection(fid2)
            fid2   = self.gibbs(fid2)
            fid2   = self.apodise(fid2, 0, self.proc.lb[0], self.proc.gb[0], self.proc.ssb[0], self.acq.groupDelay, self.acq.sw_h[0])
            fid2   = self.zeroFill(fid2, 0)
            fid2   = fftshift(fft(fid2))
            fid2   = self.phase2(fid2,self.proc.ph0[0], 360.0*self.acq.groupDelay + self.proc.ph1[0])
            fid[k] = fid2
            
        fid  = np.ndarray.transpose(fid)
        fid  = self.quad2D(fid)
        for k in range(len(fid)):
            fid2   = np.copy(fid[k])
            fid2   = self.gibbs(fid2)
            fid2   = self.apodise(fid2, 1, self.proc.lb[1], self.proc.gb[1], self.proc.ssb[1], 0.0, self.acq.sw_h[1])
            fid2   = self.zeroFill(fid2, 1)
            fid2   = fftshift(fft(fid2))
            fid2   = self.phase2(fid2,self.proc.ph0[1], self.proc.ph1[1])
            self.spc[k] = fid2
            
        self.spc = np.ndarray.transpose(self.spc)
        # end procSpc2D

    def quad2D(self, fid):
        inc  = self.acq.fnMode
        # MetLab: 6,     0,    3,      2,            1,   4
        # FnMode: 1,     2,    3,      4,            5,   6
        #         jres, QF, TPPI, States, States-TPPI , E/A
        rFid = np.copy(fid)
        rFid = np.resize(rFid,(len(fid),int(len(fid[0])/2)))
        iFid = np.copy(fid)
        iFid = np.resize(iFid,(len(fid),int(len(fid[0])/2)))
        if(inc == 1):
            fid = rFid + 1j*iFid
        
        if(inc == 2):
            fid = rFid + 1j*iFid
        
        if(inc == 3):
            fid = rFid + 1j*iFid
        
        if(inc == 4):
            fid = rFid + 1j*iFid
        
        if(inc == 5):
            fid = rFid + 1j*iFid
        
        if(inc == 6):
            for k in range(len(fid)):
                rFid[k] = np.conj(fid[k][0::2] + fid[k][1::2])
                iFid[k] = np.conj(fid[k][0::2] - fid[k][1::2])
                iFid[k] = iFid[k].imag - 1j*iFid[k].real
                rFid[k] = rFid[k].imag + 1j*iFid[k].imag
                
            fid = np.copy(rFid)
        
        return fid
        # end quad2D

    def readPipe2D(self, pName, fName):
        npd = nmrpipeData.NmrPipeData()
        npd.readPipe(pName, fName)
        self.spc = npd.spc
        self.proc.sw_h[0] = npd.fdf2sw
        self.refShift[0]  = npd.fdf2orig/self.acq.sfo1
        self.refPoint[0]  = 0
        self.proc.sw_h[1] = npd.fdf1sw
        self.refShift[1]  = npd.fdf1orig/self.acq.sfo2
        self.refPoint[1]  = 0
        self.calcPPM()
        # end readPipe2D

    def readSpc(self):
        self.acq.read(self.dataSetName + os.sep + self.dataSetNumber)
        self.proc.read(self.dataSetName + os.sep + self.dataSetNumber)
        self.origDataSet = self.dataSetName + self.dataSetNumber
        titleFile    = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'pdata' + os.sep + '1' + os.sep + 'title'
        pulProgFile = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'pulseprogram'
        if(os.path.isfile(titleFile)):
            fid = open(titleFile,"r")
            self.title = fid.read()
            fid.close()
        
        if(os.path.isfile(pulProgFile)):
            fid = open(pulProgFile,"r")
            self.pulseProgram = fid.read()
            fid.close()
        
        fidFile1   = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'fid'
        spcFile1r  = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'pdata' + os.sep + '1' + os.sep + '1r'
        spcFile1i  = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'pdata' + os.sep + '1' + os.sep + '1i'
        fidFile2   = self.dataSetName + os.sep + self.dataSetNumber + os.sep + 'ser'
        if(os.path.isfile(fidFile1)):
            #read 1D FID file
            self.fid.resize(1,int(self.acq.nDataPoints[0]/2))
            f        = open(fidFile1,'rb')
            fid = np.fromfile(f, dtype = np.int32)
            f.close()
            self.fid[0].real = fid[0::2]
            self.fid[0].imag = -fid[1::2]
            self.dim         = 1
            
        if(os.path.isfile(spcFile1r)):
            #read 1D spectrum file (real part)
            f                = open(spcFile1r,'rb')
            fid              = np.fromfile(f, dtype = np.int32)
            self.spc.resize(1,int(len(fid)))
            self.spc[0].real = fid;
            f.close()
            if(os.path.isfile(spcFile1i)):
                #read 1D spectrum file (imaginary part)
                f                = open(spcFile1i,'rb')
                fid              = np.fromfile(f, dtype = np.int32)
                f.close()
                self.spc[0].imag = fid
            
        
        elif(os.path.isfile(fidFile2)):
            #read 2D spectrum
            np1 = int(self.acq.nDataPoints[0])
            np2 = int(self.acq.nDataPoints[1])
            self.fid.resize(int(np2),int(np1/2))
            f   = open(fidFile2,'rb')
            fid = np.fromfile(f, dtype = np.int32)
            f.close()
            fid = fid.reshape(int(np2),int(np1))
            for x in range(np2):
                self.fid[x].real = fid[x][0::2]
                self.fid[x].imag = -fid[x][1::2]
                
            self.dim = 2
        
        self.proc.sw_h = np.copy(self.acq.sw_h)
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
        x     = np.linspace(0, len(fid)-1,len(fid))
        f0    = np.copy(fid)
        fid  = np.roll(fid, math.floor(-self.acq.groupDelay))
        pp   = np.polyfit(x, fid, self.proc.polyOrder)
        fid  = fid - np.polyval(pp, x)
        fid  = np.roll(fid, math.ceil(self.acq.groupDelay))
        return fid
        # end smo
    
    def waterSupp(self, fid):
        if(self.proc.waterSuppression==1):
            fid = self.conv(fid)
        
        if(self.proc.waterSuppression==2):
            fid = self.smo(fid)            
            
        return fid
        # end waterSupp
    
    def zeroFill(self,fid, dim=0):
        fid1                 = np.zeros(self.proc.nPoints[dim], dtype = 'complex')
        fid1[:int(len(fid))] = fid
        return fid1
        # end zeroFill
    
