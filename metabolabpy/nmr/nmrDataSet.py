import numpy as np

from metabolabpy.nmr import nmrData as nd
from metabolabpy.nmr import nmrPreProc as npp


class NmrDataSet:

    def __init__(self):
        self.nmrdat   = [[]]
        self.s        = 0
        self.e        = -1
        self.pp       = npp.NmrPreProc()
        self.deselect = np.array([])
        # end __init__

    def __str__(self):
        mmax     = 0
        for k in range(len(self.nmrdat)):
            mmax = max(mmax,len(self.nmrdat[k]))
            
        rString  = 'MetaboLabPy NMR Data Set (v. 0.1)\n'
        rString += '__________________________________________________________________\n'
        rString += ' Number of data sets\t\t\t: {:0.0f}\n'.format(len(self.nmrdat))
        rString += ' max number of NMR spectra per data set\t: {:0.0f}\n'.format(mmax)
        if(mmax>0):
            rString += ' Current data set/exp\t\t\t: {:0.0f}/{:0.0f}\n'.format(self.s + 1,self.e + 1)
            rString += '__________________________________________________________________\n'
            rString += 'Current title file: \n'
            rString += self.nmrdat[self.s][self.e].title
        return rString
        # end __str__

    def autobaseline1d(self):
        if(len(self.nmrdat)>0):
            if(len(self.nmrdat[self.s])>0):
                self.nmrdat[self.s][self.e].autobaseline1d()
                
                
        # end autobaseline1d

    def autobaseline1dAll(self):
        nExp    = len(self.nmrdat[self.s])
        origExp = self.e
        for k in range(nExp):
            self.e = k
            self.autobaseline1d()
    
        self.e = origExp
        return "Finished autobaseline1dAll"
        # end autobaseline1dAll

    def autophase1d(self):
        if(len(self.nmrdat)>0):
            if(len(self.nmrdat[self.s])>0):
                self.nmrdat[self.s][self.e].autophase1d()
                
                
        # end autophase1d

    def autophase1dAll(self):
        nExp    = len(self.nmrdat[self.s])
        origExp = self.e
        for k in range(nExp):
            self.e = k
            self.autophase1d()
    
        self.e = origExp
        return "Finished autophase1dAll"
        # end autophase1dAll

    def autoref(self,tmsp=True):
        if(self.nmrdat[self.s][self.e].dim == 1):
            self.nmrdat[self.s][self.e].autoRef(tmsp)
            #self.nmrdat[self.s][self.e].setRef(np.array([0.0]), np.array([14836]))
            
        else:
            refShifts                = np.array([1.33, 19.3])
            refPoints                = np.array([262, 2150])
            self.nmrdat[self.s][self.e].setRef(refShifts, refPoints)
        
        return "Finished autoref"
        # end autoref

    def autorefAll(self):
        nExp    = len(self.nmrdat[self.s])
        origExp = self.e
        for k in range(nExp):
            self.e = k
            self.autoref()
        
        self.e = origExp
        return "Finished autorefAll"
        # end autorefAll
        
    def baseline1d(self):
        if(self.nmrdat[self.s][self.e].dim == 1):
            self.nmrdat[self.s][self.e].baseline1d()
            
        # end baseline1d
        
    def baseline1dAll(self):
        origExp = self.e
        for k in range(len(self.nmrdat[self.s])):
            self.e = k
            self.baseline1d()
            
        self.e = origExp
        # end baseline1dAll
        
    def dataPreProcessing(self):
        self.ftAll()
        self.baseline1dAll()
        self.autorefAll()
        self.shiftRef()
        self.deselect = np.zeros(len(self.nmrdat[self.s][0].spc[0]))
        print(sum(self.deselect))
        if(self.pp.flagExcludeRegion == True):
            self.excludeRegion()
            
        print(sum(self.deselect))
        if(self.pp.flagSegmentalAlignment == True):
            print("segmentalAlignment")
        
        print(sum(self.deselect))
        if(self.pp.flagNoiseFiltering == True):
            self.noiseFiltering()
            
        print(sum(self.deselect))
        idx = np.where(self.deselect == 1)
        for k in range(len(self.nmrdat[self.s])):
            #print(idx)
            self.nmrdat[self.s][k].spc[0][idx].real = np.zeros(len(idx))
        
        print(sum(self.deselect))
        if(self.pp.flagBucketSpectra == True):
            print("bucketSpectra")
        
        if(self.pp.flagCompressBuckets == True):
            print("compressBuckets")
            
        if(self.pp.flagScaleSpectra == True):
            print("scaleSpectra")
        
        if(self.pp.flagVarianceStabilisation == True):
            print("varianceStabilisation")
            
        if(self.pp.flagExportDataSet == True):
            print("exportDataSet")
        
        # end dataPreProcessing
        
    def excludeRegion(self):
        if(len(self.pp.excludeStart) != len(self.pp.excludeEnd)):
            return
        
        for k in range(len(self.pp.excludeStart)):
            idx                                = np.where((self.nmrdat[self.s][0].ppm1 > self.pp.excludeStart[k]) & (self.nmrdat[self.s][0].ppm1 < self.pp.excludeEnd[k]))
            self.deselect[idx]                 = np.ones(len(idx))
        
        # end excludeRegion
        
    def ft(self):
        if(self.nmrdat[self.s][self.e].dim == 1):
            self.nmrdat[self.s][self.e].procSpc1D()
            #self.nmrdat[self.s][self.e].setRef(np.array([0.0]), np.array([14836]))
            
        else:
            self.nmrdat[self.s][self.e].procSpc()
            refShifts                = np.array([1.33, 19.3])
            refPoints                = np.array([262, 2150])
            self.nmrdat[self.s][self.e].setRef(refShifts, refPoints)
        
        # end ft

    def ftAll(self):
        nExp    = len(self.nmrdat[self.s])
        origExp = self.e
        for k in range(nExp):
            self.e = k
            self.ft()
        
        self.e = origExp
        return "Finished ftAll"
        # end ftAll
        
    def noiseFiltering(self):
        spcIdx = np.where((self.nmrdat[self.s][0].ppm1>self.pp.noiseStart) & (self.nmrdat[self.s][0].ppm1<self.pp.noiseEnd))
        stdVal = np.std(self.nmrdat[self.s][0].spc[0][spcIdx].real)
        val    = self.pp.noiseThreshold*stdVal
        for k in range(len(self.nmrdat[self.s])):
            idx = np.where(self.nmrdat[self.s][k].spc[0].real < val)
            self.deselect[idx] = np.ones(len(idx))
            
        # end noiseFiltering
        
    def preProcInit(self):
        self.pp.init(len(self.nmrdat[self.s]))
        # end preProcInit
        
    def readSpc(self, dataSetName, dataSetNumber):
        self.e = len(self.nmrdat[self.s])
        nd1    = nd.NmrData()
        nd1.dataSetName   = dataSetName
        nd1.dataSetNumber = dataSetNumber
        nd1.readSpc()
        self.nmrdat[self.s].append(nd1)
        # end readSpc
        
    def readSpcs(self, dataPath, dataExp):
        for k in range(len(dataExp)):
            self.readSpc(dataPath, str(dataExp[k]))
    
    # end readSpcs

    def setGb(self, gb):
        nExp = len(self.nmrdat[self.s])
        for k in range(nExp):
            for l in range(len(gb)):
                self.nmrdat[self.s][k].proc.gb[l] = gb[l]
            
        
    
        return "setGb"
    # end setGb
        
    def setLb(self, lb):
        nExp = len(self.nmrdat[self.s])
        for k in range(nExp):
            for l in range(len(lb)):
                self.nmrdat[self.s][k].proc.lb[l] = lb[l]
            
        
    
        return "setLb"
    # end setLb
        
    def setPhFromExp(self, exp = -1):
        if(exp == -1):
            exp = self.e
            
        nExp = len(self.nmrdat[self.s])
        for k in range(nExp):
            if(k != exp):
                self.nmrdat[self.s][k].proc.ph0 = self.nmrdat[self.s][exp].proc.ph0
                self.nmrdat[self.s][k].proc.ph1 = self.nmrdat[self.s][exp].proc.ph1
                
        
        return "setPhFromExp"
    # end setPhFromExp
        
    def setPh0(self, ph0):
        nExp = len(self.nmrdat[self.s])
        for k in range(nExp):
            for l in range(len(ph0)):
                self.nmrdat[self.s][k].proc.ph0[l] = ph0[l]
            
        
    
        return "setPh0"
    # end setPh0
        
    def setPh1(self, ph1):
        nExp = len(self.nmrdat[self.s])
        for k in range(nExp):
            for l in range(len(ph1)):
                self.nmrdat[self.s][k].proc.ph1[l] = ph1[l]
            
        
    
        return "setPh1"
    # end setPh1
        
    def setSsb(self, ssb):
        nExp = len(self.nmrdat[self.s])
        for k in range(nExp):
            for l in range(len(ssb)):
                self.nmrdat[self.s][k].proc.ssb[l] = ssb[l]
            
        
    
        return "setSsb"
    # end setSsb
        
    def setWindowType(self, wt):
        nExp = len(self.nmrdat[self.s])
        for k in range(nExp):
            for l in range(len(wt)):
                self.nmrdat[self.s][k].proc.windowType[l] = wt[l]
            
        
    
        return "setWindowType"
    # end setWindowType
        
    def setZeroFill(self, zf):
        nExp = len(self.nmrdat[self.s])
        for k in range(nExp):
            for l in range(len(zf)):
                self.nmrdat[self.s][k].proc.nPoints[l] = zf[l]
            
        
    
        return "setZeroFill"
    # end setZeroFill
    
    def shiftRef(self):
        for k in range(len(self.nmrdat[self.s])):
            shiftDelta = self.nmrdat[self.s][k].refPoint[0] - self.nmrdat[self.s][0].refPoint[0]
            self.nmrdat[self.s][k].spc[0] = np.roll(self.nmrdat[self.s][k].spc[0],shiftDelta)
            self.nmrdat[self.s][k].ppm1   = self.nmrdat[self.s][0].ppm1
            self.nmrdat[self.s][k].refShift = self.nmrdat[self.s][0].refShift
    # end shiftRef
        
