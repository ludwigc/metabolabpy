import numpy as np
import os
import pickle

from metabolabpy.nmr import nmrData as nd
from metabolabpy.nmr import nmrPreProc as npp


class NmrDataSet:

    def __init__(self):
        self.nmrdat   = [[]]
        self.s        = 0
        self.e        = -1
        self.pp       = npp.NmrPreProc()
        self.deselect  = np.array([])
        self.deselect2 = np.array([])
        self.cmdBuffer = np.array([])
        self.cmdIdx    = -1
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
        
    def bucketSpectra(self):
        idx1 = np.arange(len(self.nmrdat[self.s][0].ppm1))
        idx2 = idx1[::int(self.pp.bucketPoints)]
        idx2 = np.append(idx2, len(idx1))
        ppm  = np.array([])
        for k in range(len(idx2)-1):
            ppm = np.append(ppm, np.mean(self.nmrdat[self.s][0].ppm1[idx2[k]:idx2[k+1]]))
            
        for k in range(len(self.nmrdat[self.s])):
            spc = np.array([])
            for l in range(len(idx2)-1):
                spc = np.append(spc, np.sum(self.nmrdat[self.s][k].spc[0][idx2[l]:idx2[l+1]].real/self.pp.bucketPoints))
            
            self.nmrdat[self.s][k].ppm1   = ppm
            self.nmrdat[self.s][k].spc    = np.resize(self.nmrdat[self.s][k].spc,(1,len(spc)))
            self.nmrdat[self.s][k].spc[0] = np.copy(spc)
        
        # end bucketSpectra
        
    def dataPreProcessing(self):
        self.ftAll()
        self.baseline1dAll()
        self.autorefAll()
        self.shiftRef()
        self.noiseFilteringInit()
        self.deselect  = np.zeros(len(self.nmrdat[self.s][0].spc[0]))
        self.deselect2 = np.zeros(len(self.nmrdat[self.s][0].spc[0]))
        if(self.pp.flagExcludeRegion == True):
            self.excludeRegion()
            
        if(self.pp.flagSegmentalAlignment == True):
            print("segmentalAlignment")
        
        if(self.pp.flagNoiseFiltering == True):
            self.noiseFiltering()
            
        idx  = np.where(self.deselect  == 1)
        idx2 = np.where(self.deselect2 == len(self.nmrdat[self.s]))
        for k in range(len(self.nmrdat[self.s])):
            self.nmrdat[self.s][k].spc[0][idx]  = np.zeros(len(idx))
            self.nmrdat[self.s][k].spc[0][idx2] = np.zeros(len(idx2))
        
        if(self.pp.flagBucketSpectra == True):
            self.bucketSpectra()
        
        if(self.pp.flagCompressBuckets == True):
            print("compressBuckets")
            
        if(self.pp.flagScaleSpectra == True):
            print("scaleSpectra")
        
        if(self.pp.flagVarianceStabilisation == True):
            print("varianceStabilisation")
            
        if(self.pp.flagExportDataSet == True):
            self.exportDataSet()
        
        # end dataPreProcessing
        
    def excludeRegion(self):
        if(len(self.pp.excludeStart) != len(self.pp.excludeEnd)):
            return
        
        for k in range(len(self.pp.excludeStart)):
            idx                                = np.where((self.nmrdat[self.s][0].ppm1 > self.pp.excludeStart[k]) & (self.nmrdat[self.s][0].ppm1 < self.pp.excludeEnd[k]))
            self.deselect[idx]                 = np.ones(len(idx))
        
        # end excludeRegion
        
    def exportDataSet(self):
        fName = os.path.join(self.pp.exportPathName, self.pp.exportFileName)
        try:
            f = open(fName,'w')
        except:
            return
        
        if(self.pp.exportDelimiterTab == True):
            delim = '\t'
        else:
            delim = self.pp.exportCharacter
            
        spc = np.zeros(len(self.nmrdat[self.s][0].spc[0]))
        for k in range(len(self.nmrdat[self.s])):
            spc += self.nmrdat[self.s][k].spc[0].real
        
        deselect      = np.zeros(len(self.nmrdat[self.s][0].spc[0]))
        idx           = np.where(spc == 0)
        deselect[idx] = np.ones(len(idx))
        if(self.pp.exportSamplesInRowsCols == 0): # samples in rows
            f.write("ppm" + delim + " ")
            for k in range(len(self.nmrdat[self.s][0].ppm1)):
                if(deselect[k] == 0):
                    f.write(delim + str(self.nmrdat[self.s][0].ppm1[k]))
                
                
            f.write("\n")
            for k in range(len(self.nmrdat[self.s])):
                dse = os.path.split(self.nmrdat[self.s][k].origDataSet)
                ds  = os.path.split(dse[0])
                f.write(ds[1] + " " + dse[1] + delim + self.pp.classSelect[k])
                for l in range(len(self.nmrdat[self.s][k].spc[0])):
                    if(deselect[l] == 0):
                        f.write(delim + str(self.nmrdat[self.s][k].spc[0][l].real))
                        
                f.write("\n")
                
        else:                                     # samples in cols
            f.write("ppm")
            for k in range(len(self.nmrdat[self.s])):
                dse = os.path.split(self.nmrdat[self.s][k].origDataSet)
                ds  = os.path.split(dse[0])
                f.write(delim + ds[1] + " " + dse[1])
                
            f.write("\n")
            f.write(" ")
            for k in range(len(self.pp.classSelect)):
                f.write(delim + self.pp.classSelect[k])
                
            f.write("\n")
            for k in range(len(self.nmrdat[self.s][0].spc[0])):
                if(deselect[k] == 0):
                    f.write(str(self.nmrdat[self.s][0].ppm1[k]))
                    for l in range(len(self.nmrdat[self.s])):
                        f.write(delim + str(self.nmrdat[self.s][l].spc[0][k].real))

                    f.write("\n")

            
        
        f.close()
        # end exportDataSet
        
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
        val    = self.pp.noiseThreshold*self.pp.stdVal
        for k in range(len(self.nmrdat[self.s])):
            idx = np.where(self.nmrdat[self.s][k].spc[0].real < val)
            self.deselect2[idx] += np.ones(len(idx))
            
        # end noiseFiltering
        
    def noiseFilteringInit(self):
        spcIdx         = np.where((self.nmrdat[self.s][0].ppm1>self.pp.noiseStart) & (self.nmrdat[self.s][0].ppm1<self.pp.noiseEnd))
        self.pp.stdVal = np.std(self.nmrdat[self.s][0].spc[0][spcIdx].real)
        # end noiseFilteringInit
        
    def preProcInit(self):
        self.pp.init(len(self.nmrdat[self.s]))
        # end preProcInit
        
    def readSpc(self, dataSetName, dataSetNumber):
        self.e            = len(self.nmrdat[self.s])
        nd1               = nd.NmrData()
        nd1.dataSetName   = dataSetName
        nd1.dataSetNumber = dataSetNumber
        nd1.readSpc()
        self.nmrdat[self.s].append(nd1)
        # end readSpc
        
    def readSpcs(self, dataPath, dataExp):
        for k in range(len(dataExp)):
            self.readSpc(dataPath, str(dataExp[k]))
    
    # end readSpcs

    def load(self,dataSetName):
        dataSets = np.array([])
        lDir     = os.listdir(dataSetName)
        curDataNotFound = True
        for k in range(len(lDir)):
            if(lDir[k] == 'curPars.dat'):
                curDataNotFound = False
            
        
        if(curDataNotFound == True):
            return
        
        fName          = os.path.join(dataSetName, 'curPars.dat')
        f              = open(fName, 'rb')
        curPars        = pickle.load(f)
        f.close()
        self.s         = curPars[0]
        self.e         = curPars[1]
        self.pp        = curPars[2]
        self.deselect  = curPars[3]
        self.deselect2 = curPars[4]
        self.cmdBuffer = curPars[5]
        self.cmdIdx    = curPars[6]
        for k in range(len(lDir)):
            if(os.path.isdir(os.path.join(dataSetName,lDir[k]))):
                dataSets = np.append(dataSets, lDir[k])
            
        
        dataSets = np.sort(dataSets)
        dataSetExps = []
        for k in range(len(dataSets)):
            dirName  = os.listdir(os.path.join(dataSetName,dataSets[k]))
            dataExps = np.array([])
            for l in range(len(dirName)):
                if(os.path.isdir(os.path.join(os.path.join(dataSetName,dataSets[k]), dirName[l]))):
                    if(os.path.isfile(os.path.join(os.path.join(os.path.join(dataSetName,dataSets[k]),dirName[l]),'titleFile.txt')) and
                       os.path.isfile(os.path.join(os.path.join(os.path.join(dataSetName,dataSets[k]),dirName[l]),'acqusText.txt')) and
                       os.path.isfile(os.path.join(os.path.join(os.path.join(dataSetName,dataSets[k]),dirName[l]),'acqu2sText.txt')) and
                       os.path.isfile(os.path.join(os.path.join(os.path.join(dataSetName,dataSets[k]),dirName[l]),'acqu3sText.txt')) and
                       os.path.isfile(os.path.join(os.path.join(os.path.join(dataSetName,dataSets[k]),dirName[l]),'procsText.txt')) and
                       os.path.isfile(os.path.join(os.path.join(os.path.join(dataSetName,dataSets[k]),dirName[l]),'proc2sText.txt')) and
                       os.path.isfile(os.path.join(os.path.join(os.path.join(dataSetName,dataSets[k]),dirName[l]),'proc3sText.txt')) and
                       os.path.isfile(os.path.join(os.path.join(os.path.join(dataSetName,dataSets[k]),dirName[l]),'nmrDataSet.dat'))):
                        dataExps = np.append(dataExps, dirName[l])
                    
                
            dataExps = np.sort(dataExps)
            dataSetExps.append(dataExps)
        
        self.nmrdat = []
        for k in range(len(dataSetExps)):
            self.nmrdat.append([])
            for l in range(len(dataSetExps[k])):
                fName = os.path.join(os.path.join(os.path.join(dataSetName,dataSets[k]),dataSetExps[k][l]),'nmrDataSet.dat')
                f     = open(fName, 'rb')
                nd    = pickle.load(f)
                f.close()
                self.nmrdat[k].append(nd)
                fName = os.path.join(os.path.join(os.path.join(dataSetName,dataSets[k]),dataSetExps[k][l]),'titleFile.txt')
                f     = open(fName, 'r')
                self.nmrdat[k][l].title = f.read()
                f.close()
                fName = os.path.join(os.path.join(os.path.join(dataSetName,dataSets[k]),dataSetExps[k][l]),'acqusText.txt')
                f     = open(fName, 'r')
                self.nmrdat[k][l].acq.acqusText = f.read()
                f.close()
                fName = os.path.join(os.path.join(os.path.join(dataSetName,dataSets[k]),dataSetExps[k][l]),'acqu2sText.txt')
                f     = open(fName, 'r')
                self.nmrdat[k][l].acq.acqu2sText = f.read()
                f.close()
                fName = os.path.join(os.path.join(os.path.join(dataSetName,dataSets[k]),dataSetExps[k][l]),'acqu3sText.txt')
                f     = open(fName, 'r')
                self.nmrdat[k][l].acq.acqu3sText = f.read()
                f.close()
                fName = os.path.join(os.path.join(os.path.join(dataSetName,dataSets[k]),dataSetExps[k][l]),'procsText.txt')
                f     = open(fName, 'r')
                self.nmrdat[k][l].proc.procsText = f.read()
                f.close()
                fName = os.path.join(os.path.join(os.path.join(dataSetName,dataSets[k]),dataSetExps[k][l]),'proc2sText.txt')
                f     = open(fName, 'r')
                self.nmrdat[k][l].proc.proc2sText = f.read()
                f.close()
                fName = os.path.join(os.path.join(os.path.join(dataSetName,dataSets[k]),dataSetExps[k][l]),'proc3sText.txt')
                f     = open(fName, 'r')
                self.nmrdat[k][l].proc.proc2sText = f.read()
                f.close()
            
        
        # end load
        
    def save(self,dataSetName):
        if(len(dataSetName) == 0):
            return
        
        try:
            os.makedirs(dataSetName)
        except:
            pass
        
        fName = os.path.join(dataSetName,'curPars.dat')
        f     = open(fName,'wb')
        pickle.dump([self.s,self.e,self.pp,self.deselect,self.deselect2,self.cmdBuffer,self.cmdIdx],f)
        f.close()
        for k in range(len(self.nmrdat)):
            setPath = os.path.join(dataSetName,str(k+1))
            try:
                os.makedirs(setPath)
            except:
                pass
            
            for l in range(len(self.nmrdat[k])):
                expPath = os.path.join(setPath,str(l+1))
                try:
                    os.makedirs(expPath)
                except:
                    pass
                
                fName   = os.path.join(expPath,'nmrDataSet.dat')
                f       = open(fName,'wb')
                pickle.dump(self.nmrdat[k][l],f)
                f.close()
                fName   = os.path.join(expPath,'titleFile.txt')
                f       = open(fName,'w')
                f.write(self.nmrdat[k][l].title)
                f.close()
                fName   = os.path.join(expPath, 'procsText.txt')
                f       = open(fName,'w')
                f.write(self.nmrdat[k][l].proc.procsText)
                f.close()
                fName   = os.path.join(expPath, 'proc2sText.txt')
                f       = open(fName,'w')
                f.write(self.nmrdat[k][l].proc.proc2sText)
                f.close()
                fName   = os.path.join(expPath, 'proc3sText.txt')
                f       = open(fName,'w')
                f.write(self.nmrdat[k][l].proc.proc3sText)
                f.close()
                fName   = os.path.join(expPath, 'acqusText.txt')
                f       = open(fName,'w')
                f.write(self.nmrdat[k][l].acq.acqusText)
                f.close()
                fName   = os.path.join(expPath, 'acqu2sText.txt')
                f       = open(fName,'w')
                f.write(self.nmrdat[k][l].acq.acqu2sText)
                f.close()
                fName   = os.path.join(expPath, 'acqu3sText.txt')
                f       = open(fName,'w')
                f.write(self.nmrdat[k][l].acq.acqu3sText)
                f.close()
        # end save
        
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
        
