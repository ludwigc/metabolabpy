import numpy as np
import os
import pickle

from metabolabpy.nmr import nmrData as nd
from metabolabpy.nmr import nmrPreProc as npp
import matplotlib.pyplot as pl
import math

class NmrDataSet:

    def __init__(self):
        self.nmrdat            = [[]]
        self.s                 = 0
        self.e                 = -1
        self.pp                = npp.NmrPreProc()
        self.deselect          = np.array([])
        self.deselect2         = np.array([])
        self.cmdBuffer         = np.array([])
        self.cmdIdx            = -1
        self.script            = ""
        self.console           = ""
        self.fileFormatVersion = 0.1
        # end __init__

    def __str__(self): # pragma: no cover
        rString  = 'MetaboLabPy NMR Data Set (v. 0.1)\n'
        rString += '__________________________________________________________________\n'
        rString += ' Number of data sets\t\t\t: {:0.0f}\n'.format(len(self.nmrdat))
        rString += ' Number of NMR spectra in current data set\t: {:0.0f}\n'.format(len(self.nmrdat[self.s]))
        if(len(self.nmrdat[self.s])>0):
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
            self.nmrdat[self.s][self.e].autoRef()
            #refShifts                = np.array([1.33, 19.3])
            #refPoints                = np.array([262, 2150])
            #self.nmrdat[self.s][self.e].setRef(refShifts, refPoints)
        
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
        
    def clear(self):
        self.nmrdat = [[]]
        self.s = 0
        self.e = -1
        return "Workspace cleared"
        # end clear

    def dataPreProcessing(self):
        if(self.nmrdat[self.s][0].projectedJres == False):
            self.ftAll()
            self.baseline1dAll()
            self.autorefAll()
            self.shiftRef()

        else:
            s = self.s
            e = self.e
            self.s = self.nmrdat[s][e].origJresSet
            self.e = self.nmrdat[s][e].origJresExp
            self.pjres(s+1, self.nmrdat[s][e].pjresMode)
            self.s = s
            self.e = e

        self.noiseFilteringInit()
        self.deselect  = np.zeros(len(self.nmrdat[self.s][0].spc[0]))
        self.deselect2 = np.zeros(len(self.nmrdat[self.s][0].spc[0]))
        if(self.pp.flagExcludeRegion == True):
            self.excludeRegion()
            
        if(self.pp.flagSegmentalAlignment == True):
            self.segmentalAlignment()
        
        if(self.pp.flagNoiseFiltering == True):
            self.noiseFiltering()
            
        idx  = np.where(self.deselect  == 1)
        idx2 = np.where(self.deselect2 == len(self.nmrdat[self.s]))
        for k in range(len(self.nmrdat[self.s])):
            self.nmrdat[self.s][k].spc[0][idx]  = np.zeros(len(idx))
            self.nmrdat[self.s][k].spc[0][idx2] = np.zeros(len(idx2))
        
        mVal = 0
        for k in range(len(self.nmrdat[self.s])):
            mVal = min(mVal,np.min(self.nmrdat[self.s][k].spc[0].real))

        if mVal < 0 and self.pp.rDolphinExport is False:
            for k in range(len(self.nmrdat[self.s])):
                self.nmrdat[self.s][k].spc[0] -= mVal

        
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
        if os.path.isdir(self.pp.exportPathName) is False:
            os.makedirs(self.pp.exportPathName)

        fName = os.path.join(self.pp.exportPathName, self.pp.exportFileName)
        f = open(fName, 'w')
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
            if self.pp.rDolphinExport is False:
                f.write("ppm" + delim + " ")
                for k in range(len(self.nmrdat[self.s][0].ppm1)):
                    if(deselect[k] == 0):
                        f.write(delim + str(self.nmrdat[self.s][0].ppm1[k]))


            else:
                f.write(str(self.nmrdat[self.s][0].ppm1[0]))
                for k in range(1,len(self.nmrdat[self.s][0].ppm1)):
                    if (deselect[k] == 0):
                        f.write(delim + str(self.nmrdat[self.s][0].ppm1[k]))


            f.write("\n")
            for k in range(len(self.nmrdat[self.s])):
                dse = os.path.split(self.nmrdat[self.s][k].origDataSet)
                ds  = os.path.split(dse[0])
                if self.pp.rDolphinExport is False:
                    f.write(ds[1] + " " + dse[1] + delim + self.pp.classSelect[k])
                    for l in range(len(self.nmrdat[self.s][k].spc[0])):
                        if(deselect[l] == 0):
                            f.write(delim + str(self.nmrdat[self.s][k].spc[0][l].real))


                else:
                    f.write(str(self.nmrdat[self.s][k].spc[0][0].real))
                    for l in range(1,len(self.nmrdat[self.s][k].spc[0])):
                        if (deselect[l] == 0):
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
        if self.pp.rDolphinExport is True:
            fName = os.path.join(self.pp.exportPathName, "Parameters.csv")
            f = open(fName, 'w')
            f.write("Parameter,Value\n")
            f.write("nmr folder path,\n")
            f.write("1D data index,\n")
            f.write("proc_no,\n")
            f.write("spectra dataset path (csv format)," + self.pp.exportPathName + "/" + self.pp.exportFileName + "\n")
            f.write("Metadata path (csv format)," + self.pp.exportPathName + "/Metadata.csv\n")
            f.write("ROI patters file," + self.pp.exportPathName + "/ROI_profile.csv\n")
            f.write("Normalization (0=No;1=Eretic; 2=TSP; 3=Creatinine; 4=Spectra Sum; 5=PQN),2\n")
            f.write("Alignment (0=No;1=Glucose; 2=TSP; 3=Formate),2\n")
            f.write("Suppression,12-9.5;6.1-5.6;5.1-4.5\n")
            f.write("Spectrometer Frequency (MHz)," + str(self.nmrdat[self.s][0].acq.sfo1) + "\n")
            f.write("Bucket resolution," + str(self.pp.bucketPPM) + "\n")
            f.write("Biofluid,Urine\n")
            f.write("2D-Path,\n")
            f.write("Specific program parameters,\n")
            f.close()
            fName = os.path.join(self.pp.exportPathName, "Metadata.csv")
            f = open(fName, 'w')
            f.write("Sample,Individual,Sample Type\n")
            for k in range(len(self.nmrdat[self.s])):
                f.write(os.path.split(self.nmrdat[self.s][k].origDataSet)[1] + "," + str(k+1) + ",1\n")

            f.close()
    # end exportDataSet
        
    def ft(self):
        if(self.nmrdat[self.s][self.e].dim == 1):
            self.nmrdat[self.s][self.e].procSpc1D()

        else:
            self.nmrdat[self.s][self.e].procSpc()

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

    def load(self, dataSetName):
        dataSets = np.array([])
        lDir = os.listdir(dataSetName)
        curDataNotFound = True
        for k in range(len(lDir)):
            if (lDir[k] == 'curPars.dat'):
                curDataNotFound = False

        if (curDataNotFound == True):
            return

        fName = os.path.join(dataSetName, 'curPars.dat')
        f = open(fName, 'rb')
        curPars = pickle.load(f)
        f.close()
        self.fileFormatVersion = curPars[0]
        self.s = curPars[1]
        self.e = curPars[2]
        #self.pp = curPars[3]
        c = curPars[3]
        if hasattr(c, 'excludeStart'):
            self.pp.excludeStart = c.excludeStart

        if hasattr(c, 'excludeEnd'):
            self.pp.excludeEnd = c.excludeEnd

        if hasattr(c, 'segStart'):
            self.pp.segStart = c.segStart

        if hasattr(c, 'segEnd'):
            self.pp.segEnd = c.segEnd

        if hasattr(c, 'noiseThreshold'):
            self.pp.noiseThreshold = c.noiseThreshold

        if hasattr(c, 'noiseStart'):
            self.pp.noiseStart = c.noiseStart

        if hasattr(c, 'noiseEnd'):
            self.pp.noiseEnd = c.noiseEnd

        if hasattr(c, 'bucketPoints'):
            self.pp.bucketPoints = c.bucketPoints

        if hasattr(c, 'bucketPPM'):
            self.pp.bucketPPM = c.bucketPPM

        if hasattr(c, 'compressBuckets'):
            self.pp.compressBuckets = c.compressBuckets

        if hasattr(c, 'scaleSpc'):
            self.pp.scaleSpc= c.scaleSpc

        if hasattr(c, 'varianceStabilisation'):
            self.pp.varianceStabilisation = c.varianceStabilisation

        if hasattr(c, 'varLambda'):
            self.pp.varLambda = c.varLambda

        if hasattr(c, 'varX0'):
            self.pp.varX0 = c.varX0

        if hasattr(c, 'pName'):
            self.pp.pName = c.pName

        if hasattr(c, 'fName'):
            self.pp.fName = c.fName

        if hasattr(c, 'samplesInRows'):
            self.pp.samplesInRows = c.samplesInRows

        if hasattr(c, 'plotSelect'):
            self.pp.plotSelect = c.plotSelect

        if hasattr(c, 'classSelect'):
            self.pp.classSelect = c.classSelect

        if hasattr(c, 'plotColours'):
            self.pp.plotColours = c.plotColours

        if hasattr(c, 'preProcFill'):
            self.pp.preProcFill = c.preProcFill

        if hasattr(c, 'alpha'):
            self.pp.alpha = c.alpha

        if hasattr(c, 'colour'):
            self.pp.colour = c.colour

        if hasattr(c, 'thColour'):
            self.pp.thColour = c.thColour

        if hasattr(c, 'thLineWidth'):
            self.pp.thLineWidth = c.thLineWidth

        if hasattr(c, 'flagExcludeRegion'):
            self.pp.flagExcludeRegion = c.flagExcludeRegion

        if hasattr(c, 'flagSegmentalAlignment'):
            self.pp.flagSegmentalAlignment = c.flagSegmentalAlignment

        if hasattr(c, 'flagNoiseFiltering'):
            self.pp.flagNoiseFiltering = c.flagNoiseFiltering

        if hasattr(c, 'flagBucketSpectra'):
            self.pp.flagBucketSpectra = c.flagBucketSpectra

        if hasattr(c, 'flagCompressBuckets'):
            self.pp.flagCompressBuckets = c.flagCompressBuckets

        if hasattr(c, 'flagScaleSpectra'):
            self.pp.flagScaleSpectra = c.flagScaleSpectra

        if hasattr(c, 'flagVarianceStabilisation'):
            self.pp.flagVarianceStabilisation = c.flagVarianceStabilisation

        if hasattr(c, 'flagExportDataSet'):
            self.pp.flagExportDataSet = c.flagExportDataSet

        if hasattr(c, 'stdVal'):
            self.pp.stdVal = c.stdVal

        if hasattr(c, 'exportPathName'):
            self.pp.exportPathName= c. exportPathName

        if hasattr(c, 'exportFileName'):
            self.pp.exportFileName = c.exportFileName

        if hasattr(c, 'exportDelimiterTab'):
            self.pp.exportDelimiterTab = c.exportDelimiterTab

        if hasattr(c, 'exportCharacter'):
            self.pp.exportCharacter = c. exportCharacter

        if hasattr(c, 'exportSamplesInRowsCols'):
            self.pp.exportSamplesInRowsCols = c. exportSamplesInRowsCols

        if hasattr(c, 'rDolphinExport'):
            self.pp.rDolphinExport = c.rDolphinExport

        self.deselect = curPars[4]
        self.deselect2 = curPars[5]
        self.cmdBuffer = curPars[6]
        self.cmdIdx = curPars[7]
        self.script = curPars[8]
        self.console = curPars[9]
        for k in range(len(lDir)):
            if (os.path.isdir(os.path.join(dataSetName, lDir[k]))):
                dataSets = np.append(dataSets, lDir[k])

        dataSets = np.sort(dataSets)
        dataSetExps = []
        for k in range(len(dataSets)):
            dirName = os.listdir(os.path.join(dataSetName, dataSets[k]))
            dataExps = np.array([])
            for l in range(len(dirName)):
                if (os.path.isdir(os.path.join(os.path.join(dataSetName, dataSets[k]), dirName[l]))):
                    if (os.path.isfile(os.path.join(os.path.join(os.path.join(dataSetName, dataSets[k]), dirName[l]),
                                                    'titleFile.txt')) and
                            os.path.isfile(
                                os.path.join(os.path.join(os.path.join(dataSetName, dataSets[k]), dirName[l]),
                                             'acqusText.txt')) and
                            os.path.isfile(
                                os.path.join(os.path.join(os.path.join(dataSetName, dataSets[k]), dirName[l]),
                                             'acqu2sText.txt')) and
                            os.path.isfile(
                                os.path.join(os.path.join(os.path.join(dataSetName, dataSets[k]), dirName[l]),
                                             'acqu3sText.txt')) and
                            os.path.isfile(
                                os.path.join(os.path.join(os.path.join(dataSetName, dataSets[k]), dirName[l]),
                                             'procsText.txt')) and
                            os.path.isfile(
                                os.path.join(os.path.join(os.path.join(dataSetName, dataSets[k]), dirName[l]),
                                             'proc2sText.txt')) and
                            os.path.isfile(
                                os.path.join(os.path.join(os.path.join(dataSetName, dataSets[k]), dirName[l]),
                                             'proc3sText.txt')) and
                            os.path.isfile(
                                os.path.join(os.path.join(os.path.join(dataSetName, dataSets[k]), dirName[l]),
                                             'nmrDataSet.dat'))):
                        dataExps = np.append(dataExps, dirName[l])

            dataExps2 = []
            for k in range(len(dataExps)):
                try:
                    dataExps2.append(int(dataExps[k]))
                except:
                    pass

            dataExps2.sort()
            dataExps = list(map(str, dataExps2))
            dataSetExps.append(dataExps)

        self.nmrdat = []
        for k in range(len(dataSetExps)):
            self.nmrdat.append([])
            for l in range(len(dataSetExps[k])):
                fName = os.path.join(os.path.join(os.path.join(dataSetName, dataSets[k]), dataSetExps[k][l]),
                                     'nmrDataSet.dat')
                f = open(fName, 'rb')
                n = pickle.load(f)
                f.close()
                nd2 = nd.NmrData()
                if hasattr(n, 'fid'):
                    nd2.fid = n.fid

                if hasattr(n, 'spc'):
                    nd2.spc = n.spc

                if hasattr(n, 'ppm1'):
                    nd2.ppm1 = n.ppm1

                if hasattr(n, 'ppm2'):
                    nd2.ppm2 = n.ppm2

                if hasattr(n, 'ppm3'):
                    nd2.ppm3 = n.ppm3

                if hasattr(n, 'dim'):
                    nd2.dim = n.dim

                if hasattr(n, 'title'):
                    nd2.title = n.title

                if hasattr(n, 'origDataSet'):
                    nd2.origDataSet = n.origDataSet

                if hasattr(n, 'phCorrMode'):
                    nd2.phCorrMode = n.phCorrMode

                if hasattr(n, 'acq'):
                    a = n.acq
                    aq = nd2.acq
                    if hasattr(a, 'acqusText'):
                        aq.acqusText = a.acqusText

                    if hasattr(a, 'acqu2sText'):
                        aq.acqu2sText = a.acqu2sText

                    if hasattr(a, 'acqu3sText'):
                        aq.acqu3sText = a.acqu3sText

                    if hasattr(a, 'byteOrder'):
                        aq.byteOrder = a.byteOrder

                    if hasattr(a, 'sw'):
                        aq.sw = a.sw

                    if hasattr(a, 'sw_h'):
                        aq.sw_h = a.sw_h

                    if hasattr(a, 'sfo1'):
                        aq.sfo1 = a.sfo1

                    if hasattr(a, 'sfo2'):
                        aq.sfo2 = a.sfo2

                    if hasattr(a, 'sfo3'):
                        aq.sfo3 = a.sfo3

                    if hasattr(a, 'bf1'):
                        aq.bf1 = a.bf1

                    if hasattr(a, 'bf2'):
                        aq.bf2 = a.bf2

                    if hasattr(a, 'bf3'):
                        aq.bf3 = a.bf3

                    if hasattr(a, 'o1'):
                        aq.o1 = a.o1

                    if hasattr(a, 'o2'):
                        aq.o2 = a.o2

                    if hasattr(a, 'o3'):
                        aq.o3 = a.o3

                    if hasattr(a, 'nDataPoints'):
                        aq.nDataPoints = a.nDataPoints

                    if hasattr(a, 'aqMode'):
                        aq.aqMode = a.aqMode

                    if hasattr(a, 'decim'):
                        aq.decim = a.decim

                    if hasattr(a, 'dspfvs'):
                        aq.dspfvs = a.dspfvs

                    if hasattr(a, 'groupDelay'):
                        aq.groupDelay = a.groupDelay

                    if hasattr(a, 'digMod'):
                        aq.digMod = a.digMod

                    if hasattr(a, 'transients'):
                        aq.transients = a.transients

                    if hasattr(a, 'steadyStateScans'):
                        aq.steadyStateScans = a.steadyStateScans

                    if hasattr(a, 'relaxationDelay'):
                        aq.relaxationDelay = a.relaxationDelay

                    if hasattr(a, 'spinRate'):
                        aq.spinRate = a.spinRate

                    if hasattr(a, 'nndp'):
                        aq.nndp = a.nndp

                    if hasattr(a, 'pulseProgram'):
                        aq.pulseProgram = a.pulseProgram

                    if hasattr(a, 'pulProgName'):
                        aq.pulProgName = a.pulProgName

                    if hasattr(a, 'instrument'):
                        aq.instrument = a.instrument

                    if hasattr(a, 'dataType'):
                        aq.dataType = a.dataType

                    if hasattr(a, 'solvent'):
                        aq.solvent = a.solvent

                    if hasattr(a, 'probe'):
                        aq.probe = a.probe

                    if hasattr(a, 'title'):
                        aq.title = a.title

                    if hasattr(a, 'origin'):
                        aq.origin = a.origin

                    if hasattr(a, 'owner'):
                        aq.owner = a.owner

                    if hasattr(a, 'metaInfo'):
                        aq.metaInfo = a.metaInfo

                    if hasattr(a, 'aunm'):
                        aq.aunm = a.aunm

                    if hasattr(a, 'temperature'):
                        aq.temperature = a.temperature

                    if hasattr(a, 'cnst'):
                        aq.cnst = a.cnst

                    if hasattr(a, 'delay'):
                        aq.delay = a.delay

                    if hasattr(a, 'pulse'):
                        aq.pulse = a.pulse

                    if hasattr(a, 'pcpd'):
                        aq.pcpd = a.pcpd

                    if hasattr(a, 'powerLevel'):
                        aq.powerLevel = a.powerLevel

                    if hasattr(a, 'powerLevelWatt'):
                        aq.powerLevelWatt = a.powerLevelWatt

                    if hasattr(a, 'powerLevelMax'):
                        aq.powerLevelMax = a.powerLevelMax

                    if hasattr(a, 'shapedPower'):
                        aq.shapedPower = a.shapedPower

                    if hasattr(a, 'shapedPowerWatt'):
                        aq.shapedPowerWatt = a.shapedPowerWatt

                    if hasattr(a, 'spoal'):
                        aq.spoal = a.spoal

                    if hasattr(a, 'spoffs'):
                        aq.spoffs = a.spoffs

                    if hasattr(a, 'cpdProg'):
                        aq.cpdProg = a.cpdProg

                    if hasattr(a, 'gpName'):
                        aq.gpName = a.gpName

                    if hasattr(a, 'vcList'):
                        aq.vcList = a.vcList

                    if hasattr(a, 'vdList'):
                        aq.vdList = a.vdList

                    if hasattr(a, 'vpList'):
                        aq.vpList = a.vpList

                    if hasattr(a, 'vaList'):
                        aq.vaList = a.vaList

                    if hasattr(a, 'vtList'):
                        aq.vtList = a.vtList

                    if hasattr(a, 'nuc1'):
                        aq.nuc1 = a.nuc1

                    if hasattr(a, 'nuc2'):
                        aq.nuc2 = a.nuc2

                    if hasattr(a, 'nuc3'):
                        aq.nuc3 = a.nuc3

                    if hasattr(a, 'nuc4'):
                        aq.nuc4 = a.nuc4

                    if hasattr(a, 'nuc5'):
                        aq.nuc5 = a.nuc5

                    if hasattr(a, 'nuc6'):
                        aq.nuc6 = a.nuc6

                    if hasattr(a, 'nuc7'):
                        aq.nuc7 = a.nuc7

                    if hasattr(a, 'nuc8'):
                        aq.nuc8 = a.nuc8

                    if hasattr(a, 'gpx'):
                        aq.gpx = a.gpx

                    if hasattr(a, 'gpy'):
                        aq.gpy = a.gpy

                    if hasattr(a, 'gpz'):
                        aq.gpz = a.gpz

                    if hasattr(a, 'increments'):
                        aq.increments = a.increments

                    if hasattr(a, 'nusList'):
                        aq.nusList = a.nusList

                    if hasattr(a, 'nusAmount'):
                        aq.nusAmount = a.nusAmount

                    if hasattr(a, 'nusSeed'):
                        aq.nusSeed = a.nusSeed

                    if hasattr(a, 'nusJsp'):
                        aq.nusJsp = a.nusJsp

                    if hasattr(a, 'nusT2'):
                        aq.nusT2 = a.nusT2

                    if hasattr(a, 'nusTD'):
                        aq.nusTD = a.nusTD

                    if hasattr(a, 'overFlow'):
                        aq.overFlow = a.overFlow

                    if hasattr(a, 'pynm'):
                        aq.pynm = a.pynm

                    if hasattr(a, 'spcFrequency'):
                        aq.spcFrequency = a.spcFrequency

                    if hasattr(a, 'spcSFreq'):
                        aq.spcSFreq = a.spcSFreq

                    if hasattr(a, 'spcNucleus'):
                        aq.spcNucleus = a.spcNucleus

                    if hasattr(a, 'spcOffset'):
                        aq.spcOffset = a.spcOffset

                    if hasattr(a, 'acqT0'):
                        aq.acqT0 = a.acqT0

                    if hasattr(a, 'fnMode'):
                        aq.fnMode = a.fnMode

                    if hasattr(a, 'inf'):
                        aq.inf = a.inf

                    nd2.acq = aq

                if hasattr(n, 'proc'):
                    #nd2.proc = n.proc
                    p = n.proc
                    pc = nd2.proc
                    if hasattr(p, 'procsText'):
                        pc.procsText = p.procsText

                    if hasattr(p, 'proc2sText'):
                        pc.proc2sText = p.proc2sText

                    if hasattr(p, 'proc3sText'):
                        pc.proc3sText = p.proc3sText

                    if hasattr(p, 'tilt'):
                        pc.tilt = p.tilt

                    if hasattr(p, 'symj'):
                        pc.symj = p.symj

                    if hasattr(p, 'ph0'):
                        pc.ph0 = p.ph0

                    if hasattr(p, 'ph1'):
                        pc.ph1 = p.ph1

                    if hasattr(p, 'phCorr'):
                        pc.phCorr = p.phCorr

                    if hasattr(p, 'refShift'):
                        pc.refShift = p.refShift

                    if hasattr(p, 'refPoint'):
                        pc.refPoint = p.refPoint

                    if hasattr(p, 'nPoints'):
                        pc.nPoints = p.nPoints

                    if hasattr(p, 'pivot'):
                        pc.pivot = p.pivot

                    if hasattr(p, 'lb'):
                        pc.lb = p.lb

                    if hasattr(p, 'gb'):
                        pc.gb = p.gb

                    if hasattr(p, 'ssb'):
                        pc.ssb = p.ssb

                    if hasattr(p, 'axisNucleus'):
                        pc.axisNucleus = p.axisNucleus

                    if hasattr(p, 'aunmp'):
                        pc.aunmp = p.aunmp

                    if hasattr(p, 'polyOrder'):
                        pc.polyOrder = p.polyOrder

                    if hasattr(p, 'waterSuppression'):
                        pc.waterSuppression = p.waterSuppression

                    if hasattr(p, 'gibbs'):
                        pc.gibbs = p.gibbs

                    if hasattr(p, 'convExtrapolationSize'):
                        pc.convExtrapolationSize = p.convExtrapolationSize

                    if hasattr(p, 'convWindowSize'):
                        pc.convWindowSize = p.convWindowSize

                    if hasattr(p, 'windowType'):
                        pc.windowType = p.windowType

                    if hasattr(p, 'convWindowType'):
                        pc.convWindowType = p.convWindowType

                    if hasattr(p, 'sw_h'):
                        pc.sw_h = p.sw_h

                    if hasattr(p, 'fidOffsetCorrection'):
                        pc.fidOffsetCorrection = p.fidOffsetCorrection

                    nd2.proc = pc

                if hasattr(n, 'disp'):
                    d = n.disp
                    dp = nd2.disp
                    if hasattr(d, 'posCol'):
                        dp.posCol = d.posCol

                    if hasattr(d, 'posColRGB'):
                        dp.posColRGB = d.posColRGB

                    if hasattr(d, 'negCol'):
                        dp.negCol = d.negCol

                    if hasattr(d, 'negColRGB'):
                        dp.negColRGB = d.negColRGB

                    if hasattr(d, 'phRefCol'):
                        dp.phRefCol = d.phRefCol

                    if hasattr(d, 'phRefDS'):
                        dp.phRefDS = d.phRefDS

                    if hasattr(d, 'phRefExp'):
                        dp.phRefExp = d.phRefExp

                    if hasattr(d, 'nLevels'):
                        dp.nLevels = d.nLevels

                    if hasattr(d, 'minLevel'):
                        dp.minLevel = d.minLevel

                    if hasattr(d, 'maxLevel'):
                        dp.maxLevel = d.maxLevel

                    if hasattr(d, 'axisType1'):
                        dp.axisType1 = d.axisType1

                    if hasattr(d, 'axisType2'):
                        dp.axisType2 = d.axisType2

                    if hasattr(d, 'displaySpc'):
                        dp.displaySpc = d.displaySpc

                    if hasattr(d, 'spcOffset'):
                        dp.spcOffset = d.spcOffset

                    if hasattr(d, 'spcScale'):
                        dp.spcScale = d.spcScale

                    if hasattr(d, 'xLabel'):
                        dp.xLabel = d.xLabel

                    if hasattr(d, 'yLabel'):
                        dp.yLabel = d.yLabel

                    if hasattr(d, 'spcLabel'):
                        dp.spcLabel = d.spcLabel

                    nd2.disp = dp

                if hasattr(n, 'fidOffsetCorr'):
                    nd2.fidOffsetCorr = n.fidOffsetCorr

                if hasattr(n, 'dataSetName'):
                    nd2.dataSetName = n.dataSetName

                if hasattr(n, 'dataSetNumber'):
                    nd2.dataSetNumber = n.dataSetNumber

                if hasattr(n, 'title'):
                    nd2.title = n.title

                if hasattr(n, 'pulseProgram'):
                    nd2.pulseProgram = n.pulseProgram

                if hasattr(n, 'windowFunction'):
                    nd2.windowFunction = n.windowFunction

                if hasattr(n, 'refShift'):
                    nd2.refShift = n.refShift

                if hasattr(n, 'refPoint'):
                    nd2.refPoint = n.refPoint

                if hasattr(n, 'refsw'):
                    nd2.refsw = n.refsw

                if hasattr(n, 'refTmspRange'):
                    nd2.refTmspRange = n.refTmspRange

                if hasattr(n, 'apc'):
                    ab = n.apc
                    ac = nd2.apc
                    if hasattr(ab, 'correctBaseline'):
                        ac.correctBaseline = ab.correctBaseline

                    if hasattr(ab, 'npts'):
                        ac.npts = ab.npts

                    if hasattr(ab, 'nMax'):
                        ac.nMax = ab.nMax

                    if hasattr(ab, 'nOrder'):
                        ac.nOrder = ab.nOrder

                    if hasattr(ab, 'nbins'):
                        ac.nbins = ab.nbins

                    if hasattr(ab, 'mFact0'):
                        ac.mFact0 = ab.mFact0

                    if hasattr(ab, 'mFact1'):
                        ac.mFact1 = ab.mFact1

                    if hasattr(ab, 'startPts'):
                        ac.startPts = ab.startPts

                    if hasattr(ab, 'endPts'):
                        ac.endPts = ab.endPts

                    if hasattr(ab, 'pars'):
                        ac.pars = ab.pars

                    if hasattr(ab, 'rSpc'):
                        ac.rSpc = ab.rSpc

                    if hasattr(ab, 'iSpc'):
                        ac.iSpc = ab.iSpc

                    nd2.apc = ac

                if hasattr(n, 'projectedJres'):
                    nd2.projectedJres = n.projectedJres

                if hasattr(n, 'origJresSet'):
                    nd2.origJresSet = n.origJresSet

                if hasattr(n, 'origJresExp'):
                    nd2.origJresExp = n.origJresExp

                if hasattr(n, 'pjresMode'):
                    nd2.pjresMode = n.pjresMode

                self.nmrdat[k].append(nd2)
                nd2 = []
                fName = os.path.join(os.path.join(os.path.join(dataSetName, dataSets[k]), dataSetExps[k][l]),
                                     'titleFile.txt')
                f = open(fName, 'r')
                self.nmrdat[k][l].title = f.read()
                f.close()
                fName = os.path.join(os.path.join(os.path.join(dataSetName, dataSets[k]), dataSetExps[k][l]),
                                     'acqusText.txt')
                f = open(fName, 'r')
                self.nmrdat[k][l].acq.acqusText = f.read()
                f.close()
                fName = os.path.join(os.path.join(os.path.join(dataSetName, dataSets[k]), dataSetExps[k][l]),
                                     'acqu2sText.txt')
                f = open(fName, 'r')
                self.nmrdat[k][l].acq.acqu2sText = f.read()
                f.close()
                fName = os.path.join(os.path.join(os.path.join(dataSetName, dataSets[k]), dataSetExps[k][l]),
                                     'acqu3sText.txt')
                f = open(fName, 'r')
                self.nmrdat[k][l].acq.acqu3sText = f.read()
                f.close()
                fName = os.path.join(os.path.join(os.path.join(dataSetName, dataSets[k]), dataSetExps[k][l]),
                                     'procsText.txt')
                f = open(fName, 'r')
                self.nmrdat[k][l].proc.procsText = f.read()
                f.close()
                fName = os.path.join(os.path.join(os.path.join(dataSetName, dataSets[k]), dataSetExps[k][l]),
                                     'proc2sText.txt')
                f = open(fName, 'r')
                self.nmrdat[k][l].proc.proc2sText = f.read()
                f.close()
                fName = os.path.join(os.path.join(os.path.join(dataSetName, dataSets[k]), dataSetExps[k][l]),
                                     'proc3sText.txt')
                f = open(fName, 'r')
                self.nmrdat[k][l].proc.proc2sText = f.read()
                f.close()

        # end load

    def noiseFiltering(self):
        val    = self.pp.noiseThreshold*self.pp.stdVal
        for k in range(len(self.nmrdat[self.s])):
            idx =  np.where(self.nmrdat[self.s][k].spc[0].real < val)
            self.deselect2[idx] += np.ones(len(idx))
            if(False):
                idx2 = np.where(self.nmrdat[self.s][k].spc[0].real < 0)
                self.nmrdat[self.s][k].spc[0][idx2] = np.zeros(len(idx2))


        # end noiseFiltering
        
    def noiseFilteringInit(self):
        spcIdx         = np.where((self.nmrdat[self.s][0].ppm1>self.pp.noiseStart) & (self.nmrdat[self.s][0].ppm1<self.pp.noiseEnd))
        self.pp.stdVal = np.std(self.nmrdat[self.s][0].spc[0][spcIdx].real)
        # end noiseFilteringInit
        
    def pjres(self, set = -1, mode = 'skyline'):
        if(set<1):
            set = self.s + 2

        if(len(self.nmrdat)<set):
            self.nmrdat.append([])

        self.nmrdat[set-1] = []
        for k in range(len(self.nmrdat[self.s])):
            nd1 = nd.NmrData()
            self.nmrdat[set-1].append(nd1)
            self.nmrdat[set-1][k].dim = 1
            self.nmrdat[set-1][k].fid = [[]]
            self.nmrdat[set-1][k].acq = self.nmrdat[self.s][k].acq
            self.nmrdat[set-1][k].proc = self.nmrdat[self.s][k].proc
            self.nmrdat[set-1][k].disp = self.nmrdat[self.s][k].disp
            self.nmrdat[set-1][k].ppm1 = np.resize(self.nmrdat[set-1][k].ppm1, (len(self.nmrdat[self.s][k].ppm1)))
            self.nmrdat[set-1][k].ppm1 = np.copy(self.nmrdat[self.s][k].ppm1)
            self.nmrdat[set-1][k].spc = np.resize(self.nmrdat[set-1][k].spc, (1, len(self.nmrdat[self.s][k].spc[0])))
            self.nmrdat[set-1][k].refShift = self.nmrdat[self.s][k].refShift
            self.nmrdat[set-1][k].refPoint = self.nmrdat[self.s][k].refPoint
            self.nmrdat[set-1][k].title = "pJres spectrum\n" + self.nmrdat[self.s][k].title
            self.nmrdat[set-1][k].pjresMode = mode
            self.nmrdat[set-1][k].projectedJres = True
            self.nmrdat[set-1][k].origJresSet = self.s
            self.nmrdat[set-1][k].origJresExp = k
            if(mode == 'skyline'):
                self.nmrdat[set-1][k].spc[0] = np.max(self.nmrdat[self.s][k].spc,0)
            else:
                self.nmrdat[set-1][k].spc[0] = np.sum(self.nmrdat[self.s][k].spc,0)


        # end jresProject

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
        if len(dataExp) > 1:
            for k in range(len(dataExp)):
                self.readSpc(dataPath[0], str(dataExp[k]))


        else:
            for k in range(len(dataPath)):
                self.readSpc(dataPath[k], str(dataExp[0]))


    # end readSpcs

    def readNMRPipeSpc(self, dataSetName, dataSetNumber, procDataName = 'test.dat'):
        self.e            = len(self.nmrdat[self.s])
        nd1               = nd.NmrData()
        nd1.dataSetName   = dataSetName
        nd1.dataSetNumber = dataSetNumber
        nd1.readSpc()
        nd1.readPipe2D(dataSetName + os.sep + dataSetNumber + '.proc', procDataName)
        self.nmrdat[self.s].append(nd1)
        # end readSpc

    def readNMRPipeSpcs(self, dataPath, dataExp, procDataName = 'test.dat'):
        for k in range(len(dataExp)):
            self.readNMRPipeSpc(dataPath, str(dataExp[k]), procDataName)

    # end readSpcs

    def resetDataPreProcessing(self):
        if not self.nmrdat[self.s][0].projectedJres:
            self.ftAll()
            self.baseline1dAll()
            self.autorefAll()
            self.shiftRef()
        else:
            s = self.s
            e = self.e
            self.s = self.nmrdat[s][e].origJresSet
            self.e = self.nmrdat[s][e].origJresExp
            self.pjres(s + 1, self.nmrdat[s][e].pjresMode)
            self.s = s
            self.e = e

        # end dataPreProcessing

    def save(self,dataSetName):
        if(len(dataSetName) == 0):
            return
        
        try:
            os.makedirs(dataSetName)
        except:
            pass
        
        fName = os.path.join(dataSetName,'curPars.dat')
        f     = open(fName,'wb')
        pickle.dump([self.fileFormatVersion, self.s, self.e, self.pp, self.deselect, self.deselect2, self.cmdBuffer, self.cmdIdx, self.script, self.console], f)
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
        
    def segmentalAlignment(self):
        nExp = len(self.nmrdat[self.s])
        nPts = len(self.nmrdat[self.s][0].spc[0])
        spcs = np.zeros((nExp, nPts))
        for k in range(nExp):
            spcs[k] = np.copy(self.nmrdat[self.s][k].spc[0].real)

        segStart = self.nmrdat[self.s][0].ppm2points(self.pp.segStart, 0)
        segEnd   = self.nmrdat[self.s][0].ppm2points(self.pp.segEnd,   0)
        npts     = len(self.nmrdat[self.s][0].spc[0])
        nSpc     = len(self.nmrdat[self.s])
        excludeStart = np.zeros((len(segStart), nSpc))
        excludeEnd = np.zeros((len(segStart), nSpc))
        if self.pp.segAlignRefSpc > 0:
            refSpc = self.nmrdat[self.s][self.pp.segAlignRefSpc - 1].spc[0]
        else:
            refSpcs = np.zeros((nSpc, npts))
            for k in range(nSpc):
                refSpcs[k] = self.nmrdat[self.s][k].spc[0].real

            if self.pp.segAlignRefSpc == 0:
                refSpc = np.mean(refSpcs, 0)
            else:
                refSpc = np.median(refSpcs, 0)

            refSpcs = np.array([[]])

        posShift = np.zeros((nSpc, len(segStart)))
        negShift = np.zeros((nSpc, len(segStart)))
        for k in range(nSpc):
            if k != self.pp.segAlignRefSpc - 1:
                for l in range(len(segStart)):
                    startPts = npts - segEnd[l]
                    endPts   = npts - segStart[l]
                    corrSpc1 = np.copy(refSpc[startPts:endPts].real)
                    corrSpc2 = self.nmrdat[self.s][k].spc[0][startPts:endPts].real
                    maxShift = len(corrSpc1) - 1
                    zeros = np.zeros(len(corrSpc1))
                    corrSpc1 = np.append(np.append(zeros, corrSpc1), zeros)
                    corrSpc2 = np.append(np.append(zeros, corrSpc2), zeros)
                    shifts = np.linspace(-maxShift, maxShift, 2 * maxShift + 1, dtype='int')
                    corrVect = np.zeros(2*maxShift + 1)
                    spcShift = 0
                    for m in shifts:
                        corrVect[m + maxShift] = np.corrcoef(corrSpc1, np.roll(corrSpc2, m))[0][1]

                    maxCorrShifts = shifts[np.where(corrVect == np.max(corrVect))]
                    minShift = np.where(np.abs(maxCorrShifts) == np.min(np.abs(maxCorrShifts)))
                    if np.max(corrVect) > 0.8:
                        spcShift = maxCorrShifts[minShift][0]
                        corrSpc2 = self.nmrdat[self.s][k].spc[0][startPts:endPts]
                        corrSpc2 = np.append(np.append(zeros, corrSpc2), zeros)
                        corrSpc2 = np.roll(corrSpc2, spcShift)[maxShift + 1:2*maxShift + 2]
                        self.nmrdat[self.s][k].spc[0][startPts:endPts] = np.copy(corrSpc2)
                        exSta = 0
                        exEnd = 0
                        if spcShift < 0:
                            negShift[k][l] = 0 - spcShift
                            exEnd = self.nmrdat[self.s][self.pp.segAlignRefSpc - 1].points2ppm(npts - endPts - spcShift, 0)
                            exSta = self.nmrdat[self.s][self.pp.segAlignRefSpc - 1].points2ppm(npts - endPts, 0)
                        elif spcShift > 0:
                            posShift[k][l] = spcShift
                            exEnd = self.nmrdat[self.s][self.pp.segAlignRefSpc - 1].points2ppm(npts - startPts, 0)
                            exSta = self.nmrdat[self.s][self.pp.segAlignRefSpc - 1].points2ppm(npts - startPts - spcShift, 0)

                        excludeStart[l][k] = exSta
                        excludeEnd[l][k] = exEnd





        ps = np.max(posShift, 0)
        ns = np.max(negShift, 0)
        ps2 = np.transpose(posShift)
        ns2 = np.transpose(negShift)
        for k in range(len(posShift[0])):
            l = np.where(ps2[k] == ps[k])[0][0]
            if ps[k] > 0:
                sVal = math.floor(1e4 * excludeStart[k][l]) / 1e4
                eVal = math.floor(1e4 * excludeEnd[k][l]) / 1e4
                if sVal not in self.pp.excludeStart and eVal not in self.pp.excludeEnd:
                    self.pp.excludeStart = np.append(self.pp.excludeStart, sVal)
                    self.pp.excludeEnd   = np.append(self.pp.excludeEnd,   eVal)

            l = np.where(ns2[k] == ns[k])[0][0]
            if ns[k] > 0:
                sVal = math.floor(1e4 * excludeStart[k][l]) / 1e4
                eVal = math.floor(1e4 * excludeEnd[k][l]) / 1e4
                if sVal not in self.pp.excludeStart and eVal not in self.pp.excludeEnd:
                    self.pp.excludeStart = np.append(self.pp.excludeStart, math.floor(1e4 * excludeStart[k][l]) / 1e4)
                    self.pp.excludeEnd = np.append(self.pp.excludeEnd, math.floor(1e4 * excludeEnd[k][l]) / 1e4)



        self.excludeRegion()

    # end segmentalAlignment

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
        
