import numpy as np # pragma: no cover
import os # pragma: no cover
import pickle # pragma: no cover

from metabolabpy.nmr import nmrData as nd # pragma: no cover
from metabolabpy.nmr import nmrPreProc as npp # pragma: no cover
import math # pragma: no cover
from openpyxl import Workbook # pragma: no cover
from string import ascii_uppercase # pragma: no cover
import itertools # pragma: no cover
import webbrowser # pragma: no cover
import matplotlib # pragma: no cover
import matplotlib.pyplot as pl  # pragma: no cover

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
        self.keepZoom          = True
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

    def compressBuckets(self):
        if (len(self.pp.compressStart) != len(self.pp.compressEnd)):
            return

        for k in range(len(self.pp.compressStart)):
            idx = np.where((self.nmrdat[self.s][0].ppm1 > self.pp.compressStart[k]) & (
                        self.nmrdat[self.s][0].ppm1 < self.pp.compressEnd[k]))
            self.deselect[idx] = np.ones(len(idx))
            self.deselect[int(np.round(np.mean(idx)))] = 0
            #print(idx)
            for l in range(len(self.nmrdat[self.s])):
                val = np.sum(self.nmrdat[self.s][l].spc[0][idx])
                self.nmrdat[self.s][l].spc[0][idx] = np.zeros(len(idx))
                self.nmrdat[self.s][l].spc[0][int(np.round(np.mean(idx)))] = val


        # end compressBuckets

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
        
        for k in range(len(self.nmrdat[self.s])):
            self.nmrdat[self.s][k].spc[0][idx]  = np.zeros(len(idx))
            self.nmrdat[self.s][k].spc[0][idx2] = np.zeros(len(idx2))

        if(self.pp.flagBucketSpectra == True):
            self.bucketSpectra()
        
        if(self.pp.flagCompressBuckets == True):
            self.compressBuckets()
            
        if(self.pp.flagScaleSpectra == True):
            self.scaleSpectra()

        spc = np.zeros(len(self.nmrdat[self.s][0].spc[0]))
        nSpc = len(self.nmrdat[self.s])
        for k in range(nSpc):
            spc += self.nmrdat[self.s][k].spc[0].real

        idx = np.where(spc == 0)
        mVal = 0
        for k in range(len(self.nmrdat[self.s])):
            mVal = min(mVal,np.min(self.nmrdat[self.s][k].spc[0].real))

        if mVal < 0 and self.pp.avoidNegativeValues:
            for k in range(len(self.nmrdat[self.s])):
                self.nmrdat[self.s][k].spc[0] -= mVal

        for k in range(len(self.nmrdat[self.s])):
            self.nmrdat[self.s][k].spc[0][idx] = np.zeros(len(idx))

        if(self.pp.flagVarianceStabilisation == True):
            self.varianceStabilisation()
            
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
        if self.pp.exportMethod == 0:
            print("export Excel format")
            if os.path.isdir(self.pp.exportExcelPath) is False:
                os.makedirs(self.pp.exportExcelPath)

            fName = os.path.join(self.pp.exportExcelPath, self.pp.exportExcel)
            wb = Workbook()
            wsNMRdata = wb.active
            wsNMRdata.title = "NMR data"
            wsMetaData = wb.create_sheet("Metadata")
            spc = np.zeros(len(self.nmrdat[self.s][0].spc[0]))
            npts = len(self.nmrdat[self.s][0].spc[0])
            nSpc = len(self.nmrdat[self.s])
            for k in range(nSpc):
                spc += self.nmrdat[self.s][k].spc[0].real

            deselect = np.zeros(npts)
            idx = np.where(spc == 0)
            deselect[idx] = np.ones(len(idx))
            select = np.where(deselect == 0)
            if (self.pp.exportSamplesInRowsCols == 0):  # samples in rows
                colString = []
                for s in itertools.islice(self.iter_all_strings(), len(select[0]) + 2):
                    colString.append(s)

                wsNMRdata['A1'] = 'Name'
                wsNMRdata['B1'] = 'Class / ppm -->'
                for k in range(len(select[0])):
                    wsNMRdata[colString[k+2] + '1'] = str(self.nmrdat[self.s][0].ppm1[select[0][k]])

                for k in range(len(self.nmrdat[self.s])):
                    dse = os.path.split(self.nmrdat[self.s][k].origDataSet)
                    ds = os.path.split(dse[0])
                    wsNMRdata['A' + str(k+2)] = ds[1] + " " + dse[1]
                    wsNMRdata['B' + str(k+2)] = str(self.pp.classSelect[k])
                    for l in range(len(select[0])):
                        wsNMRdata[colString[l+2] + str(k+2)] = str(self.nmrdat[self.s][k].spc[0][select[0][l]].real)



            else:  # samples in cols
                wsNMRdata['A1'] = 'Name'
                wsNMRdata['A2'] = 'Class / ppm -v'
                colString = []
                for s in itertools.islice(self.iter_all_strings(), nSpc + 1):
                    colString.append(s)

                for k in range(len(self.nmrdat[self.s])):
                    dse = os.path.split(self.nmrdat[self.s][k].origDataSet)
                    ds = os.path.split(dse[0])
                    wsNMRdata[colString[k+1] + '1'] = ds[1] + " " + dse[1]
                    wsNMRdata[colString[k+1] + '2'] = str(self.pp.classSelect[k])

                for l in range(len(select[0])):
                    wsNMRdata['A' + str(l+3)] = str(self.nmrdat[self.s][0].ppm1[select[0][l]])
                    for k in range(len(self.nmrdat[self.s])):
                        wsNMRdata[colString[k+1] + str(l+3)] = str(self.nmrdat[self.s][k].spc[0][select[0][l]].real)



            wb.save(fName)

        elif self.pp.exportMethod == 1:
            print("export CSV format")
            if os.path.isdir(self.pp.exportPathName) is False:
                os.makedirs(self.pp.exportPathName)

            fName = os.path.join(self.pp.exportPathName, self.pp.exportFileName)
            f = open(fName, 'w')
            if (self.pp.exportDelimiterTab == True):
                delim = '\t'
            else:
                delim = self.pp.exportCharacter

            spc = np.zeros(len(self.nmrdat[self.s][0].spc[0]))
            for k in range(len(self.nmrdat[self.s])):
                spc += self.nmrdat[self.s][k].spc[0].real

            deselect = np.zeros(len(self.nmrdat[self.s][0].spc[0]))
            idx = np.where(spc == 0)
            deselect[idx] = np.ones(len(idx))
            if (self.pp.exportSamplesInRowsCols == 0):  # samples in rows
                f.write("ppm" + delim + " ")
                for k in range(len(self.nmrdat[self.s][0].ppm1)):
                    if (deselect[k] == 0):
                        f.write(delim + str(self.nmrdat[self.s][0].ppm1[k]))


                f.write("\n")
                for k in range(len(self.nmrdat[self.s])):
                    dse = os.path.split(self.nmrdat[self.s][k].origDataSet)
                    ds = os.path.split(dse[0])
                    f.write(ds[1] + " " + dse[1] + delim + self.pp.classSelect[k])
                    for l in range(len(self.nmrdat[self.s][k].spc[0])):
                        if (deselect[l] == 0):
                            f.write(delim + str(self.nmrdat[self.s][k].spc[0][l].real))


                    f.write("\n")


            else:  # samples in cols
                f.write("ppm")
                for k in range(len(self.nmrdat[self.s])):
                    dse = os.path.split(self.nmrdat[self.s][k].origDataSet)
                    ds = os.path.split(dse[0])
                    f.write(delim + ds[1] + " " + dse[1])

                f.write("\n")
                f.write(" ")
                for k in range(len(self.pp.classSelect)):
                    f.write(delim + self.pp.classSelect[k])

                f.write("\n")
                for k in range(len(self.nmrdat[self.s][0].spc[0])):
                    if (deselect[k] == 0):
                        f.write(str(self.nmrdat[self.s][0].ppm1[k]))
                        for l in range(len(self.nmrdat[self.s])):
                            f.write(delim + str(self.nmrdat[self.s][l].spc[0][k].real))

                        f.write("\n")



            f.close()

        elif self.pp.exportMethod == 2:
            print("export MetaboAnalyst")
            if os.path.isdir(self.pp.exportMetaboAnalystPath) is False:
                os.makedirs(self.pp.exportMetaboAnalystPath)

            fName = os.path.join(self.pp.exportMetaboAnalystPath, self.pp.exportMetaboAnalyst)
            f = open(fName, 'w')
            delim = ','
            spc = np.zeros(len(self.nmrdat[self.s][0].spc[0]))
            for k in range(len(self.nmrdat[self.s])):
                spc += self.nmrdat[self.s][k].spc[0].real

            deselect = np.zeros(len(self.nmrdat[self.s][0].spc[0]))
            idx = np.where(spc == 0)
            deselect[idx] = np.ones(len(idx))
            f.write("Sample" + delim + " Class")
            for k in range(len(self.nmrdat[self.s][0].ppm1)):
                if (deselect[k] == 0):
                    f.write(delim + " Bin." + str(self.nmrdat[self.s][0].ppm1[k]))


            f.write("\n")
            for k in range(len(self.nmrdat[self.s])):
                dse = os.path.split(self.nmrdat[self.s][k].origDataSet)
                ds = os.path.split(dse[0])
                f.write(ds[1] + " " + dse[1] + delim + " " + self.pp.classSelect[k])
                for l in range(len(self.nmrdat[self.s][k].spc[0])):
                    if (deselect[l] == 0):
                        f.write(delim + " " + str(self.nmrdat[self.s][k].spc[0][l].real))


                f.write("\n")


            f.close()

        elif self.pp.exportMethod == 3:
            print("export rDolphin")
            if os.path.isdir(self.pp.exportrDolphinPath) is False:
                os.makedirs(self.pp.exportrDolphinPath)

            fName = os.path.join(self.pp.exportrDolphinPath, self.pp.exportrDolphin)
            f = open(fName, 'w')
            delim = ','
            spc = np.zeros(len(self.nmrdat[self.s][0].spc[0]))
            for k in range(len(self.nmrdat[self.s])):
                spc += self.nmrdat[self.s][k].spc[0].real

            deselect = np.zeros(len(self.nmrdat[self.s][0].spc[0]))
            idx = np.where(spc == 0)
            deselect[idx] = np.ones(len(idx))
            f.write(str(self.nmrdat[self.s][0].ppm1[0]))
            for k in range(1, len(self.nmrdat[self.s][0].ppm1)):
                if (deselect[k] == 0):
                    f.write(delim + str(self.nmrdat[self.s][0].ppm1[k]))


            f.write("\n")
            for k in range(len(self.nmrdat[self.s])):
                dse = os.path.split(self.nmrdat[self.s][k].origDataSet)
                ds = os.path.split(dse[0])
                f.write(str(self.nmrdat[self.s][k].spc[0][0].real))
                for l in range(1, len(self.nmrdat[self.s][k].spc[0])):
                    if (deselect[l] == 0):
                        f.write(delim + str(self.nmrdat[self.s][k].spc[0][l].real))


                f.write("\n")

            f.close()
            fName = os.path.join(self.pp.exportrDolphinPath, "Parameters.csv")
            f = open(fName, 'w')
            f.write("Parameter,Value\n")
            f.write("nmr folder path,\n")
            f.write("1D data index,\n")
            f.write("proc_no,\n")
            f.write("spectra dataset path (csv format)," + self.pp.exportrDolphinPath + "/" + self.pp.exportrDolphin + "\n")
            f.write("Metadata path (csv format)," + self.pp.exportrDolphinPath + "/Metadata.csv\n")
            f.write("ROI patters file," + self.pp.exportrDolphinPath + "/ROI_profile.csv\n")
            f.write("Normalization (0=No;1=Eretic; 2=TSP; 3=Creatinine; 4=Spectra Sum; 5=PQN),2\n")
            f.write("Alignment (0=No;1=Glucose; 2=TSP; 3=Formate),2\n")
            f.write("Suppression,12-9.5;6.1-5.6;5.1-4.5\n")
            f.write("Spectrometer Frequency (MHz)," + str(self.nmrdat[self.s][0].acq.sfo1) + "\n")
            f.write("Bucket resolution," + str(self.pp.bucketPPM) + "\n")
            f.write("Biofluid,Urine\n")
            f.write("2D-Path,\n")
            f.write("Specific program parameters,\n")
            f.close()
            fName = os.path.join(self.pp.exportrDolphinPath, "Metadata.csv")
            f = open(fName, 'w')
            f.write("Sample,Individual,Sample Type\n")
            for k in range(len(self.nmrdat[self.s])):
                f.write(os.path.split(self.nmrdat[self.s][k].origDataSet)[1] + "," + str(k + 1) + ",1\n")

            f.close()

        elif self.pp.exportMethod == 4:
            print("export Batman")
            if os.path.isdir(self.pp.exportBatmanPath) is False:
                os.makedirs(self.pp.exportBatmanPath)

            fName = os.path.join(self.pp.exportBatmanPath, self.pp.exportBatman)
            f = open(fName, 'w')
            delim = '\t'
            spc = np.zeros(len(self.nmrdat[self.s][0].spc[0]))
            for k in range(len(self.nmrdat[self.s])):
                spc += self.nmrdat[self.s][k].spc[0].real

            deselect = np.zeros(len(self.nmrdat[self.s][0].spc[0]))
            idx = np.where(spc == 0)
            deselect[idx] = np.ones(len(idx))
            f.write("ppm")
            for k in range(len(self.nmrdat[self.s])):
                f.write(delim + "X" + str(k + 1))

            f.write("\n")
            for k in range(len(self.nmrdat[self.s][0].spc[0])):
                if (deselect[k] == 0):
                    f.write(str(self.nmrdat[self.s][0].ppm1[k]))
                    for l in range(len(self.nmrdat[self.s])):
                        f.write(delim + str(self.nmrdat[self.s][l].spc[0][k].real))

                    f.write("\n")


            f.close()

        else:
            print("export Bruker")
            if os.path.isdir(self.pp.exportBrukerPath + os.sep + self.pp.exportBruker) is False:
                os.makedirs(self.pp.exportBrukerPath + os.sep + self.pp.exportBruker)

            mmax = 0
            for k in range(len(self.nmrdat[self.s])):
                mmax = max(np.max(self.nmrdat[0][0].spc[0].real), mmax)

            scaleFactor = 2 * mmax / 2147483647
            for k in range(len(self.nmrdat[self.s])):
                self.nmrdat[self.s][k].exportBruker1d(self.pp.exportBrukerPath + os.sep + self.pp.exportBruker, str(k+1), scaleFactor)


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

    def help(self):
        fName = os.path.join(os.path.dirname(__file__), "web", "introduction", "index.html")
        url = "file://" + fName
        webbrowser.open(url, new=2)
        # end help

    def iter_all_strings(self):
        for size in itertools.count(1):
            for s in itertools.product(ascii_uppercase, repeat=size):
                yield "".join(s)

        # end iter_all_strings

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
        for k in self.pp.__dict__.keys():
            if hasattr(c, k):
                exec('self.pp.' + k + '=c.' + k)


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
                for kk in nd2.__dict__.keys():
                    if kk != 'acq' and kk != 'proc' and kk != 'display' and kk != 'apc':
                        if hasattr(n, kk):
                            exec('nd2.' + kk + '=n.' + kk)

                    elif kk == 'acq':
                        if hasattr(n, kk):
                            a = n.acq
                            aq = nd2.acq
                            for kkk in aq.__dict__.keys():
                                if kkk != 'regEx':
                                    if hasattr(a, kkk):
                                        exec('aq.' + kkk + '=a.' + kkk)

                                else:
                                    r  = a.regEx
                                    re = aq.regEx
                                    for kkkk in re.__dict__.keys():
                                        if hasattr(r, kkkk):
                                            exec('re.' + kkkk + '=r.' + kkkk)

                                    aq.regEx = re


                            nd2.acq = aq


                    elif kk == 'proc':
                        if hasattr(n, kk):
                            p = n.proc
                            pc = nd2.proc
                            for kkk in pc.__dict__.keys():
                                if kkk != 'regEx':
                                    if hasattr(p, kkk):
                                        exec('pc.' + kkk + '=p.' + kkk)
                                else:
                                    r = p.regEx
                                    re = pc.regEx
                                    for kkkk in re.__dict__.keys():
                                        if hasattr(r, kkkk):
                                            exec('re.' + kkkk + '=r.' + kkkk)

                                    aq.regEx = re


                        nd2.proc = pc


                    elif kk == 'display':
                        if hasattr(n, kk):
                            d = n.display
                            dp = nd2.display
                            for kkk in dp.__dict__.keys():
                                if hasattr(d, kkk):
                                    exec('dp.' + kkk+ '=d.' + kkk)


                            nd2.display = dp


                    elif kk == 'apc':
                        if hasattr(n, kk):
                            ab = n.apc
                            ac = nd2.apc
                            for kkk in ac.__dict__.keys():
                                if hasattr(ab, kkk):
                                    exec('ac.' + kkk + '=ab.' + kkk)


                            nd2.apc = ac



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
            self.nmrdat[set-1][k].display = self.nmrdat[self.s][k].display
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

    def plotSpc(self):
        pl.plot()
        ax = pl.gca()
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ax.clear()
        if (len(self.nmrdat[self.s]) == 0):
            return

        if (len(self.nmrdat[self.s][self.e].spc) == 0):
            return

        d = self.nmrdat[self.s][self.e].display
        if (d.posCol == "RGB"):
            posCol = d.posColRGB
        else:
            posCol = d.posCol

        if (d.negCol == "RGB"):
            negCol = d.negColRGB
        else:
            negCol = d.negCol

        posCol = matplotlib.colors.to_hex(posCol)
        negCol = matplotlib.colors.to_hex(negCol)
        xlabel = d.xLabel + " [" + d.axisType1 + "]"
        ylabel = d.yLabel + " [" + d.axisType2 + "]"
        if (self.nmrdat[self.s][self.e].dim == 1):
            for k in range(len(self.nmrdat[self.s])):
                if ((k != self.e) and (self.nmrdat[self.s][k].display.displaySpc == True)):
                    d = self.nmrdat[self.s][k].display
                    if (d.posCol == "RGB"):
                        posCol = d.posColRGB
                    else:
                        posCol = d.posCol

                    if (d.negCol == "RGB"):
                        negCol = d.negColRGB
                    else:
                        negCol = d.negCol

                    posCol = matplotlib.colors.to_hex(posCol)
                    negCol = matplotlib.colors.to_hex(negCol)
                    pl.plot(self.nmrdat[self.s][k].ppm1, self.nmrdat[self.s][k].spc[0].real, color=posCol)

            d = self.nmrdat[self.s][self.e].display
            if (d.posCol == "RGB"):
                posCol = d.posColRGB
            else:
                posCol = d.posCol

            if (d.negCol == "RGB"):
                negCol = d.negColRGB
            else:
                negCol = d.negCol

            posCol = matplotlib.colors.to_hex(posCol)
            negCol = matplotlib.colors.to_hex(negCol)
            xlabel = d.xLabel + " [" + d.axisType1 + "]"
            ylabel = d.yLabel + " [" + d.axisType2 + "]"
            pl.plot(self.nmrdat[self.s][self.e].ppm1, self.nmrdat[self.s][self.e].spc[0].real, color=posCol)
            ax = pl.gca()
            ax.set_xlabel(xlabel)
            ax.autoscale()
            ax.invert_xaxis()
            if (self.keepZoom == True):
                if xlim[0] != -0.05 and xlim[1] != 1.05 and ylim[0] != -0.05 and ylim[1] != 1.05:
                    ax.set_xlim(xlim)
                    ax.set_ylim(ylim)


        else:
            mm = np.max(np.abs(self.nmrdat[self.s][self.e].spc.real))
            posLev = np.linspace(d.minLevel * mm, d.maxLevel * mm, d.nLevels)
            negLev = np.linspace(-d.maxLevel * mm, -d.minLevel * mm, d.nLevels)
            pl.contour(self.nmrdat[self.s][self.e].ppm1,
                                                 self.nmrdat[self.s][self.e].ppm2,
                                                 self.nmrdat[self.s][self.e].spc.real, posLev, colors=posCol,
                                                 linestyles='solid', antialiased=True)
            pl.contour(self.nmrdat[self.s][self.e].ppm1,
                                                 self.nmrdat[self.s][self.e].ppm2,
                                                 self.nmrdat[self.s][self.e].spc.real, negLev, colors=negCol,
                                                 linestyles='solid', antialiased=True)
            ax = pl.gca()
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.autoscale()
            ax.invert_xaxis()
            ax.invert_yaxis()

        # end plotSpc

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
        if len(dataExp) > 1:
            for k in range(len(dataExp)):
                self.readNMRPipeSpc(dataPath[0], str(dataExp[k]), procDataName)


        else:
            for k in range(len(dataPath)):
                self.readNMRPipeSpc(dataPath[k], str(dataExp[0]), procDataName)


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
        
    def scaleSpectra(self):
        nSpc     = len(self.nmrdat[self.s])
        npts     = len(self.nmrdat[self.s][0].spc[0])
        if self.pp.scaleSpectraRefSpc > 0:
            refSpc = self.nmrdat[self.s][self.pp.scaleSpectraRefSpc - 1].spc[0].real
        else:
            refSpcs = np.zeros((nSpc, npts))
            for k in range(nSpc):
                refSpcs[k] = self.nmrdat[self.s][k].spc[0].real

            if self.pp.segAlignRefSpc == 0:
                refSpc = np.mean(refSpcs, 0)
            else:
                refSpc = np.median(refSpcs, 0)

            refSpcs = np.array([[]])

        scale = np.ones(nSpc)
        if self.pp.scalePQN is True:
            for k in range(nSpc):
                scaleVect = self.nmrdat[self.s][k].spc[0][np.where(refSpc != 0)].real/refSpc[np.where(refSpc != 0)]
                self.nmrdat[self.s][k].spc[0] /= np.median(scaleVect[np.where(scaleVect != 0)])
        else:
            if self.pp.preserveOverallScale is True:
                for k in range(nSpc):
                    scale[k] = np.sum(self.nmrdat[self.s][k].spc[0]).real

            for k in range(nSpc):
                self.nmrdat[self.s][k].spc[0] /= np.sum(self.nmrdat[self.s][k].spc[0]).real
                self.nmrdat[self.s][k].spc[0] *= np.max(scale)


        # end scaleSpectra

    def segmentalAlignment(self):
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
        excludeStart2 = np.array([])
        excludeEnd2   = np.array([])
        for k in range(len(posShift[0])):
            l = np.where(ps2[k] == ps[k])[0][0]
            if ps[k] > 0:
                sVal = math.floor(1e4 * excludeStart[k][l]) / 1e4
                eVal = math.floor(1e4 * excludeEnd[k][l]) / 1e4
                if sVal not in self.pp.excludeStart and eVal not in self.pp.excludeEnd:
                    excludeStart2 = np.append(excludeStart2, math.floor(1e4 * sVal) / 1e4)
                    excludeEnd2   = np.append(excludeEnd2,   math.floor(1e4 * eVal) / 1e4)
                    #self.pp.excludeStart = np.append(self.pp.excludeStart, sVal)
                    #self.pp.excludeEnd   = np.append(self.pp.excludeEnd,   eVal)

            l = np.where(ns2[k] == ns[k])[0][0]
            if ns[k] > 0:
                sVal = math.floor(1e4 * excludeStart[k][l]) / 1e4
                eVal = math.floor(1e4 * excludeEnd[k][l]) / 1e4
                if sVal not in self.pp.excludeStart and eVal not in self.pp.excludeEnd:
                    excludeStart2 = np.append(excludeStart2, math.floor(1e4 * sVal) / 1e4)
                    excludeEnd2   = np.append(excludeEnd2,   math.floor(1e4 * eVal) / 1e4)
                    #self.pp.excludeStart = np.append(self.pp.excludeStart, math.floor(1e4 * sVal) / 1e4)
                    #self.pp.excludeEnd = np.append(self.pp.excludeEnd, math.floor(1e4 * eVal) / 1e4)



        #self.excludeRegion()
        if self.pp.segAlignRefSpc > 0:
            spcIdx = self.pp.segAlignRefSpc - 1
        else:
            spcIdx = 0

        for k in range(len(excludeStart2)):
            idx = np.where((self.nmrdat[self.s][spcIdx].ppm1 > excludeStart2[k]) & (
                        self.nmrdat[self.s][spcIdx].ppm1 < excludeEnd2[k]))
            self.deselect[idx] = np.ones(len(idx))

    # end segmentalAlignment

    def selectPlotAll(self):
        for k in range(len(self.nmrdat[self.s])):
            self.nmrdat[self.s][k].display.displaySpc = True

        self.plotSpc()
        # end selectPlotAll

    def selectPlotClear(self):
        for k in range(len(self.nmrdat[self.s])):
            self.nmrdat[self.s][k].display.displaySpc = False

        self.plotSpc()
        # end selectPlotClear

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

    def varianceStabilisation(self):
        if self.pp.autoScaling:
            self.varStabAutoScale()

        if self.pp.paretoScaling:
            self.varStabParetoScale()

        if self.pp.gLogTransform:
            self.varStabgLogTransform()

    # end varianceStabilisation

    def varStabAutoScale(self):
        npts = len(self.nmrdat[self.s][0].spc[0])
        nSpc = len(self.nmrdat[self.s])
        spcs = np.zeros((nSpc, npts))
        for k in range(nSpc):
            spcs[k] = self.nmrdat[self.s][k].spc[0].real

        spcMean = np.mean(spcs, 0)
        for k in range(nSpc):
            spcs[k] -= spcMean

        spcVar = np.var(spcs, 0)
        for k in range(nSpc):
            idx = np.where(spcVar != 0)
            spcs[k][idx] = spcs[k][idx]/np.sqrt(spcVar[idx])
            self.nmrdat[self.s][k].spc[0] = spcs[k]

    # end varStabAutoScale

    def varStabgLogTransform(self):
        nSpc = len(self.nmrdat[self.s])
        lMin = np.zeros(nSpc)
        lMax = np.zeros(nSpc)
        for k in range(nSpc):
            lMax[k] = np.max(self.nmrdat[self.s][k].spc[0].real)

        mMax = np.max(lMax)
        for k in range(nSpc):
            spc = np.copy(self.nmrdat[self.s][k].spc[0].real)
            spc /= mMax
            spc = spc - self.pp.varY0 + np.sqrt((spc -self.pp.varY0)**2 + self.pp.varLambda)
            idx = np.where(spc <= 0)
            spc[idx] = 1e-100
            self.nmrdat[self.s][k].spc[0] = np.copy(np.log(spc))
            lMin[k] = np.min(self.nmrdat[self.s][k].spc[0].real)

        for k in range(nSpc):
            self.nmrdat[self.s][k].spc[0] -= np.min(lMin)

    # end varStabgLogTransform

    def varStabParetoScale(self):
        npts = len(self.nmrdat[self.s][0].spc[0])
        nSpc = len(self.nmrdat[self.s])
        spcs = np.zeros((nSpc, npts))
        for k in range(nSpc):
            spcs[k] = self.nmrdat[self.s][k].spc[0].real

        spcVar = np.sqrt(np.std(spcs, 0))
        for k in range(nSpc):
            idx = np.where(spcVar != 0)
            spcs[k][idx] = spcs[k][idx]/spcVar[idx]
            self.nmrdat[self.s][k].spc[0] = spcs[k]

    # end varStabParetoScale
