'''
NMR spectrum processing parameters
'''

import numpy as np
import os
from metabolabpy.nmr import procRegEx
from metabolabpy.nmr import procProcparRegEx
import math

class ProcPars:
    
    def __init__(self):
        self.procsText             = str('')
        self.proc2sText            = str('')
        self.proc3sText            = str('')
        self.tilt                  = False
        self.symj                  = False
        self.ph0                   = np.array([0.0, 0.0, 0.0])
        self.ph1                   = np.array([0.0, 0.0, 0.0])
        self.phCorr                = np.array([1, 1, 0], dtype='int')
        self.refShift              = np.array([0.0, 0.0, 0.0])
        self.refPoint              = np.array([0, 0, 0])
        self.nPoints               = np.array([0, 0, 0])
        self.pivot                 = np.array([0, 0, 0])
        self.lb                    = np.array([0.3, 0.0, 0.0])
        self.gb                    = np.array([0.0, 0.0, 0.0])
        self.ssb                   = np.array([0.0, 0.0, 0.0])
        self.axisNucleus           = np.array(['      ', '      ', '      '], dtype='str')
        self.aunmp                 = str('')
        self.polyOrder             = 4
        self.waterSuppression      = 0
        self.gibbs                 = np.array([True, True, False])
        self.convExtrapolationSize = np.array([32, 0, 0])
        self.convWindowSize        = np.array([32, 0, 0])
        self.windowType            = np.array([1, 4, 0], dtype='int')
        self.convWindowType        = np.array([0, 0, 0], dtype='int')
        self.regEx                 = procRegEx.ProcRegEx()
        self.regExVarian           = procProcparRegEx.ProcProcparRegEx()
        self.sw_h                  = np.array([0.0, 0.0, 0.0])
        self.sf                    = np.array([0.0, 0.0, 0.0])
        self.fidOffsetCorrection   = 0
        self.stripStart            = 0
        self.stripEnd              = 0
        self.phaseInversion        = True
        self.multFactor            = [0, 0]
        self.invertMatrix          = [False, False, False]
        self.windowFunctions       = {
            0             : "None",
            1             : "Exponential",
            2             : "Gaussian",
            3             : "Sine",
            4             : "QSine",
            5             : "SEM",
            "None"        : 0,
            "Exponential" : 1,
            "Gaussian"    : 2,
            "Sine"        : 3,
            "QSine"       : 4,
            "SEM"         : 5
        }
        self.phaseCorrection       = {
            0           : "None",
            1           : "Manual",
            2           : "Automatic",
            "None"      : 0,
            "Manual"    : 1,
            "Automatic" : 2
        }
        self.winType = {
            0          : "Gaussian",
            1          : "Sine",
            "Gaussian" : 0,
            "Sine"     : 1
        }
        self.gibbsP = {
            0     : False,
            1     : True,
            False : 0,
            True  : 1
        }
        self.waterSupp = {
            0         : "None",
            1         : "Conv",
            2         : "Poly",
            3         : "Wavewat",
            "None"    : 0,
            "Conv"    : 1,
            "Poly"    : 2,
            "Wavewat" : 3
        }
        self.dataType = ''
        # end __init__

    def __str__(self): # pragma: no cover
        return self.procsText
        # end __str__
        
    def parseRegEx(self):
        self.ph0[0]         = float(self.regEx.ph0.findall(self.procsText)[0][0])
        self.ph1[0]         = float(self.regEx.ph1.findall(self.procsText)[0][0])
        self.nPoints[0]     = int(self.regEx.nPoints.findall(self.procsText)[0])
        self.lb[0]          = float(self.regEx.lb.findall(self.procsText)[0][0])
        self.gb[0]          = float(self.regEx.gb.findall(self.procsText)[0][0])
        self.ssb[0]         = float(self.regEx.ssb.findall(self.procsText)[0][0])
        self.sf[0]          = float(self.regEx.sf.findall(self.procsText)[0][0])
        #self.stripStart     = int(self.regEx.stsr.findall(self.procsText)[0][0])
        #self.stripEnd       = int(self.regEx.stsi.findall(self.procsText)[0][0])
        if self.ssb[0] == 2.0:
            self.ssb[0] = 90.0

        if self.procsText.find("$AXNUC=") > -1:
            self.axisNucleus[0] = self.regEx.axisNucleus.findall(self.procsText)[0]
            self.axisNucleus[0] = self.axisNucleus[0].replace('<', '')
            self.axisNucleus[0] = self.axisNucleus[0].replace('>', '')

        if self.proc2sText.find("$AXNUC=") > -1:
            self.axisNucleus[1] = self.regEx.axisNucleus.findall(self.proc2sText)[0]
            self.axisNucleus[1] = self.axisNucleus[1].replace('<', '')
            self.axisNucleus[1] = self.axisNucleus[1].replace('>', '')

        self.windowType[0]  = int(self.regEx.wdw.findall(self.procsText)[0])
        self.aunmp          = self.regEx.aunmp.findall(self.procsText)[0]
        try:
            self.sf[1] = float(self.regEx.sf.findall(self.proc2sText)[0][0])
        except:
            pass

        try:
            self.ph0[1]         = float(self.regEx.ph0.findall(self.proc2sText)[0][0])
            self.ph1[1]         = float(self.regEx.ph1.findall(self.proc2sText)[0][0])
            self.nPoints[1]     = int(self.regEx.nPoints.findall(self.proc2sText)[0])
            self.lb[1]          = float(self.regEx.lb.findall(self.proc2sText)[0][0])
            self.gb[1]          = float(self.regEx.gb.findall(self.proc2sText)[0][0])
            self.ssb[1]         = float(self.regEx.ssb.findall(self.proc2sText)[0][0])
            if self.ssb[1] == 2.0:
                self.ssb[1] = 90.0

            self.wdw[1]         = float(self.regEx.wdw.findall(self.proc2sText)[0][0])
            #self.axisNucleus[1] = self.regEx.axisNucleus.findall(self.proc2sText)[0]
            self.windowType[1]  = int(self.regEx.wdw.findall(self.proc2sText)[0])

        except:
            a = 2+3
            #print("1D data")
        
        # end parseRegEx

    def parseRegExVarian(self):
        dd        = self.regExVarian.td.search(self.procsText)
        if hasattr(dd, 'span'):
            dd        = self.procsText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            td        = int(dd)
            self.nPoints[0] = pow(2, math.ceil(math.log(td)/math.log(2)))

        ni = 1
        dd        = self.regExVarian.ni.search(self.procsText)
        if hasattr(dd, 'span'):
            dd        = self.procsText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            ni = int(dd)

        ni = max(ni,1)
        ni2 = 1
        dd        = self.regExVarian.ni2.search(self.procsText)
        if hasattr(dd, 'span'):
            dd        = self.procsText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            ni2 = int(dd)

        ni2 = max(ni2,1)
        phase = np.array([1], dtype='int')
        dd        = self.regExVarian.phase.search(self.procsText)
        if hasattr(dd, 'span'):
            dd         = self.procsText[dd.span()[0] + 1:]
            dd         = dd[dd.find('\n') + 1:]
            dd         = dd[dd.find(' '):dd.find('\n')]
            phase = np.array(dd.split(), dtype='int')

        phase2 = np.array([1], dtype='int')
        dd        = self.regExVarian.phase2.search(self.procsText)
        if hasattr(dd, 'span'):
            dd         = self.procsText[dd.span()[0] + 1:]
            dd         = dd[dd.find('\n') + 1:]
            dd         = dd[dd.find(' '):dd.find('\n')]
            phase2 = np.array(dd.split(), dtype='int')

        self.nPoints[1] = pow(2, math.ceil(math.log(ni*ni2*len(phase)*len(phase2))/math.log(2)))
        dd        = self.regExVarian.nuc1.search(self.procsText)
        if hasattr(dd, 'span'):
            dd        = self.procsText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            dd        = dd.replace(" ", "")
            dd        = dd.replace("\"", "")
            self.axisNucleus[0] = dd

        dd        = self.regExVarian.nuc2.search(self.procsText)
        if hasattr(dd, 'span'):
            dd        = self.procsText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            dd        = dd.replace(" ", "")
            dd        = dd.replace("\"", "")
            self.axisNucleus[1] = dd

        dd        = self.regExVarian.nuc3.search(self.procsText)
        if hasattr(dd, 'span'):
            dd        = self.procsText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            dd        = dd.replace(" ", "")
            dd        = dd.replace("\"", "")
            self.axisNucleus[2] = dd

        # end parseRegExVarian

    def read(self, spcDir):
        procsName    = spcDir + os.sep + 'pdata' + os.sep + '1' + os.sep + 'procs'
        proc2sName   = spcDir + os.sep + 'pdata' + os.sep + '1' + os.sep + 'proc2s'
        proc3sName   = spcDir + os.sep + 'pdata' + os.sep + '1' + os.sep + 'proc3s'
        procparName  = spcDir + os.sep + 'procpar'
        procparName2 = spcDir + os.sep + 'PROCPAR'
        if(os.path.isfile(procsName)):
            try:
                f = open(procsName, "r")
                self.procsText = f.read()
                f.close()

            except:
                f = open(procsName, "r", encoding='latin-1')
                self.procsText = f.read()
                f.close()

            self.dataType = 'Bruker'

        if(os.path.isfile(proc2sName)):
            try:
                f = open(proc2sName,"r")
                self.proc2sText = f.read()
                f.close()

            except:
                f = open(proc2sName,"r", encoding='latin-1')
                self.proc2sText = f.read()
                f.close()

            self.dataType = 'Bruker'

        if(os.path.isfile(proc3sName)):
            try:
                f = open(proc3sName,"r")
                self.proc3sText = f.read()
                f.close()

            except:
                f = open(proc3sName,"r", encoding='latin-1')
                self.proc3sText = f.read()
                f.close()

            self.dataType = 'Bruker'

        if (os.path.isfile(procparName)):
            f = open(procparName, "r")
            self.procsText = f.read()
            f.close()
            self.dataType = 'Varian'

        if (os.path.isfile(procparName2)):
            f = open(procparName2, "r")
            self.procsText = f.read()
            f.close()
            self.dataType = 'Varian'

        if self.dataType == 'Bruker':
            self.parseRegEx()
        else:
            self.parseRegExVarian()

        # end read
    
