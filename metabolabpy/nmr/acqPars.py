'''
NMR spectrum acquisition parameters
'''

import numpy as np
import os
from metabolabpy.nmr import acqRegEx

class AcqPars:
    
    def __init__(self):
        self.acqusText        = str('')
        self.acqu2sText       = str('')
        self.acqu3sText       = str('')
        self.byteOrder        = -1
        self.sw               = np.array([0, 0, 0],dtype = 'float64')
        self.sw_h             = np.array([0, 0, 0],dtype = 'float64')
        self.sfo1             = 0.0
        self.sfo2             = 0.0
        self.sfo3             = 0.0
        self.bf1              = 0.0
        self.bf2              = 0.0
        self.bf3              = 0.0
        self.o1               = 0.0
        self.o2               = 0.0
        self.o3               = 0.0
        self.nDataPoints      = np.array([0, 0, 0],dtype = 'int')
        self.aqMode           = 0
        self.decim            = 0
        self.dspfvs           = 0
        self.groupDelay       = 0.0
        self.digMod           = 0
        self.transients       = 0
        self.steadyStateScans = 0
        self.relaxationDelay  = 0.0
        self.spinRate         = 0.0
        self.nndp             = 0.0
        self.pulseProgram     = str('')
        self.pulProgName      = str('')
        self.instrument       = str('')
        self.dataType         = 0
        self.solvent          = str('')
        self.probe            = str('')
        self.title            = str('')
        self.origin           = str('')
        self.owner            = str('')
        self.metaInfo         = str('')
        self.aunm             = str('')
        self.temperature      = -273.15
        self.cnst             = np.array([], 'float64')
        self.delay            = np.array([], 'float64')
        self.pulse            = np.array([], 'float64')
        self.pcpd             = np.array([], 'float64')
        self.powerLevel       = np.array([], 'float64')
        self.powerLevelWatt   = np.array([], 'float64')
        self.powerLevelMax    = np.array([], 'float64')
        self.shapedPower      = np.array([], 'float64')
        self.shapedPowerWatt  = np.array([], 'float64')
        self.spoal            = np.array([], 'float64')
        self.spoffs           = np.array([], 'float64')
        self.cpdProg          = np.array([], dtype = 'str')
        self.gpName           = np.array([], dtype = 'str')
        self.vcList           = str('')
        self.vdList           = str('')
        self.vpList           = str('')
        self.vaList           = str('')
        self.vtList           = str('')
        self.nuc1             = str('')
        self.nuc2             = str('')
        self.nuc3             = str('')
        self.nuc4             = str('')
        self.nuc5             = str('')
        self.nuc6             = str('')
        self.nuc7             = str('')
        self.nuc8             = str('')
        self.gpx              = np.array([], dtype = 'float64')
        self.gpy              = np.array([], dtype = 'float64')
        self.gpz              = np.array([], dtype = 'float64')
        self.increments       = np.array([], dtype = 'float64')
        self.nusList          = str('')
        self.nusAmount        = 0.0
        self.nusSeed          = 0
        self.nusJsp           = 0
        self.nusT2            = 0.0
        self.nusTD            = 0
        self.overFlow         = 0
        self.pynm             = str('')
        self.spcFrequency     = np.array([0, 0], dtype = 'float64')
        self.spcSFreq         = np.array([0, 0], dtype = 'float64')
        self.spcNucleus       = np.array(['     ', '     '], dtype = 'str')
        self.spcOffset        = np.array([0, 0], dtype = 'float64')
        self.acqT0            = np.array([0, 0], dtype = 'int')
        self.fnMode           = 0
        self.inf              = np.array([], dtype = 'float64')
        self.regEx            = acqRegEx.AcqRegEx()
        # end __init__
        
    def __str__(self):
        return self.acqusText
        # end __str__
        
    def parseRegEx(self):
        self.sfo1             = float(self.regEx.sfo1.findall(self.acqusText)[0])
        self.sfo2             = float(self.regEx.sfo2.findall(self.acqusText)[0])
        self.sfo3             = float(self.regEx.sfo3.findall(self.acqusText)[0])
        self.sfo4             = float(self.regEx.sfo4.findall(self.acqusText)[0])
        self.sfo5             = float(self.regEx.sfo5.findall(self.acqusText)[0])
        self.sfo6             = float(self.regEx.sfo6.findall(self.acqusText)[0])
        self.sfo7             = float(self.regEx.sfo7.findall(self.acqusText)[0])
        self.sfo8             = float(self.regEx.sfo8.findall(self.acqusText)[0])
        self.bf1              = float(self.regEx.bf1.findall(self.acqusText)[0])
        self.bf2              = float(self.regEx.bf2.findall(self.acqusText)[0])
        self.bf3              = float(self.regEx.bf3.findall(self.acqusText)[0])
        self.bf4              = float(self.regEx.bf4.findall(self.acqusText)[0])
        self.bf5              = float(self.regEx.bf5.findall(self.acqusText)[0])
        self.bf6              = float(self.regEx.bf6.findall(self.acqusText)[0])
        self.bf7              = float(self.regEx.bf7.findall(self.acqusText)[0])
        self.bf8              = float(self.regEx.bf8.findall(self.acqusText)[0])
        self.o1               = float(self.regEx.o1.findall(self.acqusText)[0])
        self.o2               = float(self.regEx.o2.findall(self.acqusText)[0])
        self.o3               = float(self.regEx.o3.findall(self.acqusText)[0])
        self.o4               = float(self.regEx.o4.findall(self.acqusText)[0])
        try:
            self.o5               = float(self.regEx.o5.findall(self.acqusText)[0])
            self.o6               = float(self.regEx.o6.findall(self.acqusText)[0])
            self.o7               = float(self.regEx.o7.findall(self.acqusText)[0])
            self.o8               = float(self.regEx.o8.findall(self.acqusText)[0])
        except:
            pass
        
        self.sw[0]            = float(self.regEx.sw.findall(self.acqusText)[0])
        self.sw_h[0]          = float(self.regEx.sw_h.findall(self.acqusText)[0])
        self.nDataPoints[0]   = int(self.regEx.td.findall(self.acqusText)[0])
        self.decim            = int(self.regEx.decim.findall(self.acqusText)[0])
        self.dspfvs           = int(self.regEx.dspfvs.findall(self.acqusText)[0])
        self.groupDelay       = float(self.regEx.grpdly.findall(self.acqusText)[0][0])
        self.byteOrder        = int(self.regEx.byteOrder.findall(self.acqusText)[0])
        self.aqMode           = int(self.regEx.aqMode.findall(self.acqusText)[0])
        self.digMod           = int(self.regEx.digMod.findall(self.acqusText)[0])
        self.transients       = int(self.regEx.transients.findall(self.acqusText)[0])
        self.steadyStateScans = int(self.regEx.steadyStateScans.findall(self.acqusText)[0])
        self.relaxationDelay  = float(self.regEx.relaxationDelay.findall(self.acqusText)[0]) if (len(self.regEx.relaxationDelay.findall(self.acqusText))>0) else 0.0
        self.spinRate         = int(self.regEx.spinRate.findall(self.acqusText)[0])
        self.pulProgName      = self.regEx.pulProg.findall(self.acqusText)[0]
        self.aunm             = self.regEx.aunm.findall(self.acqusText)[0]
        self.instrument       = self.regEx.instrument.findall(self.acqusText)[0]
        self.dataType         = int(self.regEx.dataType.findall(self.acqusText)[0])
        self.solvent          = self.regEx.solvent.findall(self.acqusText)[0]
        self.probe            = self.regEx.probe.findall(self.acqusText)[0]
        self.probe            = self.probe.replace('<','')
        self.probe            = self.probe.replace('>','')
        self.title            = self.regEx.title.findall(self.acqusText)[0]
        self.origin           = self.regEx.origin.findall(self.acqusText)[0]
        self.owner            = self.regEx.owner.findall(self.acqusText)[0]
        self.metaInfo         = self.regEx.metaInfo.findall(self.acqusText)[0]
        self.temperature      = float(self.regEx.temperature.findall(self.acqusText)[0])
        dd                    = self.regEx.cnst.search(self.acqusText)
        dd                    = self.acqusText[dd.span()[0]:]
        dd                    = dd[dd.find('\n')+1:]
        dd                    = dd[:dd.find('##$')]
        self.cnst             = np.array(dd.split(), dtype = 'float64')
        dd                    = self.regEx.delay.search(self.acqusText)
        dd                    = self.acqusText[dd.span()[0]:]
        dd                    = dd[dd.find('\n')+1:]
        dd                    = dd[:dd.find('##$')]
        self.delay            = np.array(dd.split(), dtype = 'float64')
        try:
            dd                    = self.regEx.cpdProg.search(self.acqusText)
            dd                    = self.acqusText[dd.span()[0]:]
            dd                    = dd[dd.find('\n')+1:]
            dd                    = dd[:dd.find('##$')]
            self.cpdProg          = np.array(dd.split(), dtype = 'str')
        except:
            pass
            
        try:
            dd                    = self.regEx.gpName.search(self.acqusText)
            dd                    = self.acqusText[dd.span()[0]:]
            dd                    = dd[dd.find('\n')+1:]
            dd                    = dd[:dd.find('##$')]
            self.gpName           = np.array(dd.split(), dtype = 'str')
        except:
            pass
            
        dd                    = self.regEx.gpx.search(self.acqusText)
        dd                    = self.acqusText[dd.span()[0]:]
        dd                    = dd[dd.find('\n')+1:]
        dd                    = dd[:dd.find('##$')]
        self.gpx              = np.array(dd.split(), dtype = 'float64')
        dd                    = self.regEx.gpy.search(self.acqusText)
        dd                    = self.acqusText[dd.span()[0]:]
        dd                    = dd[dd.find('\n')+1:]
        dd                    = dd[:dd.find('##$')]
        self.gpy              = np.array(dd.split(), dtype = 'float64')
        dd                    = self.regEx.gpz.search(self.acqusText)
        dd                    = self.acqusText[dd.span()[0]:]
        dd                    = dd[dd.find('\n')+1:]
        dd                    = dd[:dd.find('##$')]
        self.gpz              = np.array(dd.split(), dtype = 'float64')
        dd                    = self.regEx.increments.search(self.acqusText)
        dd                    = self.acqusText[dd.span()[0]:]
        dd                    = dd[dd.find('\n')+1:]
        dd                    = dd[:dd.find('##$')]
        self.increments       = np.array(dd.split(), dtype = 'float64')
        dd                    = self.regEx.pulse.search(self.acqusText)
        dd                    = self.acqusText[dd.span()[0]:]
        dd                    = dd[dd.find('\n')+1:]
        dd                    = dd[:dd.find('##$')]
        self.pulse            = np.array(dd.split(), dtype = 'float64')
        dd                    = self.regEx.pcpd.search(self.acqusText)
        dd                    = self.acqusText[dd.span()[0]:]
        dd                    = dd[dd.find('\n')+1:]
        dd                    = dd[:dd.find('##$')]
        self.pcpd            = np.array(dd.split(), dtype = 'float64')
        dd                    = self.regEx.powerLevel.search(self.acqusText)
        dd                    = self.acqusText[dd.span()[0]:]
        dd                    = dd[dd.find('\n')+1:]
        dd                    = dd[:dd.find('##$')]
        self.powerLevel       = np.array(dd.split(), dtype = 'float64')
        dd                    = self.regEx.powerLevelWatt.search(self.acqusText)
        dd                    = self.acqusText[dd.span()[0]:]
        dd                    = dd[dd.find('\n')+1:]
        dd                    = dd[:dd.find('##$')]
        self.powerLevelWatt   = np.array(dd.split(), dtype = 'float64')
        dd                    = self.regEx.powerLevelMax.search(self.acqusText)
        dd                    = self.acqusText[dd.span()[0]:]
        dd                    = dd[dd.find('\n')+1:]
        dd                    = dd[:dd.find('##$')]
        self.powerLevelMax    = np.array(dd.split(), dtype = 'float64')
        dd                    = self.regEx.shapedPower.search(self.acqusText)
        dd                    = self.acqusText[dd.span()[0]:]
        dd                    = dd[dd.find('\n')+1:]
        dd                    = dd[:dd.find('##$')]
        self.shapedPower      = np.array(dd.split(), dtype = 'float64')
        dd                    = self.regEx.shapedPowerWatt.search(self.acqusText)
        dd                    = self.acqusText[dd.span()[0]:]
        dd                    = dd[dd.find('\n')+1:]
        dd                    = dd[:dd.find('##$')]
        self.shapedPowerWatt  = np.array(dd.split(), dtype = 'float64')
        dd                    = self.regEx.spoal.search(self.acqusText)
        dd                    = self.acqusText[dd.span()[0]:]
        dd                    = dd[dd.find('\n')+1:]
        dd                    = dd[:dd.find('##$')]
        self.spoal            = np.array(dd.split(), dtype = 'float64')
        dd                    = self.regEx.spoffs.search(self.acqusText)
        dd                    = self.acqusText[dd.span()[0]:]
        dd                    = dd[dd.find('\n')+1:]
        dd                    = dd[:dd.find('##$')]
        self.spoffs           = np.array(dd.split(), dtype = 'float64')
        dd                    = self.regEx.inf.search(self.acqusText)
        dd                    = self.acqusText[dd.span()[0]:]
        dd                    = dd[dd.find('\n')+1:]
        dd                    = dd[:dd.find('##$')]
        self.inf              = np.array(dd.split(), dtype = 'float64')
        self.vcList           = self.regEx.vcList.findall(self.acqusText)[0]
        self.vdList           = self.regEx.vdList.findall(self.acqusText)[0]
        self.vpList           = self.regEx.vpList.findall(self.acqusText)[0]
        self.vaList           = self.regEx.vaList.findall(self.acqusText)[0]
        self.vtList           = self.regEx.vtList.findall(self.acqusText)[0]
        self.nuc1             = self.regEx.nuc1.findall(self.acqusText)[0]
        self.nuc2             = self.regEx.nuc2.findall(self.acqusText)[0]
        self.nuc3             = self.regEx.nuc3.findall(self.acqusText)[0]
        self.nuc4             = self.regEx.nuc4.findall(self.acqusText)[0]
        self.nuc5             = self.regEx.nuc5.findall(self.acqusText)[0]
        self.nuc6             = self.regEx.nuc6.findall(self.acqusText)[0]
        self.nuc7             = self.regEx.nuc7.findall(self.acqusText)[0]
        self.nuc8             = self.regEx.nuc8.findall(self.acqusText)[0]
        self.nusList          = self.regEx.nusList.findall(self.acqusText)[0]
        try:
            self.nusAmount        = float(self.regEx.nusAmount.findall(self.acqusText)[0])
            self.nusSeed          = int(self.regEx.nusSeed.findall(self.acqusText)[0])
            self.nusJsp           = int(self.regEx.nusJsp.findall(self.acqusText)[0])
            self.nusT2            = float(self.regEx.nusT2.findall(self.acqusText)[0][0])
        except:
            pass
            
        self.overFlow         = int(self.regEx.overFlow.findall(self.acqusText)[0])
        self.pynm             = self.regEx.pynm.findall(self.acqusText)[0]
        self.spcFrequency[0]  = float(self.regEx.bf1.findall(self.acqusText)[0])
        self.spcSFreq[0]      = float(self.regEx.sfo1.findall(self.acqusText)[0])
        self.spcNucleus[0]    = self.regEx.nuc1.findall(self.acqusText)[0]
        self.spcOffset[0]     = float(self.regEx.o1.findall(self.acqusText)[0])
        try:
            self.sw[1]            = float(self.regEx.sw.findall(self.acqu2sText)[0])
        except:
            pass

        try:
            self.sw_h[1]          = float(self.regEx.sw_h.findall(self.acqu2sText)[0])
        except:
            pass

        try:
            self.spcFrequency[1]  = float(self.regEx.bf1.findall(self.acqu2sText)[0])
        except:
            pass

        try:
            self.spcSFreq[1]      = float(self.regEx.sfo1.findall(self.acqu2sText)[0])
        except:
            pass

        try:
            self.spcNucleus[1]    = self.regEx.nuc1.findall(self.acqu2sText)[0]
        except:
            pass

        try:
            self.spcOffset[1]     = float(self.regEx.o1.findall(self.acqu2sText)[0])
        except:
            pass

        try:
            self.acqT0[1]         = int(self.regEx.acqT0.findall(self.acqu2sText)[0])
        except:
            pass

        try:
            self.nDataPoints[1]   = int(self.regEx.td.findall(self.acqu2sText)[0])
        except:
            pass

        try:
            self.nusTD            = int(self.regEx.nusTD.findall(self.acqu2sText)[0])
        except:
            pass

        try:
            self.fnMode           = int(self.regEx.fnMode.findall(self.acqu2sText)[0])
        except:
            pass

        try:
            self.acqT0[0]         = int(self.regEx.acqT0.findall(self.acqusText)[0])            
        except:
            pass

        try:
            self.sw[2]            = float(self.regEx.sw.findall(self.acqu3sText)[0]) if (len(self.regEx.sw.findall(self.acqu3sText))>0) else 0.0
            self.sw_h[2]          = float(self.regEx.sw_h.findall(self.acqu3sText)[0]) if (len(self.regEx.sw_h.findall(self.acqu3sText))>0) else 0.0
            self.nDataPoints[2]   = int(self.regEx.td.findall(self.acqu3sText)[0]) if (len(self.regEx.td.findall(self.acqu3sText))>0) else 0

        except:
            pass
        
        # end parseRegEx

    def read(self, spcDir):
        acqusName  = spcDir + os.sep + 'acqus'
        acqu2sName = spcDir + os.sep + 'acqu2s'
        acqu3sName = spcDir + os.sep + 'acqu3s'
        if(os.path.isfile(acqusName)):
            try:
                f = open(acqusName,"r")
                self.acqusText = f.read()
                f.close()

            except:
                f = open(acqusName,"r", encoding='latin-1')
                self.acqusText = f.read()
                f.close()

            
        if(os.path.isfile(acqu2sName)):
            try:
                f = open(acqu2sName,"r")
                self.acqu2sText = f.read()
                f.close()

            except:
                f = open(acqu2sName,"r", encoding='latin-1')
                self.acqu2sText = f.read()
                f.close()


        if(os.path.isfile(acqu3sName)):
            try:
                f = open(acqu3sName,"r")
                self.acqu3sText = f.read()
                f.close()

            except:
                f = open(acqu3sName,"r", encoding='latin-1')
                self.acqu3sText = f.read()
                f.close()

        
        self.parseRegEx()
        # end read
    

