'''
NMR spectrum acquisition parameters
'''

import numpy as np
import os
from metabolabpy.nmr import acqRegEx
from metabolabpy.nmr import acqProcparRegEx

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
        self.regExVarian      = acqProcparRegEx.AcqProcparRegEx()
        self.manufacturer     = ''
        self.ni               = 0
        self.ni2              = 0
        self.np               = 1
        self.np2              = 1
        self.phase            = np.array([], dtype = 'int')
        self.phase2           = np.array([], dtype = 'int')
        self.autopos          = ''
        # end __init__
        
    def __str__(self): # pragma: no cover
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
        if self.acqusText.find("$GRPDLY=") > -1:
            self.groupDelay       = float(self.regEx.grpdly.findall(self.acqusText)[0][0])
            if self.groupDelay < 0:
                self.groupDelay = 0

        self.byteOrder        = int(self.regEx.byteOrder.findall(self.acqusText)[0])
        self.aqMode           = int(self.regEx.aqMode.findall(self.acqusText)[0])
        self.digMod           = int(self.regEx.digMod.findall(self.acqusText)[0])
        self.transients       = int(self.regEx.transients.findall(self.acqusText)[0])
        self.steadyStateScans = int(self.regEx.steadyStateScans.findall(self.acqusText)[0])
        self.relaxationDelay  = float(self.regEx.relaxationDelay.findall(self.acqusText)[0]) if (len(self.regEx.relaxationDelay.findall(self.acqusText))>0) else 0.0
        self.spinRate         = int(self.regEx.spinRate.findall(self.acqusText)[0])
        self.pulProgName      = self.regEx.pulProg.findall(self.acqusText)[0]
        self.aunm             = self.regEx.aunm.findall(self.acqusText)[0]
        self.autopos          = self.regEx.autopos.findall(self.acqusText)[0]
        self.instrument       = self.regEx.instrument.findall(self.acqusText)[0]
        self.dataType         = int(self.regEx.dataType.findall(self.acqusText)[0])
        self.solvent          = self.regEx.solvent.findall(self.acqusText)[0]
        self.probe            = self.regEx.probe.findall(self.acqusText)[0]
        self.probe            = self.probe.replace('<','')
        self.probe            = self.probe.replace('>','')
        self.title            = self.regEx.title.findall(self.acqusText)[0]
        self.origin           = self.regEx.origin.findall(self.acqusText)[0]
        self.owner            = self.regEx.owner.findall(self.acqusText)[0]
        try:
            self.metaInfo         = self.regEx.metaInfo.findall(self.acqusText)[0]
        except:
            pass

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
        if self.acqusText.find("$PLW=") > -1:
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
        if self.acqusText.find("$SPW=") > -1:
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
        if self.acqusText.find("$INF=") > -1:
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
        if self.acqusText.find("$NUSLIST=") > -1:
            self.nusList          = self.regEx.nusList.findall(self.acqusText)[0]
            try:
                self.nusAmount        = float(self.regEx.nusAmount.findall(self.acqusText)[0])
                self.nusSeed          = int(self.regEx.nusSeed.findall(self.acqusText)[0])
                self.nusJsp           = int(self.regEx.nusJsp.findall(self.acqusText)[0])
                self.nusT2            = float(self.regEx.nusT2.findall(self.acqusText)[0][0])
            except:
                pass


        self.overFlow         = int(self.regEx.overFlow.findall(self.acqusText)[0])
        if self.acqusText.find("$PYNM=") > -1:
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

    def parseRegExVarian(self):
        dd        = self.regExVarian.sfo1.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd        = self.acqusText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            self.sfo1 = float(dd)

        dd        = self.regExVarian.sfo2.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd        = self.acqusText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            self.sfo2 = float(dd)

        dd        = self.regExVarian.sfo3.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd        = self.acqusText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            self.sfo3 = float(dd)

        dd        = self.regExVarian.o1.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd        = self.acqusText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            self.o1 = float(dd)

        dd        = self.regExVarian.o2.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd        = self.acqusText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            self.o2 = float(dd)

        dd        = self.regExVarian.o3.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd        = self.acqusText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            self.o3 = float(dd)

        self.bf1 = self.sfo1 - self.o1/1000.0
        self.bf2 = self.sfo2 - self.o2/1000.0
        self.bf3 = self.sfo3 - self.o3/1000.0
        self.o1 = 0
        self.o2 = 0
        self.o3 = 0
        dd        = self.regExVarian.sw_h.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd        = self.acqusText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            self.sw_h[0] = float(dd)
            self.sw[0]   = self.sw_h[0]/self.sfo1

        dd        = self.regExVarian.sw2_h.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd        = self.acqusText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            self.sw_h[1] = float(dd)
            self.sw[1]   = self.sw_h[1]/self.sfo2

        dd        = self.regExVarian.sw3_h.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd        = self.acqusText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            dd2 = self.regExVarian.sw2_h.search(self.acqusText)
            if hasattr(dd2, 'span'):
                self.sw_h[2] = float(dd)
                self.sw[2]   = self.sw_h[2]/self.sfo3
            else:
                self.sw_h[1] = float(dd)
                self.sw[1]   = self.sw_h[1]/self.sfo3


        dd        = self.regExVarian.td.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd        = self.acqusText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            self.nDataPoints[0] = int(dd)  # int(int(dd)/2)

        dd        = self.regExVarian.phase.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd      = self.acqusText[dd.span()[0] + 1:]
            dd      = dd[dd.find('\n') + 1:]
            dd      = dd[:dd.find(' ')]
            self.np = int(dd)

        dd        = self.regExVarian.phase2.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd       = self.acqusText[dd.span()[0] + 1:]
            dd       = dd[dd.find('\n') + 1:]
            dd       = dd[:dd.find(' ')]
            self.np2 = int(dd)

        dd        = self.regExVarian.ni.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd        = self.acqusText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            self.ni = int(dd)

        self.ni = max(self.ni, 1)
        dd        = self.regExVarian.ni2.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd        = self.acqusText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            self.ni2 = int(dd)

        self.ni2 = max(self.ni2, 1)
        dd        = self.regExVarian.phase.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd         = self.acqusText[dd.span()[0] + 1:]
            dd         = dd[dd.find('\n') + 1:]
            dd         = dd[dd.find(' '):dd.find('\n')]
            self.phase = np.array(dd.split(), dtype = 'int')

        dd        = self.regExVarian.phase2.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd         = self.acqusText[dd.span()[0] + 1:]
            dd         = dd[dd.find('\n') + 1:]
            dd         = dd[dd.find(' '):dd.find('\n')]
            self.phase2 = np.array(dd.split(), dtype = 'int')

        dd        = self.regExVarian.transients.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd        = self.acqusText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            self.transients = int(dd)

        dd        = self.regExVarian.steadyStateScans.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd        = self.acqusText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            self.steadyStateScans = int(dd)

        dd        = self.regExVarian.nucleus.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd        = self.acqusText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            dd        = dd.replace(" ", "")
            dd        = dd.replace("\"", "")
            self.spcNucleus = dd

        dd        = self.regExVarian.pulProgName.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd        = self.acqusText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            dd        = dd.replace(" ", "")
            dd        = dd.replace("\"", "")
            self.pulProgName = dd

        dd        = self.regExVarian.temperature.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd        = self.acqusText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            self.temperature = float(dd) + 273.15

        dd        = self.regExVarian.nuc1.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd        = self.acqusText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            dd        = dd.replace(" ", "")
            dd        = dd.replace("\"", "")
            self.nuc1 = dd

        dd        = self.regExVarian.nuc2.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd        = self.acqusText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            dd        = dd.replace(" ", "")
            dd        = dd.replace("\"", "")
            self.nuc2 = dd

        dd        = self.regExVarian.nuc3.search(self.acqusText)
        if hasattr(dd, 'span'):
            dd        = self.acqusText[dd.span()[0] + 1:]
            dd        = dd[dd.find('\n') + 1:]
            dd        = dd[dd.find(' '):dd.find('\n')]
            dd        = dd.replace(" ", "")
            dd        = dd.replace("\"", "")
            self.nuc3 = dd

        # end parseRegExVarian

    def read(self, spcDir):
        acqusName    = spcDir + os.sep + 'acqus'
        acqu2sName   = spcDir + os.sep + 'acqu2s'
        acqu3sName   = spcDir + os.sep + 'acqu3s'
        procparName  = spcDir + os.sep + 'procpar'
        procparName2 = spcDir + os.sep + 'PROCPAR'
        if(os.path.isfile(acqusName)):
            try:
                f = open(acqusName,"r")
                self.acqusText = f.read()
                f.close()

            except:
                f = open(acqusName,"r", encoding='latin-1')
                self.acqusText = f.read()
                f.close()

            self.manufacturer = 'Bruker'
            
        if(os.path.isfile(acqu2sName)):
            try:
                f = open(acqu2sName,"r")
                self.acqu2sText = f.read()
                f.close()

            except:
                f = open(acqu2sName,"r", encoding='latin-1')
                self.acqu2sText = f.read()
                f.close()

            self.manufacturer = 'Bruker'


        if(os.path.isfile(acqu3sName)):
            try:
                f = open(acqu3sName,"r")
                self.acqu3sText = f.read()
                f.close()

            except:
                f = open(acqu3sName,"r", encoding='latin-1')
                self.acqu3sText = f.read()
                f.close()

            self.manufacturer = 'Bruker'

        if (os.path.isfile(procparName)):
            f = open(procparName, "r")
            self.acqusText = f.read()
            f.close()
            self.byteOrder = 0
            self.manufacturer = 'Varian'

        if (os.path.isfile(procparName2)):
            f = open(procparName2, "r")
            self.acqusText = f.read()
            f.close()
            self.byteOrder = 0
            self.manufacturer = 'Varian'

        if self.manufacturer == 'Bruker':
            self.parseRegEx()
            if(self.groupDelay == 0.0):
                self.setGroupDelay()


        if self.manufacturer == 'Varian':
            self.parseRegExVarian()
            self.nDataPoints[1] = self.ni*self.ni2*len(self.phase)*len(self.phase2)*2
            self.groupDelay = 0.0
            self.decim = 0
            self.dspfvs = 0

        # end read
    
    def setGroupDelay(self):
        decims = np.array([2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048])
        dspfirm10 = np.array([179, 201, 533, 709, 1097, 1449, 2225, 2929, 4481, 5889, 8993, 11809, 18017, 23649, 36065, 47329, 72161, 94689, 144353, 189409, 288737])
        dspfirm11 = np.array([184, 219, 384, 602, 852, 1668, 2312, 3368, 4656, 6768, 9344, 13568, 18560, 27392, 36992, 55040, 73856, 110336, 147584, 220928, 295040])
        dspfirm12 = np.array([184, 219, 384, 602, 852, 1668, 2292, 3368, 4616, 6768, 9264, 13568, 18560, 27392, 36992, 55040, 73856, 110336, 147584, 220928, 295040])
        dspfirm13 = np.array([11, 17, 23, 35, 47, 71, 95, 143, 191, 287, 383, 575])
        dspfirm14 = np.array([60, 90, 118, 179, 244, 360, 492, 724, 980, 1444, 1958, 2886, 3912, 5768, 7820, 11532])
        dspfirm15 = np.array([0, 0, 58, 152, 202, 318, 418, 642, 842, 1290, 1690, 2586, 3386])
        dspfirm = [dspfirm10, dspfirm11, dspfirm12, dspfirm13, dspfirm14, dspfirm15]
        addr = np.where(decims == self.decim)[0][0]
        dly = dspfirm[self.dspfvs-10][addr]
        if self.decim == 3 and self.dspfvs== 15 and self.sw_h > 104000.0:
            self.groupDelay = 55.0
        else:
            self.groupDelay = dly / self.decim / 2.0

        # end setGroupDelay


