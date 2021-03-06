'''
acqus regular expression class
'''

import re

class AcqRegEx:
    
    def __init__(self):
        self.sfo1             = re.compile(r'##\$SFO1= (\d+\.\d+)')
        self.sfo2             = re.compile(r'##\$SFO2= (\d+\.\d+)')
        self.sfo3             = re.compile(r'##\$SFO3= (\d+\.\d+)')
        self.sfo4             = re.compile(r'##\$SFO4= (\d+\.\d+)')
        self.sfo5             = re.compile(r'##\$SFO5= (\d+\.\d+)')
        self.sfo6             = re.compile(r'##\$SFO6= (\d+\.\d+)')
        self.sfo7             = re.compile(r'##\$SFO7= (\d+\.\d+)')
        self.sfo8             = re.compile(r'##\$SFO8= (\d+\.\d+)')
        self.bf1              = re.compile(r'##\$BF1= (\d+\.\d+)')
        self.bf2              = re.compile(r'##\$BF2= (\d+\.\d+)')
        self.bf3              = re.compile(r'##\$BF3= (\d+\.\d+)')
        self.bf4              = re.compile(r'##\$BF4= (\d+\.\d+)')
        self.bf5              = re.compile(r'##\$BF5= (\d+\.\d+)')
        self.bf6              = re.compile(r'##\$BF6= (\d+\.\d+)')
        self.bf7              = re.compile(r'##\$BF7= (\d+\.\d+)')
        self.bf8              = re.compile(r'##\$BF8= (\d+\.\d+)')
        self.o1               = re.compile(r'##\$O1= (\d+\.?\d?\d?\d?\d?\d?\d?)')
        self.o2               = re.compile(r'##\$O2= (\d+\.?\d?\d?\d?\d?\d?\d?)')
        self.o3               = re.compile(r'##\$O3= (\d+\.?\d?\d?\d?\d?\d?\d?)')
        self.o4               = re.compile(r'##\$O4= (\d+\.?\d?\d?\d?\d?\d?\d?)')
        self.o5               = re.compile(r'##\$O5= (\d+\.?\d?\d?\d?\d?\d?\d?)')
        self.o6               = re.compile(r'##\$O6= (\d+\.?\d?\d?\d?\d?\d?\d?)')
        self.o7               = re.compile(r'##\$O7= (\d+\.?\d?\d?\d?\d?\d?\d?)')
        self.o8               = re.compile(r'##\$O8= (\d+\.?\d?\d?\d?\d?\d?\d?)')
        self.sw               = re.compile(r'##\$SW= (\d+\.\d+)')
        self.sw_h             = re.compile(r'##\$SW_h= (\d+\.?\d?\d?\d?\d?\d?\d?\d?\d?\d?\d?\d?\d?\d?\d?\d?\d?\d?\d?)')
        self.td               = re.compile(r'##\$TD= (\d+)')
        self.decim            = re.compile(r'##\$DECIM= (\d+)')
        self.dspfvs           = re.compile(r'##\$DSPFVS= (\d+)')
        self.grpdly           = re.compile(r'##\$GRPDLY= (-?\d+(\.\d+)?)')
        self.byteOrder        = re.compile(r'##\$BYTORDA= (\d+)')
        self.aqMode           = re.compile(r'##\$AQ_mod= (\d+)')
        self.digMod           = re.compile(r'##\$DIGMOD= (\d+)')
        self.transients       = re.compile(r'##\$NS= (\d+)')
        self.steadyStateScans = re.compile(r'##\$DS= (\d+)')
        self.relaxationDelay  = re.compile(r'##\$RD= (\d+\.?\d?)')
        self.spinRate         = re.compile(r'##\$MASR= (\d+)')
        self.pulProg          = re.compile(r'##\$PULPROG= (.+)')
        self.aunm             = re.compile(r'##\$AUNM= (.+)')
        self.autopos          = re.compile(r'##\$AUTOPOS= (.+)')
        self.nucleus          = re.compile(r'##\$NUC\d= (.+)')
        self.nucleusIndex     = re.compile(r'##\$NUC(\d)=.+')
        self.instrument       = re.compile(r'##\$INSTRUM= (.+)')
        self.dataType         = re.compile(r'##\$DTYPA= (\d+)')
        self.solvent          = re.compile(r'##\$SOLVENT= (.+)')
        self.probe            = re.compile(r'##\$PROBHD= (.+)')
        self.title            = re.compile(r'##TITLE= (.+), (.+)(\t)?(.+)')
        self.origin           = re.compile(r'##ORIGIN= (.+)')
        self.owner            = re.compile(r'##OWNER= (.+)')
        self.metaInfo         = re.compile(r'\$\$ (.+)')
        self.temperature      = re.compile(r'##\$TE= (\d+\.?\d?)')
        self.cnst             = re.compile(r'##\$CNST= \(\d+..\d+\)')
        self.delay            = re.compile(r'##\$D= \(\d+..\d+\)')
        self.cpdProg          = re.compile(r'##\$CPDPRG= \(\d+..\d+\)')
        self.gpName           = re.compile(r'##\$GPNAM= \(\d+..\d+\)')
        self.gpx              = re.compile(r'##\$GPX= \(\d+..\d+\)')
        self.gpy              = re.compile(r'##\$GPY= \(\d+..\d+\)')
        self.gpz              = re.compile(r'##\$GPZ= \(\d+..\d+\)')
        self.increments       = re.compile(r'##\$IN= \(\d+..\d+\)')
        self.pulse            = re.compile(r'##\$P= \(\d+..\d+\)')
        self.pcpd             = re.compile(r'##\$PCPD= \(\d+..\d+\)')
        self.powerLevel       = re.compile(r'##\$PL= \(\d+..\d+\)')
        self.powerLevelWatt   = re.compile(r'##\$PLW= \(\d+..\d+\)')
        self.powerLevelMax    = re.compile(r'##\$PLWMAX= \(\d+..\d+\)')
        self.shapedPower      = re.compile(r'##\$SP= \(\d+..\d+\)')
        self.shapedPowerWatt  = re.compile(r'##\$SPW= \(\d+..\d+\)')
        self.spoal            = re.compile(r'##\$SPOAL= \(\d+..\d+\)')
        self.spoffs           = re.compile(r'##\$SPOFFS= \(\d+..\d+\)')
        self.vcList           = re.compile(r'##\$VCLIST= (.+)')
        self.vdList           = re.compile(r'##\$VDLIST= (.+)')
        self.vpList           = re.compile(r'##\$VPLIST= (.+)')
        self.vaList           = re.compile(r'##\$VALIST= (.+)')
        self.vtList           = re.compile(r'##\$VTLIST= (.+)')
        self.nuc1             = re.compile(r'##\$NUC1= (.+)')
        self.nuc2             = re.compile(r'##\$NUC2= (.+)')
        self.nuc3             = re.compile(r'##\$NUC3= (.+)')
        self.nuc4             = re.compile(r'##\$NUC4= (.+)')
        self.nuc5             = re.compile(r'##\$NUC5= (.+)')
        self.nuc6             = re.compile(r'##\$NUC6= (.+)')
        self.nuc7             = re.compile(r'##\$NUC7= (.+)')
        self.nuc8             = re.compile(r'##\$NUC8= (.+)')
        self.nusList          = re.compile(r'##\$NUSLIST= (.+)')
        self.nusAmount        = re.compile(r'##\$NusAMOUNT= (\d+)')
        self.nusSeed          = re.compile(r'##\$NusSEED= (\d+)')
        self.nusJsp           = re.compile(r'##\$NusJSP= (\d+)')
        self.nusT2            = re.compile(r'##\$NusT2= (\d+(\.\d+)?)')
        self.nusTD            = re.compile(r'##\$NusTD= (\d+)')
        self.overFlow         = re.compile(r'##\$OVERFLW= (\d+)')
        self.pynm             = re.compile(r'##\$PYNM= (.+)')
        self.acqT0            = re.compile(r'##\$ACQT0= (\d+)')
        self.fnMode           = re.compile(r'##\$FnMODE= (\d+)')
        self.inf              = re.compile(r'##\$INF= \(\d+..\d+\)')
        # end __init__
