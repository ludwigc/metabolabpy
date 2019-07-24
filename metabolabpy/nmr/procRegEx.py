'''
procs regular expression class
'''

import re

class ProcRegEx:
    
    def __init__(self):
        self.aunmp            = re.compile(r'##\$AUNMP= (.+)')
        self.offset           = re.compile(r'##\$OFFSET= (\d+(\.\d+)?)')
        self.ph0              = re.compile(r'##\$PHC0= (-?\d+(\.\d+)?)')
        self.ph1              = re.compile(r'##\$PHC1= (-?\d+(\.\d+)?)')
        self.nPoints          = re.compile(r'##\$SI= (\d+)')
        self.lb               = re.compile(r'##\$LB= (-?\d+(\.\d+)?)')
        self.gb               = re.compile(r'##\$GB= (-?\d+(\.\d+)?)')
        self.ssb              = re.compile(r'##\$SSB= (-?\d+(\.\d+)?)')
        self.wdw              = re.compile(r'##\$WDW= (\d+)')
        self.axisNucleus      = re.compile(r'##\$AXNUC= (.+)')
        # end __init__

