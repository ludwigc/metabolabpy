"""
procs regular expression class
"""

import re


class ProcRegEx:

    def __init__(self):
        self.aunmp = re.compile(r'##\$AUNMP= (.+)')
        self.offset = re.compile(r'##\$OFFSET= (\d+(\.\d+)?)')
        self.ph0 = re.compile(r'##\$PHC0= (-?\d+(\.\d+)?)')
        self.ph1 = re.compile(r'##\$PHC1= (-?\d+(\.\d+)?)')
        self.n_points = re.compile(r'##\$SI= (\d+)')
        self.lb = re.compile(r'##\$LB= (-?\d+(\.\d+)?)')
        self.gb = re.compile(r'##\$GB= (-?\d+(\.\d+)?)')
        self.ssb = re.compile(r'##\$SSB= (-?\d+(\.\d+)?)')
        self.wdw = re.compile(r'##\$WDW= (\d+)')
        self.axis_nucleus = re.compile(r'##\$AXNUC= (.+)')
        self.si = re.compile(r'##\$SI= (\d+)')
        self.ft_mod = re.compile(r'##\$FT_mod= (\d+)')
        self.nc_proc = re.compile(r'##\$NC_proc= (-?\d+)')
        self.sf = re.compile(r'##\$SF= (\d+(\d.\d+)?)')
        self.stsr = re.compile(r'##\$STSR= (\d+(\d+)?)')
        self.stsi = re.compile(r'##\$STSI= (\d+(\d+)?)')
        # end __init__
