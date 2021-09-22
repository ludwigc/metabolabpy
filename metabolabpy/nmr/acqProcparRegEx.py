"""
procpar regular expression class for acquisition parameters
"""

import re


class AcqProcparRegEx:

    def __init__(self):
        self.sfo1 = re.compile(r'\nsfrq ')
        self.sfo2 = re.compile(r'\ndfrq ')
        self.sfo3 = re.compile(r'\ndfrq2 ')
        self.o1 = re.compile(r'\ntof ')
        self.o2 = re.compile(r'\ndof ')
        self.o3 = re.compile(r'\ndof2 ')
        self.sw_h = re.compile(r'\nsw ')
        self.sw2_h = re.compile(r'\nsw1 ')
        self.sw3_h = re.compile(r'\nsw2 ')
        self.td = re.compile(r'\nnp ')
        self.ni = re.compile(r'\nni ')
        self.ni2 = re.compile(r'\nni2 ')
        self.phase = re.compile(r'\nphase ')
        self.phase2 = re.compile(r'\nphase2 ')
        self.transients = re.compile(r'\nnt ')
        self.steady_state_scans = re.compile(r'\nss ')
        self.nucleus = re.compile(r'\ntn ')
        self.temperature = re.compile(r'\ntemp ')
        self.nuc1 = re.compile(r'\ntn ')
        self.nuc2 = re.compile(r'\ndn ')
        self.nuc3 = re.compile(r'\ndn2 ')
        self.pul_prog_name = re.compile(r'\npslabel ')
        # end __init__
