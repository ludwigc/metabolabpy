"""
procpar regular expression class for processing parameters
"""

import re


class ProcProcparRegEx:

    def __init__(self):
        self.td = re.compile(r'\nnp \d+ \d+ \d+ \d+ \d+ \d+ \d+ \d+ \d+ \d+')
        self.nuc1 = re.compile(r'\ntn \d+ \d+ \d+ \d+ \d+ \d+ \d+ \d+ \d+ \d+')
        self.nuc2 = re.compile(r'\ndn \d+ \d+ \d+ \d+ \d+ \d+ \d+ \d+ \d+ \d+')
        self.nuc3 = re.compile(r'\ndn2 \d+ \d+ \d+ \d+ \d+ \d+ \d+ \d+ \d+ \d+')
        self.ni = re.compile(r'\nni ')
        self.ni2 = re.compile(r'\nni2 ')
        self.phase = re.compile(r'\nphase ')
        self.phase2 = re.compile(r'\nphase2 ')
        # end __init__
