'''
interactive NMR spectrum phase correction
'''

import numpy as np

class PhCorr:
    
    def __init__(self):
        self.start          = 0.0
        self.maxPh0         = 1440.0
        self.maxPh1         = 1440.0
        self.ph0            = 0.0
        self.ph1            = 0.0
        self.pivot          = 0.0
        self.pivPoints      = 0
        self.offset         = 0.0
        self.horOffset      = 0.0
        self.spcMax         = 0.0
        self.scale          = 1.0
        self.spc            = np.array([[]], dtype = 'complex')
        self.xData          = 0.0
        self.yData          = 0.0
        # end __init__

    def __str__(self):
        return "Interactive NMR spectrum phase correction"
        # end __str__
        

