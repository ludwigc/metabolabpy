'''
NMR data pre-processing

'''

import numpy as np


class NmrPreProc:
    
    def __init__(self):
        self.excludeStart          = np.array([])
        self.excludeEnd            = np.array([])
        self.segStart              = np.array([])
        self.segEnd                = np.array([])
        self.noiseThreshold        = 4                            # [%]
        self.bucketPoints          = 0
        self.bucketPPM             = 0
        self.compressBuckets       = False
        self.scaleSpc              = ""
        self.varianceStabilisation = ""
        self.varLambda             = 1e-5
        self.varX0                 = 0.0
        self.pName                 = ""
        self.fName                 =""
        self.samplesInRows         = True
        self.plotSelect            = np.array([])
        self.classSelect           = np.array([""])
        # end __init__

    def __str__(self):
        strStr = "NMR data pre-processing"
        return strStr
        # end __str__

    def init(self, nspc):
        self.plotSelect     = np.arange(nspc)
        self.classSelect    = np.empty(nspc, dtype = 'str')
        for k in range(nspc):
            self.classSelect[k] = "1"
            
        return "pre-processing initialised"
        # end init
        
        
        
    
