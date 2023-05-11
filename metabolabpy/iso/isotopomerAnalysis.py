import os
import pandas as pd
import numpy as np
import scipy as sp
import os
import math
import time
import metabolabpy.__init__ as ml_version


class IsotopomerAnalysis:

    def __init__(self):
        self.nmr_multiplets = pd.DataFrame()
        self.ver = ml_version.__version__
        # end __init__

    def __str__(self):  # pragma: no cover
        r_string = '______________________________________________________________________________________\n'
        r_string += '\nMetaboLabPy Isotopomer Data Analysis (v. ' + self.ver + ')\n'
        r_string += '______________________________________________________________________________________\n\n'
        return r_string
        # end __str__

    def read_hsqc_multiplets(self, file_name=''):
        if len(file_name) == 0:
            return

        self.nmr_multiplets = pd.read_excel(file_name, sheet_name=None, keep_default_na=False)
        return
    # end read_hsqc_multiplets
