"""
interactive NMR spectrum phase correction
"""

import numpy as np  # pragma: no cover


class PhCorr:  # pragma: no cover

    def __init__(self):
        self.start = 0.0
        self.max_ph0 = 1440.0
        self.max_ph1 = 1440.0
        self.ph0 = 0.0
        self.ph1 = 0.0
        self.pivot = 0.0
        self.piv_points = 0
        self.offset = 0.0
        self.hor_offset = 0.0
        self.spc_max = 0.0
        self.scale = 1.0
        self.spc = np.array([[]], dtype='complex')
        self.spc2 = np.array([[]], dtype='complex')
        self.x_data = 0.0
        self.y_data = 0.0
        self.spc_row = []
        self.spc_col = []
        self.spc_row_pts = []
        self.spc_col_pts = []
        self.pivot2d = np.array([0.0, 0.0])
        self.pivot_points2d = np.array([-1, -1])
        self.ph0_2d = np.array([0.0, 0.0])
        self.ph1_2d = np.array([0.0, 0.0])
        self.ppm = np.array([], dtype='float')
        self.dim = 0
        self.n_dims = 1
        # end __init__

    def __str__(self):
        return "Interactive NMR spectrum phase correction"
        # end __str__
