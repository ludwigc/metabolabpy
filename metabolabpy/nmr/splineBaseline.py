"""
Spline baseline correction parameters
"""

import numpy as np


class SplineBaseline:

    def __init__(self):
        self.linear_spline = 200
        self.average_points = 50
        self.baseline_points = []
        self.baseline_points_pts = []
        self.baseline_values = []
        # end __init__

    def __str__(self):  # pragma: no cover
        return "Spline Baseline Correction Parameters"
        # end __str__


