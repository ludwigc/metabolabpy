"""
NMR HSQC multiplet analysis

"""

import numpy as np
from metabolabpy.nmr import hsqcData  # pragma: no cover
import os


class NmrHsqc:

    def __init__(self):
        self.range_h = 0.1           # [ppm]
        self.range_c = 0.7           # [ppm]
        self.threshold = 0.05        # [%]
        self.j_cc = 60               # [Hz]
        self.j_ch = 145              # [Hz]
        self.n_max = 4
        self.use_splitting = True
        self.tilt_hsqc = False
        self.display_library_shift = True
        self.display_peaks_of_metabolite = True
        self.highlight_double_assignments = True
        self.display_overlay_shifts = True
        self.metabolite_list = []
        self.metabolite_information = ''
        self.hsqc_data = {}
        self.co_hsqc = False
        self.cur_peak = 0
        self.cur_metabolite = ''
        self.pick_local_opt = True
        self.autoscale_j = True
        self.j_scale = -1
        self.assigned_peaks = []

        # end __init__

    def __str__(self):  # pragma: no cover
        str_str = "NMR HSQC multiplet analysis"
        return str_str
        # end __str__

    def read_metabolite_information(self, metabolite_name=''):
        print('====================')
        print(metabolite_name)
        print('====================')
        if len(metabolite_name) == 0:
            return

        file_name = os.path.join(os.path.dirname(__file__), 'metabolites', metabolite_name + '.mlInfo')
        fid = open(file_name, 'r')
        self.metabolite_information = fid.read()
        fid.close()
        # end read_metabolite_information

    def set_metabolite_list(self):
        met_list = next(os.walk(os.path.join(os.path.dirname(__file__), 'metabolites')), (None, None, []))[2]
        met_list = sorted(met_list)
        self.metabolite_list = []
        for k in met_list:
            idx = k.find('.mlInfo')
            if idx > -1:
                self.metabolite_list.append(k[:idx])

        # end set_metabolite_list

    def set_metabolite_information(self, metabolite_name='', metabolite_information=''):
        if len(metabolite_name) == 0 or len(metabolite_information) == 0:
            return

        if metabolite_name not in self.hsqc_data.keys():
            self.hsqc_data[metabolite_name] = hsqcData.HsqcData()
            self.hsqc_data[metabolite_name].init_data(metabolite_information)

        # end set_metabolite_information
