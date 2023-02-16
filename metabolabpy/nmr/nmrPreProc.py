"""
NMR data pre-processing

"""

import numpy as np
from metabolabpy.nmr import nmrConfig  # pragma: no cover
import os
import darkdetect


class NmrPreProc:

    def __init__(self):
        self.exclude_start = np.array([])
        self.exclude_end = np.array([])
        self.seg_start = np.array([])
        self.seg_end = np.array([])
        self.compress_start = np.array([])
        self.compress_end = np.array([])
        self.seg_align_ref_spc = 1  # use first spectrum
        self.noise_threshold = 4.0  # times std of noise region
        self.noise_start = 10.0
        self.noise_end = 10.5
        self.bucket_points = 0
        self.bucket_ppm = 0.005
        self.compress_buckets = False
        self.scale_spc = ""
        self.variance_stabilisation = ""
        self.auto_scaling = False
        self.pareto_scaling = True
        self.g_log_transform = False
        self.var_lambda = 1e-5
        self.var_y0 = 0.0
        self.p_name = ""
        self.f_name = ""
        self.samples_in_rows = True
        self.plot_select = np.array([])
        self.class_select = np.array([])
        self.plot_colours = []
        self.pre_proc_fill = False
        self.alpha = 1.0
        self.colour = 'gray'
        self.th_colour = 'red'
        self.th_line_width = 2.0
        self.flag_exclude_region = False
        self.flag_segmental_alignment = False
        self.flag_noise_filtering = False
        self.flag_bucket_spectra = False
        self.flag_compress_buckets = False
        self.flag_scale_spectra = False
        self.flag_variance_stabilisation = False
        self.flag_export_data_set = False
        self.scale_pqn = True
        self.scale_spectra_ref_spc = 0  # use mean spectra
        self.preserve_overall_scale = False
        self.std_val = 0.0
        self.export_path_name = os.getcwd()
        self.export_file_name = "csvExport.csv"
        self.export_excel_path = os.getcwd()
        self.export_excel = "excelExport.xlsx"
        self.export_metabo_analyst_path = os.getcwd()
        self.export_metabo_analyst = "metaboAnalystExport.csv"
        self.export_r_dolphin_path = os.getcwd()
        self.export_r_dolphin = "rDolphinExport.csv"
        self.export_batman_path = os.getcwd()
        self.export_batman = "NMRdata.txt"
        self.export_bruker_path = os.getcwd()
        self.export_bruker = "brukerExport"
        self.export_delimiter_tab = True
        self.export_character = ","
        self.export_samples_in_rows_cols = 0
        self.export_method = 0
        self.avoid_negative_values = False
        self.cf = nmrConfig.NmrConfig()
        self.cf.read_config()
        self.int1 = 0.0
        self.int2 = 0.0
        self.int3 = 0.0
        self.spc_scale = []
        # end __init__

    def __str__(self):  # pragma: no cover
        str_str = "NMR data pre-processing"
        return str_str
        # end __str__

    def init(self, nspc):
        self.plot_select = np.arange(nspc)
        self.class_select = np.empty(nspc, dtype='str')
        for k in range(nspc):
            self.class_select[k] = "1"

        self.init_plot_colours()
        return "pre-processing initialised"
        # end init

    def init_plot_colours(self):
        self.cf.read_config()
        if self.cf.mode == 'dark' or (self.cf.mode == 'system' and darkdetect.isDark()):
            self.int1 = 1.0
            self.int2 = 0.6
            self.int3 = 0.3
        else:
            self.int1 = 0.4
            self.int2 = 0.8
            self.int3 = 0.5

        int1 = self.int1
        int2 = self.int2
        int3 = self.int3
        if self.cf.mode == 'dark' or (self.cf.mode == 'system' and darkdetect.isDark()):
            self.plot_colours = [(int1, int1, 0.0),
                                 (0.0, int1, int1),
                                 (int1, 0.0, int1),
                                 (int2, int2, int1),
                                 (int1, int2, int2),
                                 (int2, int1, int2),
                                 (int1, int1, int3),
                                 (int2, int3, int3),
                                 (int3, int2, int3),
                                 (int3, int2, int2),
                                 (int2, int2, int3),
                                 (int2, int3, int2)]

        else:
            self.plot_colours = [(0.0, 0.0, int1),
                                 (int1, 0.0, 0.0),
                                 (0.0, int1, 0.0),
                                 (0.0, int1, int1),
                                 (int1, int1, 0.0),
                                 (int1, 0.0, int1),
                                 (int3, int3, int2),
                                 (int2, int3, int3),
                                 (int3, int2, int3),
                                 (int3, int2, int2),
                                 (int2, int2, int3),
                                 (int2, int3, int2)]

        # end initplot_colours

    def set_alpha(self, alpha):
        self.alpha = alpha
        return "Alpha set"
        # end set_alpha
