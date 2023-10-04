"""
NMR spectrum processing parameters
"""

import numpy as np
import os
from metabolabpy.nmr import procRegEx
from metabolabpy.nmr import procProcparRegEx
import math


class ProcPars:

    def __init__(self):
        self.autobaseline = False
        self.procs_text = str('')
        self.proc2s_text = str('')
        self.proc3s_text = str('')
        self.tilt = False
        self.symj = False
        self.ph0 = np.array([0.0, 0.0, 0.0])
        self.ph1 = np.array([0.0, 0.0, 0.0])
        self.ph_corr = np.array([1, 1, 0], dtype='int')
        self.ref_shift = np.array([0.0, 0.0, 0.0])
        self.ref_point = np.array([0, 0, 0])
        self.n_points = np.array([0, 0, 0])
        self.pivot = np.array([0, 0, 0])
        self.lb = np.array([0.3, 0.0, 0.0])
        self.gb = np.array([0.0, 0.0, 0.0])
        self.ssb = np.array([0.0, 0.0, 0.0])
        self.axis_nucleus = np.array(['      ', '      ', '      '], dtype='str')
        self.aunmp = str('')
        self.poly_order = 24
        self.water_suppression = 0
        self.gibbs = np.array([True, True, False])
        self.conv_extrapolation_size = np.array([32, 0, 0])
        self.conv_window_size = np.array([32, 0, 0])
        self.window_type = np.array([1, 4, 0], dtype='int')
        self.conv_window_type = np.array([0, 0, 0], dtype='int')
        self.reg_ex = procRegEx.ProcRegEx()
        self.reg_ex_varian = procProcparRegEx.ProcProcparRegEx()
        self.sw_h = np.array([0.0, 0.0, 0.0])
        self.sf = np.array([0.0, 0.0, 0.0])
        self.fid_offset_correction = 0
        self.strip_start = 0
        self.strip_end = 0
        self.phase_inversion = True
        self.mult_factor = [0, 0]
        self.invert_matrix = [False, False, False]
        self.window_functions = {
            0: "None",
            1: "Exponential",
            2: "Gaussian",
            3: "Sine",
            4: "QSine",
            5: "SEM",
            "None": 0,
            "Exponential": 1,
            "Gaussian": 2,
            "Sine": 3,
            "QSine": 4,
            "SEM": 5
        }
        self.phase_correction = {
            0: "None",
            1: "Manual",
            2: "Automatic",
            "None": 0,
            "Manual": 1,
            "Automatic": 2
        }
        self.win_type = {
            0: "Gaussian",
            1: "Sine",
            "Gaussian": 0,
            "Sine": 1
        }
        self.gibbs_p = {
            0: False,
            1: True,
            False: 0,
            True: 1
        }
        self.water_supp = {
            0: "None",
            1: "Conv",
            2: "Poly",
            3: "WaveWat",
            "None": 0,
            "Conv": 1,
            "Poly": 2,
            "WaveWat": 3
        }
        self.data_type = ''
        self.autobaseline_alg = 'rolling_ball'
        self.autobaseline_lam = 1e5
        self.autobaseline_max_iter = 50
        self.autobaseline_alpha = 0.1
        self.autobaseline_beta = 10
        self.autobaseline_gamma = 15
        self.autobaseline_beta_mult = 0.98
        self.autobaseline_gamma_mult = 0.94
        self.autobaseline_half_window = 4096
        self.autobaseline_quantile = 0.3
        self.autobaseline_poly_order = 4
        self.autobaseline_smooth_half_window = 16
        self.autobaseline_add_ext = 2
        self.ww_start = 9
        self.ww_zf = 16
        self.ww_wavelet_type = 'db'
        self.ww_wavelet_type_number = '10'
        self.ww_wavelet_number = 11
        self.wavelet_names = ['bior', 'coif', 'db', 'dmey', 'haar', 'rbio', 'sym']
        self.wavelet_numbers = {
            'bior': ['1.1', '1.3', '1.5', '2.2', '2.4', '2.6', '2.8', '3.1', '3.3', '3.5', '3.7', '3.9', '4.4', '5.5',
                     '6.8'], 'cgau': ['1', '2', '3', '4', '5', '6', '7', '8'],
            'coif': ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17'],
            'db': ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18',
                   '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35',
                   '36', '37', '38'], 'dmey': [''], 'haar': [''],
            'rbio': ['1.1', '1.3', '1.5', '2.2', '2.4', '2.6', '2.8', '3.1', '3.3', '3.5', '3.7', '3.9', '4.4', '5.5',
                     '6.8'],
            'sym': ['2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
                    '20']}
        self.wavelet_default = {'bior': 11, 'coif': 11, 'db': 9, 'dmey': 0,
                                'haar': 0, 'rbio': 11, 'sym': 7}
        self.wavelet_start_default = {'bior': 9, 'coif': 7, 'db': 9, 'dmey': 7,
                                'haar': 7, 'rbio': 9, 'sym': 9}

        # end __init__

    def __str__(self):  # pragma: no cover
        return self.procs_text
        # end __str__

    def parse_reg_ex(self):
        self.ph0[0] = float(self.reg_ex.ph0.findall(self.procs_text)[0][0])
        self.ph1[0] = float(self.reg_ex.ph1.findall(self.procs_text)[0][0])
        self.n_points[0] = int(self.reg_ex.n_points.findall(self.procs_text)[0])
        self.lb[0] = float(self.reg_ex.lb.findall(self.procs_text)[0][0])
        self.gb[0] = float(self.reg_ex.gb.findall(self.procs_text)[0][0])
        self.ssb[0] = float(self.reg_ex.ssb.findall(self.procs_text)[0][0])
        self.sf[0] = float(self.reg_ex.sf.findall(self.procs_text)[0][0])
        # self.strip_start     = int(self.reg_ex.stsr.findall(self.procs_text)[0][0])
        # self.strip_end       = int(self.reg_ex.stsi.findall(self.procs_text)[0][0])
        if self.ssb[0] == 2.0:
            self.ssb[0] = 90.0

        if self.procs_text.find("$AXNUC=") > -1:
            self.axis_nucleus[0] = self.reg_ex.axis_nucleus.findall(self.procs_text)[0]
            self.axis_nucleus[0] = self.axis_nucleus[0].replace('<', '')
            self.axis_nucleus[0] = self.axis_nucleus[0].replace('>', '')

        if self.proc2s_text.find("$AXNUC=") > -1:
            self.axis_nucleus[1] = self.reg_ex.axis_nucleus.findall(self.proc2s_text)[0]
            self.axis_nucleus[1] = self.axis_nucleus[1].replace('<', '')
            self.axis_nucleus[1] = self.axis_nucleus[1].replace('>', '')

        self.window_type[0] = int(self.reg_ex.wdw.findall(self.procs_text)[0])
        self.aunmp = self.reg_ex.aunmp.findall(self.procs_text)[0]
        try:
            self.sf[1] = float(self.reg_ex.sf.findall(self.proc2s_text)[0][0])
        except:
            pass

        try:
            self.ph0[1] = float(self.reg_ex.ph0.findall(self.proc2s_text)[0][0])
            self.ph1[1] = float(self.reg_ex.ph1.findall(self.proc2s_text)[0][0])
            self.n_points[1] = int(self.reg_ex.n_points.findall(self.proc2s_text)[0])
            self.lb[1] = float(self.reg_ex.lb.findall(self.proc2s_text)[0][0])
            self.gb[1] = float(self.reg_ex.gb.findall(self.proc2s_text)[0][0])
            self.ssb[1] = float(self.reg_ex.ssb.findall(self.proc2s_text)[0][0])
            if self.ssb[1] == 2.0:
                self.ssb[1] = 90.0

            #self.wdw[1] = float(self.reg_ex.wdw.findall(self.proc2s_text)[0][0])
            # self.axis_nucleus[1] = self.reg_ex.axis_nucleus.findall(self.proc2s_text)[0]
            self.window_type[1] = int(self.reg_ex.wdw.findall(self.proc2s_text)[0])

        except:
            a = 2 + 3
            # print("1D data")

        # end parse_reg_ex

    def parse_reg_ex_varian(self):
        dd = self.reg_ex_varian.td.search(self.procs_text)
        if hasattr(dd, 'span'):
            dd = self.procs_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            td = int(dd)
            self.n_points[0] = pow(2, math.ceil(math.log(td) / math.log(2)))

        ni = 1
        dd = self.reg_ex_varian.ni.search(self.procs_text)
        if hasattr(dd, 'span'):
            dd = self.procs_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            ni = int(dd)

        ni = max(ni, 1)
        ni2 = 1
        dd = self.reg_ex_varian.ni2.search(self.procs_text)
        if hasattr(dd, 'span'):
            dd = self.procs_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            ni2 = int(dd)

        ni2 = max(ni2, 1)
        phase = np.array([1], dtype='int')
        dd = self.reg_ex_varian.phase.search(self.procs_text)
        if hasattr(dd, 'span'):
            dd = self.procs_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            phase = np.array(dd.split(), dtype='int')

        phase2 = np.array([1], dtype='int')
        dd = self.reg_ex_varian.phase2.search(self.procs_text)
        if hasattr(dd, 'span'):
            dd = self.procs_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            phase2 = np.array(dd.split(), dtype='int')

        self.n_points[1] = pow(2, math.ceil(math.log(ni * ni2 * len(phase) * len(phase2)) / math.log(2)))
        dd = self.reg_ex_varian.nuc1.search(self.procs_text)
        if hasattr(dd, 'span'):
            dd = self.procs_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            dd = dd.replace(" ", "")
            dd = dd.replace("\"", "")
            self.axis_nucleus[0] = dd

        dd = self.reg_ex_varian.nuc2.search(self.procs_text)
        if hasattr(dd, 'span'):
            dd = self.procs_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            dd = dd.replace(" ", "")
            dd = dd.replace("\"", "")
            self.axis_nucleus[1] = dd

        dd = self.reg_ex_varian.nuc3.search(self.procs_text)
        if hasattr(dd, 'span'):
            dd = self.procs_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            dd = dd.replace(" ", "")
            dd = dd.replace("\"", "")
            self.axis_nucleus[2] = dd

        # end parse_reg_ex_varian

    def read(self, spc_dir):
        procs_name = spc_dir + os.sep + 'pdata' + os.sep + '1' + os.sep + 'procs'
        proc2s_name = spc_dir + os.sep + 'pdata' + os.sep + '1' + os.sep + 'proc2s'
        proc3s_name = spc_dir + os.sep + 'pdata' + os.sep + '1' + os.sep + 'proc3s'
        procpar_name = spc_dir + os.sep + 'procpar'
        procpar_name2 = spc_dir + os.sep + 'PROCPAR'
        if (os.path.isfile(procs_name)):
            try:
                f = open(procs_name, "r")
                self.procs_text = f.read()
                f.close()

            except:
                f = open(procs_name, "r", encoding='latin-1')
                self.procs_text = f.read()
                f.close()

            self.data_type = 'Bruker'

        if (os.path.isfile(proc2s_name)):
            try:
                f = open(proc2s_name, "r")
                self.proc2s_text = f.read()
                f.close()

            except:
                f = open(proc2s_name, "r", encoding='latin-1')
                self.proc2s_text = f.read()
                f.close()

            self.data_type = 'Bruker'

        if (os.path.isfile(proc3s_name)):
            try:
                f = open(proc3s_name, "r")
                self.proc3s_text = f.read()
                f.close()

            except:
                f = open(proc3s_name, "r", encoding='latin-1')
                self.proc3s_text = f.read()
                f.close()

            self.data_type = 'Bruker'

        if (os.path.isfile(procpar_name)):
            f = open(procpar_name, "r")
            self.procs_text = f.read()
            f.close()
            self.data_type = 'Varian'

        if (os.path.isfile(procpar_name2)):
            f = open(procpar_name2, "r")
            self.procs_text = f.read()
            f.close()
            self.data_type = 'Varian'

        if self.data_type == 'Bruker':
            self.parse_reg_ex()
        else:
            self.parse_reg_ex_varian()

        # end read
