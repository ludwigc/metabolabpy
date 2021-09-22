import numpy as np
# from metabolabpy.nmr import wdwf
from metabolabpy.nmr import acqPars
from metabolabpy.nmr import procPars
from metabolabpy.nmr import dispPars
import os
from scipy.fftpack import fft, ifft, fftshift
import math
import matplotlib.pyplot as pl
from scipy import signal
import time
from metabolabpy.nmr import nmrpipeData
from scipy.optimize import leastsq
from metabolabpy.nmr import apcbc
import metabolabpy.__init__ as ml_version


class NmrData:

    def __init__(self):
        self.fid = np.array([[]], dtype='complex')
        self.spc = np.array([[]], dtype='complex')
        self.ppm1 = np.array([], dtype='float64')
        self.ppm2 = np.array([], dtype='float64')
        self.ppm3 = np.array([], dtype='float64')
        self.start_peak = np.array([], dtype='float64')
        self.end_peak = np.array([], dtype='float64')
        self.start_peak_points = np.array([], dtype='int')
        self.end_peak_points = np.array([], dtype='int')
        self.peak_int = np.array([], dtype='float64')
        self.peak_max = np.array([], dtype='float64')
        self.peak_max_ppm = np.array([], dtype='float64')
        self.peak_max_points = np.array([], dtype='int')
        self.peak_label = np.array([], dtype='str')
        self.dim = 0
        self.title = np.array([], dtype='str')
        self.orig_data_set = str('')
        self.ph_corr_mode = 0
        self.acq = acqPars.AcqPars()
        self.proc = procPars.ProcPars()
        self.display = dispPars.DispPars()
        self.fid_offset_corr = 0
        self.data_set_name = ''
        self.data_set_number = ''
        self.title = str('Empty NMR data set')
        self.pulse_program = str('')
        self.window_function = {'none': 0, 'exponential': 1, 'gaussian': 2, 'sine': 3, 'qsine': 4, 'sem': 5}
        self.ref_shift = np.array([0, 0, 0], dtype='float64')
        self.ref_point = np.array([0, 0, 0], dtype='int')
        self.refsw = np.array([0, 0, 0], dtype='float64')
        self.ref_tmsp_range = 0.3  # [ppm]
        self.apc = apcbc.Apcbc()
        self.projected_j_res = False
        self.orig_j_res_set = -1
        self.orig_j_res_exp = -1
        self.pj_res_mode = ''
        self.j_res = False
        self.acqu = str('')
        self.acqu2 = str('')
        self.acqu3 = str('')
        self.audita_txt = str('')
        self.format_temp = str('')
        self.fq1list = str('')
        self.scon2 = str('')
        self.shimvalues = str('')
        self.uxnmr_par = str('')
        self.proc1 = str('')
        self.proc2 = str('')
        self.proc3 = str('')
        self.outd = str('')
        self.ver = ml_version.__version__
        # end __init__

    def __str__(self):  # pragma: no cover
        r_string = '______________________________________________________________________________________\n'
        r_string += '\nMetaboLabPy NMR Data (v. ' + self.ver + ')\n'
        r_string += '______________________________________________________________________________________\n\n'
        r_string += self.title
        return r_string
        # end __str__

    def add_peak(self, start_end=np.array([], dtype='float64'), peak_label=''):
        if len(start_end) > 0:
            start_peak = int(1e4 * self.points2ppm(self.ppm2points(max(start_end), 0), 0)) / 1e4
            end_peak = int(1e4 * self.points2ppm(self.ppm2points(min(start_end), 0), 0)) / 1e4
            start_peak_points = self.ppm2points(start_peak, 0)
            end_peak_points = self.ppm2points(end_peak, 0)
            self.start_peak = np.append(self.start_peak, start_peak)
            self.end_peak = np.append(self.end_peak, end_peak)
            self.start_peak_points = np.append(self.start_peak_points, start_peak_points)
            self.end_peak_points = np.append(self.end_peak_points, end_peak_points)
            self.peak_label = np.append(self.peak_label, peak_label)
            npts = len(self.spc[0])
            spc = self.spc[0][npts - start_peak_points:npts - end_peak_points].real
            max_pos = np.where(spc == np.max(spc))[0][0] + 1
            self.peak_max = np.append(self.peak_max, np.max(spc))
            self.peak_max_points = np.append(self.peak_max_points, start_peak_points - max_pos)
            self.peak_max_ppm = np.append(self.peak_max_ppm, self.points2ppm(start_peak_points - max_pos))
            self.peak_int = np.append(self.peak_int, sum(spc - np.mean([spc[0], spc[-1]])))
            sort_idx = np.argsort(self.start_peak)[::-1]
            self.start_peak = self.start_peak[sort_idx]
            self.end_peak = self.end_peak[sort_idx]
            self.start_peak_points = self.start_peak_points[sort_idx]
            self.end_peak_points = self.end_peak_points[sort_idx]
            self.peak_label = self.peak_label[sort_idx]
            self.peak_max = self.peak_max[sort_idx]
            self.peak_max_points = self.peak_max_points[sort_idx]
            self.peak_max_ppm = self.peak_max_ppm[sort_idx]
            self.peak_int = self.peak_int[sort_idx]

        # end add_peak

    def set_peak(self, start_peak, end_peak, peak_label):
        self.start_peak = np.array([], dtype='float64')
        self.end_peak = np.array([], dtype='float64')
        self.peak_label = np.array([], dtype='str')
        self.start_peak_points = np.array([], dtype='int')
        self.end_peak_points = np.array([], dtype='int')
        self.peak_max = np.array([], dtype='float64')
        self.peak_max_ppm = np.array([], dtype='float64')
        self.peak_max_points = np.array([], dtype='int')
        self.peak_int = np.array([], dtype='float64')
        for k in range(len(start_peak)):
            self.add_peak(np.array([start_peak[k], end_peak[k]]), peak_label[k])

        # end set_peak

    def clear_peak(self):
        self.start_peak = np.array([], dtype='float64')
        self.end_peak = np.array([], dtype='float64')
        self.start_peak_points = np.array([], dtype='int')
        self.end_peak_points = np.array([], dtype='int')
        self.peak_label = np.array([], dtype='str')

        # end add_peak

    def apodise(self, fid, dim, lb, gb, ssb, group_delay, sw_h):
        fid = np.copy(fid)
        wdwf = np.zeros(len(fid))
        if self.proc.window_type[dim] == 0:  # no window
            wdwf = np.ones(len(fid))

        if self.proc.window_type[dim] == 1:  # exponential window
            t = (np.linspace(0.0, len(fid) - 1, len(fid)) - group_delay) / sw_h
            wdwf = np.exp(-lb * t)

        if self.proc.window_type[dim] == 2:  # gaussian window
            t = (np.linspace(0.0, len(fid) - 1 - group_delay, len(fid))) / sw_h
            wdwf = np.exp(-lb * 2 * math.pi * t - 2 * math.pi * (t ** 2) / ((2 * math.pi * gb * len(fid) / sw_h)))

        if self.proc.window_type[dim] == 3:  # sine window
            if self.acq.fn_mode == 1 or self.acq.fn_mode == 2:
                npts = int(min(self.acq.n_data_points[dim], len(fid)))
            else:
                npts = int(min(self.acq.n_data_points[dim] / 2, len(fid)))

            t = (np.linspace(0.0, npts - 1, npts)) / (npts - 1)
            wdwf = np.zeros(len(fid))
            ssb2 = ssb * math.pi / 180.0
            wdwf[:npts] = np.sin(ssb2 + (math.pi - ssb2) * t)

        if self.proc.window_type[dim] == 4:  # qsine window
            if self.acq.fn_mode == 1 or self.acq.fn_mode == 2:
                npts = int(min(self.acq.n_data_points[dim], len(fid)))
            else:
                npts = int(min(self.acq.n_data_points[dim] / 2, len(fid)))

            t = (np.linspace(0.0, npts - 1, npts)) / (npts - 1)
            wdwf = np.zeros(len(fid))
            ssb2 = ssb * math.pi / 180.0
            wdwf[:npts] = np.sin(ssb2 + (math.pi - ssb2) * t) ** 2

        if self.proc.window_type[dim] == 5:  # sem window
            if self.acq.fn_mode == 1 or self.acq.fn_mode == 2:
                npts = int(min(self.acq.n_data_points[dim], len(fid)))
            else:
                npts = int(min(self.acq.n_data_points[dim] / 2, len(fid)))

            t1 = (np.linspace(0.0, npts - 1 - group_delay, npts)) / sw_h
            t2 = (np.linspace(0.0, npts - 1, npts)) / (npts - 1)
            wdwf = np.zeros(len(fid))
            ssb2 = ssb * math.pi / 180.0
            wdwf[:npts] = np.exp(-lb * t1) * np.sin(ssb2 + (math.pi - ssb2) * t2)

        fid = np.copy(wdwf * fid)
        return fid
        # end apodise

    def autobaseline1d(self):
        spc = self.spc[0]
        self.apc.npts = len(spc)
        scale_fact = np.max(np.abs(spc))
        spc /= scale_fact
        xaxis = np.linspace(-self.apc.n_max, self.apc.n_max, self.apc.npts)
        par_eval = self.apc.fit_baseline(spc, xaxis)
        spc2 = self.apc.baseline_fit_func_eval(par_eval, spc, xaxis)  # , False)
        self.apc.pars = par_eval
        self.apc.r_spc = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self.apc.i_spc = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self.apc.r_spc[:int(len(par_eval) / 2)] = par_eval[:int(len(par_eval) / 2)]
        self.apc.i_spc[:int(len(par_eval) / 2)] = par_eval[int(len(par_eval) / 2):]
        self.spc[0] = spc2 * scale_fact
        self.apc.correct_baseline = 1
        self.proc_spc1d()
        self.baseline1d()
        # end autobaseline1d

    def autophase1d(self):
        spc = self.spc[0]
        self.apc.npts = len(spc)
        scale_fact = np.max(np.abs(spc))
        spc /= scale_fact
        xaxis = np.linspace(-self.apc.n_max, self.apc.n_max, self.apc.npts)
        par_eval = self.apc.fit_phase(spc, xaxis)
        spc2 = self.apc.phase_fit_func_eval(par_eval, spc, xaxis)  # , False)
        if (np.min(spc2.real) == -np.max(np.abs(spc2.real))):
            par_eval[0] = ((par_eval[0] * self.apc.m_fact0 + 180.0) % 360) / self.apc.m_fact0

        par_eval[0] = (((par_eval[0] * self.apc.m_fact0 + 180.0) % 360) - 180.0) / self.apc.m_fact0
        spc2 = self.apc.phase_fit_func_eval(par_eval, spc, xaxis)  # , False)
        self.apc.pars = par_eval
        self.apc.r_spc = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self.apc.i_spc = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self.apc.r_spc[:int((len(par_eval) - 2) / 2)] = par_eval[2:2 + int((len(par_eval) - 2) / 2)]
        self.apc.i_spc[:int((len(par_eval) - 2) / 2)] = par_eval[2 + int((len(par_eval) - 2) / 2):]
        self.spc[0] = spc2 * scale_fact
        self.proc.ph0[0] += par_eval[0] * self.apc.m_fact0
        self.proc.ph1[0] += par_eval[1] * self.apc.m_fact1
        self.apc.correct_baseline = 1
        self.auto_ref()
        self.proc_spc1d()
        self.baseline1d()
        # end autophase1d

    def auto_ref(self, tmsp=True):
        if self.acq.o1 == 0:
            self.ref_shift[0] = 4.76
        else:
            self.ref_shift[0] = self.acq.o1 / self.acq.bf1

        self.ref_point[0] = int(len(self.spc[0]) / 2)
        if (self.dim == 2):
            self.ref_shift[1] = (self.acq.spc_frequency[1] + self.acq.spc_offset[1]) / self.acq.spc_frequency[1] - 1.0
            self.ref_point[1] = int(len(self.spc) / 2)
            if (tmsp == True):
                self.ref_point[0] = self.ppm2points(0.0, 0)
                self.ref_shift[0] = 0.0
                pts = self.ppm2points(np.array([-self.ref_tmsp_range, self.ref_tmsp_range]), 0)
                npts = len(self.spc[0])
                r = np.arange(npts - max(pts), npts - min(pts))
                spc = np.sum(self.spc, 0)
                spc = spc[r].real
                ref_p = np.where(spc == np.amax(spc))
                self.ref_point[0] -= ref_p[0][0] - int((max(pts) - min(pts)) / 2) + 1

            self.proc.ref_point[0] = self.ref_point[0]
            self.proc.ref_point[1] = self.ref_point[1]

        if (self.dim == 1):
            if (tmsp == True):
                self.ref_point[0] = self.ppm2points(0.0, 0)
                self.ref_shift[0] = 0.0
                pts = self.ppm2points(np.array([-self.ref_tmsp_range, self.ref_tmsp_range]), 0)
                npts = len(self.spc[0])
                r = np.arange(npts - max(pts), npts - min(pts))
                r = np.copy(r[np.where(r < len(self.spc[0]))])
                spc = self.spc[0][r].real
                ref_p = np.where(spc == np.amax(spc))
                self.ref_point[0] -= ref_p[0][0] - int((max(pts) - min(pts)) / 2) + 1

        self.calc_ppm()

        # end auto_ref

    def baseline1d(self):
        spc = self.spc[0]
        self.apc.npts = len(spc)
        # scale_fact         = np.max(np.abs(spc))
        # spc              /= scale_fact
        xaxis = np.linspace(-self.apc.n_max, self.apc.n_max, self.apc.npts)
        par_eval = np.copy(self.apc.r_spc)
        par_eval = np.append(par_eval, self.apc.i_spc)
        # print(par_eval)
        spc2 = self.apc.baseline_fit_func_eval(par_eval, spc, xaxis)  # , False)
        self.spc[0] = spc2  # *scale_fact
        # end baseline1d

    def calc_ppm(self):
        if self.proc.strip_start == 0:
            stsr = 1
        else:
            stsr = self.proc.strip_start

        if self.proc.strip_end > stsr:
            stsi = self.proc.strip_end
        else:
            stsi = int(len(self.spc[0]))

        if (self.display.axis_type1 == 'ppm'):
            self.ppm1 = self.points2ppm(np.linspace(self.proc.n_points[0] - 1, 0, self.proc.n_points[0]), 0)
        else:
            self.ppm1 = self.points2hz(np.linspace(self.proc.n_points[0] - 1, 0, self.proc.n_points[0]), 0)

        self.ppm1 = self.ppm1[stsr - 1:stsi]

        if (self.dim > 1):
            npts = int(len(self.spc))
            if (self.display.axis_type2 == 'ppm'):
                self.ppm2 = self.points2ppm(np.linspace(npts - 1, 0, npts), 1)
            else:
                self.ppm2 = self.points2hz(np.linspace(npts - 1, 0, npts), 1)

        # end calc_ppm

    def conv(self, fid):
        f1 = np.copy(fid)
        x = np.linspace(0, len(fid) - 1, len(fid))
        f0 = np.copy(fid)
        filt_fid = fid[list(range(round(self.acq.group_delay)))]
        fid = fid[list(range(round(self.acq.group_delay), len(fid)))]
        ws2 = int(self.proc.conv_window_size[0] / 2)
        es = self.proc.conv_extrapolation_size[0]
        wt = int(self.proc.conv_window_type[0])
        win = np.zeros(int(ws2 * 2 + 1))
        s_win = 0
        for k in range(int(ws2 * 2 + 1)):
            if (wt == 0):  # Gaussian window
                win[k] = math.exp(-4.0 * ((k - ws2) ** 2) / (ws2 ** 2))
                s_win += win[k]
            else:  # sine bell window
                win[k] = math.cos((math.pi * (k - ws2)) / (2 * ws2 + 2))
                s_win += win[k]

        fid2 = np.convolve(win, fid) / s_win
        # Extrapolation of first sw2 data points
        idx = np.linspace(2 * ws2 - 1, 0, 2 * ws2, dtype='int')
        idx2 = np.linspace(0, 2 * ws2 - 1, 2 * ws2, dtype='int')
        fid2[idx2].real = fid2[int(2 * ws2)].real + np.mean(np.diff(fid2[idx + int(2 * ws2) - 1].real)) * idx
        fid2[idx2].imag = fid2[int(2 * ws2)].imag + np.mean(np.diff(fid2[idx + int(2 * ws2) - 1].imag)) * idx
        # Extrapolation of last sw2 data points
        idx = np.linspace(0, 2 * ws2 - 1, 2 * ws2, dtype='int')
        idx2 = np.linspace(len(fid) - int(2 * ws2) - 1, len(fid) - 1, int(2 * ws2), dtype='int')
        fid2[idx2].real = fid2[len(fid) - int(2 * ws2) - 1].real - np.mean(
            np.diff(fid2[len(fid) - int(2 * ws2) - 1 - idx].real)) * idx
        fid2[idx2].imag = fid2[len(fid) - int(2 * ws2) - 1].imag - np.mean(
            np.diff(fid2[len(fid) - int(2 * ws2) - 1 - idx].imag)) * idx
        fid2 = np.delete(fid2, np.linspace(0, 2 * ws2 - 1, 2 * ws2, dtype='int'))
        fid -= fid2  # (fidRe + 1j*fidIm)
        fid = np.concatenate([filt_fid, fid])
        return fid
        # end conv

    def export_bruker_1d(self, path_name, exp_name, scale_factor=-1):
        if self.acq.manufacturer != 'Bruker':
            return

        if os.path.isdir(path_name + os.sep + exp_name) is False:
            os.makedirs(path_name + os.sep + exp_name)

        if os.path.isdir(path_name + os.sep + exp_name + os.sep + 'pdata' + os.sep + '1') is False:
            os.makedirs(path_name + os.sep + exp_name + os.sep + 'pdata' + os.sep + '1')

        fid_file1 = path_name + os.sep + exp_name + os.sep + 'fid'
        spc_file1r = path_name + os.sep + exp_name + os.sep + 'pdata' + os.sep + '1' + os.sep + '1r'
        # write 1D FID file
        fid = np.zeros(2 * len(self.fid[0]))
        fid[0::2] = self.fid[0].real
        fid[1::2] = self.fid[0].imag
        if self.acq.byte_order == 1:
            fid.astype('>i4').tofile(fid_file1)
        else:
            if self.acq.data_type == 0:
                fid.astype('<i4').tofile(fid_file1)
            else:
                if self.acq.data_type == 1:
                    fid.astype(np.float64).tofile(fid_file1)
                else:
                    fid.astype(np.float32).tofile(fid_file1)

        # write 1r file
        spc = self.spc[0].real
        if scale_factor == -1:
            scale_factor = int(2 * np.max(spc) / 2147483647)

        spc /= scale_factor
        spc.real.astype(np.int32).tofile(spc_file1r)
        # write parameter files
        # acqu file
        f_name = path_name + os.sep + exp_name + os.sep + 'acqu'
        f = open(f_name, 'w')
        f.write(self.acq.acqus_text)
        f.close()
        # acqus file
        f_name = path_name + os.sep + exp_name + os.sep + 'acqus'
        f = open(f_name, 'w')
        f.write(self.acq.acqus_text)
        f.close()
        # pulse_program file
        f_name = path_name + os.sep + exp_name + os.sep + 'pulse_program'
        f = open(f_name, 'w')
        f.write(self.pulse_program)
        f.close()
        # audita.txt file
        f_name = path_name + os.sep + exp_name + os.sep + 'audita.txt'
        f = open(f_name, 'w')
        f.write(self.audita_txt)
        f.close()
        # format.temp file
        f_name = path_name + os.sep + exp_name + os.sep + 'format.temp'
        f = open(f_name, 'w')
        f.write(self.format_temp)
        f.close()
        # fq1list file
        f_name = path_name + os.sep + exp_name + os.sep + 'fq1list'
        f = open(f_name, 'w')
        f.write(self.fq1list)
        f.close()
        # scon2 file
        f_name = path_name + os.sep + exp_name + os.sep + 'scon2'
        f = open(f_name, 'w')
        f.write(self.scon2)
        f.close()
        # shimvalues file
        f_name = path_name + os.sep + exp_name + os.sep + 'shimvalues'
        f = open(f_name, 'w')
        f.write(self.shimvalues)
        f.close()
        # uxnmr.par file
        f_name = path_name + os.sep + exp_name + os.sep + 'uxnmr.par'
        f = open(f_name, 'w')
        f.write(self.uxnmr_par)
        f.close()
        # procs file
        orig_si = self.proc.reg_ex.si.findall(self.proc.procs_text)[0]
        ft_mod = int(self.proc.reg_ex.ft_mod.findall(self.proc.procs_text)[0])
        nc_proc = int(self.proc.reg_ex.nc_proc.findall(self.proc.procs_text)[0])
        procs_text = self.proc.procs_text.replace(orig_si, str(len(self.spc[0])))
        procs_text = procs_text.replace('FT_mod= ' + str(ft_mod), 'FT_mod= 6')
        procs_text = procs_text.replace('NC_proc= ' + str(nc_proc), 'NC_proc= 0')
        f_name = path_name + os.sep + exp_name + os.sep + 'pdata' + os.sep + '1' + os.sep + 'procs'
        f = open(f_name, 'w')
        f.write(procs_text)
        f.close()
        # proc file
        f_name = path_name + os.sep + exp_name + os.sep + 'pdata' + os.sep + '1' + os.sep + 'proc'
        f = open(f_name, 'w')
        f.write(procs_text)
        f.close()
        # outd file
        f_name = path_name + os.sep + exp_name + os.sep + 'pdata' + os.sep + '1' + os.sep + 'outd'
        f = open(f_name, 'w')
        f.write(self.outd)
        f.close()
        # title file
        f_name = path_name + os.sep + exp_name + os.sep + 'pdata' + os.sep + '1' + os.sep + 'title'
        f = open(f_name, 'w')
        f.write(self.title)
        f.close()
        # end export_bruker_1d

    def fid_offset_correction(self, fid):
        fid = np.copy(fid)
        if (self.fid_offset_corr > 0):
            m_mean = np.mean(fid[:self.fid_offset_corr])
            fid -= m_mean

        return fid
        # end fid_offset_correction

    def gibbs(self, fid):
        if (self.proc.gibbs[0] == True):
            fid[0] /= 2.0

        return fid
        # end gibbs

    def hilbert(self, mat, dim=0):
        if (dim == 1):
            mat = np.ndarray.transpose(mat)

        npts = len(mat[0])
        npts1 = len(mat)
        v1 = np.ones(npts1)
        mat1 = np.array([[]], dtype='complex')
        mat1 = np.resize(mat1, (npts1, npts))
        b_mat = np.zeros(int(2 * npts), dtype='complex')
        b_mat[:(npts + 1)] = np.ones(npts + 1)
        b_mat[1:npts] += b_mat[1:npts]
        z_mat = np.zeros(int(2 * npts), dtype='complex')
        b_mat = np.outer(v1, b_mat)
        z_mat = np.outer(v1, z_mat)
        z_mat[:, :npts] = mat
        z_mat = np.fft.ifft(b_mat * np.fft.fft(z_mat))
        mat = z_mat[:, :npts]
        if (dim == 1):
            mat = np.ndarray.transpose(mat)

        return mat
        # end hilbert

    def hilbert1(self, mat, dim):
        if (self.dim == 1):
            npts = len(mat)
            b_mat = np.zeros(int(2 * npts), dtype='complex')
            z_mat = np.zeros(int(2 * npts), dtype='complex')
            z_mat[:int(len(mat))] = mat
            b_mat[:(npts + 1)] = np.ones(npts + 1)
            b_mat[1:npts] += b_mat[1:npts]
            mat2 = np.zeros(2 * npts, dtype='complex')
            mat2[:len(mat)] = mat
            mat2 = ifft(b_mat * fft(mat2))
            mat1 = np.copy(mat2[:len(mat)])

        if (self.dim == 2):
            if (dim == 1):
                mat = np.ndarray.transpose(mat)

            npts = len(mat[0])
            npts1 = len(mat)
            mat1 = np.array([[]], dtype='complex')
            mat1 = np.resize(mat1, (npts1, npts))
            for k in range(len(mat)):
                b_mat = np.zeros(int(2 * npts), dtype='complex')
                z_mat = np.zeros(int(2 * npts), dtype='complex')
                z_mat[:int(len(mat[k]))] = mat[k]
                b_mat[:(npts + 1)] = np.ones(npts + 1)
                b_mat[1:npts] += b_mat[1:npts]
                mat2 = np.zeros(2 * npts, dtype='complex')
                mat2[:len(mat[k])] = mat[k]
                mat2 = ifft(b_mat * fft(mat2))
                mat1[k] = mat2[:len(mat[k])]

            if (dim == 1):
                mat1 = np.ndarray.transpose(mat1)

        return mat1
        # end hilbert1

    def multiply(self, factor=1.0):
        self.spc = factor * self.spc

    def phase(self, ph0, ph1, npts):
        ph0 = -ph0 * math.pi / 180.0
        ph1 = -ph1 * math.pi / 180.0
        t = complex()
        frac = np.linspace(0, 1, npts)
        ph = ph0 + frac * ph1
        self.spc[0] = np.cos(ph) * self.spc[0].real + np.sin(ph) * self.spc[0].imag + 1j * (
                -np.sin(ph) * self.spc[0].real + np.cos(ph) * self.spc[0].imag)
        # end phase

    def phase2(self, mat, ph0, ph1):
        npts = len(mat)
        ph0 = -ph0 * math.pi / 180.0
        ph1 = -ph1 * math.pi / 180.0
        t = complex()
        frac = np.linspace(0, 1, npts)
        ph = ph0 + frac * ph1
        mat = np.cos(ph) * mat.real + np.sin(ph) * mat.imag + 1j * (-np.sin(ph) * mat.real + np.cos(ph) * mat.imag)
        return mat
        # end phase2

    def phase2a(self, ph0, ph1, dim=0):
        mat = np.copy(self.spc.real)
        mat = self.hilbert(mat, dim)
        if dim == 1:
            mat = np.ndarray.transpose(mat)

        npts = len(mat[0])
        npts1 = len(mat)
        v1 = np.ones(npts1)
        ph0 = -ph0 * math.pi / 180.0
        ph1 = -ph1 * math.pi / 180.0
        frac = np.outer(v1, np.linspace(0, 1, npts))
        ph = ph0 + frac * ph1
        mat = np.cos(ph) * mat.real + np.sin(ph) * mat.imag + 1j * (-np.sin(ph) * mat.real + np.cos(ph) * mat.imag)
        if dim == 1:
            mat = np.ndarray.transpose(mat)

        self.spc = np.copy(mat.real)
        # end phase2a

    def phase2d(self, ph0, ph1, dim):
        mat = self.hilbert1(self.spc, dim)
        if dim == 1:
            mat = np.ndarray.transpose(mat)

        npts = len(mat)
        for k in range(npts):
            mat[k] = self.phase2(mat[k], ph0, ph1)

        if dim == 1:
            mat = np.ndarray.transpose(mat)

        self.spc = np.copy(mat.real)
        # end phase2d

    def phase3(self, mat, ph0, ph1):
        npts = len(mat)
        ph0 = -ph0 * math.pi / 180.0
        ph1 = -ph1 * math.pi / 180.0
        t = complex()
        for k in range(int(npts)):
            frac = float(k) / float(npts)
            ph = ph0 + frac * ph1
            t = complex(math.cos(ph) * mat[k].real + math.sin(ph) * mat[k].imag,
                        -math.sin(ph) * mat[k].real + math.cos(ph) * mat[k].imag)
            mat[k] = t

        return mat
        # end phase3

    def points2hz(self, points, dim=0):
        sw = self.acq.sw_h[dim]
        if dim == 0:
            npts = int(len(self.spc[0]))

        if dim == 1:
            npts = int(len(self.spc))

        hz = sw * (points - npts / 2) / npts
        return hz
        # end points2hz

    def points2ppm(self, points, dim=0):
        sw = self.acq.sw_h[dim]
        if dim == 0:
            sfo = self.acq.sfo1
            npts = self.proc.n_points[0]  # int(len(self.spc[0]))

        if dim == 1:
            if self.acq.manufacturer == 'Bruker':
                sfo = self.proc.sf[1]
                if sfo == 0.0:
                    sfo = self.acq.sfo2

            else:
                sfo = self.acq.sfo2

            npts = int(len(self.spc))

        ppm = (sw / sfo) * (points / (npts - 1) - self.ref_point[dim] / (npts - 1)) + self.ref_shift[dim]
        return ppm
        # end points2ppm

    def ppm2points(self, ppm, dim=0):
        sw = self.acq.sw_h[dim]
        if dim == 0:
            sfo = self.acq.sfo1
            npts = int(len(self.spc[0]))

        if dim == 1:
            if self.acq.manufacturer == 'Bruker':
                sfo = self.proc.sf[1]
            else:
                sfo = self.acq.sfo2

            npts = int(len(self.spc))

        points = np.round(((ppm - self.ref_shift[dim]) * (sfo / sw) + self.ref_point[dim] / (npts - 1)) * (npts - 1))
        return points.astype(int)
        # end ppm2points

    def ppm2points2d(self, ppm):
        points = np.copy(ppm)
        for k in range(len(ppm)):
            points[k][0] = self.ppm2points(ppm[k][0], 0)
            points[k][1] = self.ppm2points(ppm[k][1], 1)

        return points
        # end ppm2points2d

    def proc_spc(self):
        if self.dim == 1:
            self.proc_spc1d()

        if self.dim == 2:
            self.proc_spc2d()

        # self.auto_ref()
        self.calc_ppm()
        # end proc_spc

    def proc_spc1d(self):
        fid = np.copy(self.fid[0])
        fid = self.water_supp(fid)
        fid = self.fid_offset_correction(fid)
        fid = self.gibbs(fid)
        fid = self.apodise(fid, 0, self.proc.lb[0], self.proc.gb[0], self.proc.ssb[0], self.acq.group_delay,
                           self.acq.sw_h[0])
        fid = self.zero_fill(fid)
        spc = fftshift(fft(fid))
        self.spc = np.resize(self.spc, (1, len(spc)))
        self.spc[0] = np.copy(spc)
        self.phase(0, 360.0 * self.acq.group_delay, self.proc.n_points[0])
        if self.proc.invert_matrix[0]:
            self.spc[0] = np.copy(np.flip(self.spc[0]))

        self.phase(self.proc.ph0[0], self.proc.ph1[0], self.proc.n_points[0])
        # end proc_spc1d

    def proc_spc2d(self, test_quad_2d=False, no_abs=False):
        fid = np.copy(self.fid)
        if self.proc.mult_factor[0] == 0:
            self.proc.mult_factor[0] = self.proc.n_points[0] / len(self.fid[0])

        if self.proc.mult_factor[1] == 0:
            self.proc.mult_factor[1] = self.proc.n_points[1] / len(self.fid)

        npts2 = len(self.spc)
        npts1 = len(self.spc[0])
        if npts1 > 0:
            if self.proc.strip_start < 1:
                stsr = 1
            else:
                stsr = self.proc.strip_start

            if self.proc.n_points[0] != npts1:
                self.ref_point[0] = int((self.proc.ref_point[0]) * self.proc.n_points[0] / (
                        self.proc.mult_factor[0] * len(self.fid[0])))  # - stsr + 1

        if npts2 > 1:
            if self.proc.n_points[1] != npts2:
                self.ref_point[1] = int(
                    self.proc.ref_point[1] * self.proc.n_points[1] / (self.proc.mult_factor[1] * len(self.fid)))

        self.spc = np.copy(np.array([[]], dtype='complex'))
        self.spc = np.copy(np.resize(self.spc, (self.proc.n_points[0], self.proc.n_points[1])))
        self.spc *= 0
        if self.proc.n_points[0] > len(fid[0]):
            fid = np.resize(fid, (self.proc.n_points[1], self.proc.n_points[0]))
            fid *= 0
            for k in range(len(self.fid)):
                fid[k][:len(self.fid[k])] = np.copy(self.fid[k][:])

        if self.proc.n_points[0] < len(fid[0]):
            fid = np.resize(fid, (len(fid), self.proc.n_points[0]))
            for k in range(len(fid)):
                fid[k] = np.copy(self.fid[k][:self.proc.n_points[0]])

        for k in range(len(self.fid)):
            fid2 = np.copy(fid[k])
            fid2 = self.water_supp(fid2)
            fid2 = self.fid_offset_correction(fid2)
            fid2 = self.gibbs(fid2)
            fid2 = self.apodise(fid2, 0, self.proc.lb[0], self.proc.gb[0], self.proc.ssb[0], self.acq.group_delay,
                                self.acq.sw_h[0])
            fid2 = self.zero_fill(fid2, 0)
            fid2 = fftshift(fft(fid2))
            if self.acq.fn_mode != 1:
                fid2 = self.phase2(fid2, 0, 360.0 * self.acq.group_delay)

            if self.proc.invert_matrix[0]:
                fid2 = np.flip(fid2)

            if self.acq.fn_mode != 1:
                fid2 = self.phase2(fid2, self.proc.ph0[0], self.proc.ph1[0])

            fid[k] = np.copy(fid2)

        fid = np.copy(np.ndarray.transpose(np.conj(fid)))
        if not test_quad_2d:
            fid = self.quad_2d(fid)
            for k in range(len(fid)):
                fid2 = np.copy(fid[k])
                fid2 = self.gibbs(fid2)
                fid2 = self.apodise(fid2, 1, self.proc.lb[1], self.proc.gb[1], self.proc.ssb[1], 0.0, self.acq.sw_h[1])
                fid2 = self.zero_fill(fid2, 1)
                fid2 = fftshift(fft(fid2))
                if self.proc.invert_matrix[1]:
                    fid2 = np.flip(fid2)

                if self.acq.fn_mode != 1:
                    fid2 = self.phase2(fid2, self.proc.ph0[1], self.proc.ph1[1])
                    self.spc[k] = fid2.real
                else:
                    self.spc[k] = np.copy(fid2)

            self.spc = np.copy(np.ndarray.transpose(self.spc))
            if (self.acq.fn_mode == 1) and (self.proc.tilt == True):
                # self.spc = np.copy(self.hilbert(self.spc, 0))
                self.tiltj_res()

            if self.acq.fn_mode == 1 and no_abs is False:
                for k in range(len(self.spc)):
                    self.spc[k] = np.copy(np.abs(self.spc[k]))

                if self.proc.symj == True:
                    self.symj_res()

        else:
            self.spc = np.copy(self.fid)

        if self.proc.strip_start == 0:
            stsr = 1
        else:
            stsr = self.proc.strip_start

        if self.proc.strip_end > stsr:
            stsi = self.proc.strip_end
            self.spc = np.ndarray.transpose(self.spc)
            self.spc = self.spc[stsr - 1:stsi]
            self.ppm1 = self.ppm1[stsr - 1:stsi]
            self.spc = np.ndarray.transpose(self.spc)

        self.spc = np.copy(self.spc.real)
        # end proc_spc2d

    def quad_2d(self, fid):
        fid = np.copy(fid)
        inc = self.acq.fn_mode
        # MetLab: 6,     0,    3,      2,            1,   4
        # FnMode: 1,     2,    3,      4,            5,   6
        #         j_res, QF, TPPI, States, States-TPPI , E/A
        rFid = np.copy(fid)
        rFid = np.resize(rFid, (len(fid), int(len(fid[0]) / 2)))
        iFid = np.copy(fid)
        iFid = np.resize(iFid, (len(fid), int(len(fid[0]) / 2)))
        if inc == 1 or inc == 2:
            fid = fid  # rFid + 1j * iFid

        if inc == 5:
            for k in range(len(fid)):
                npts = int(len(fid[k]) / 2)
                rFid[k] = fid[k][0::2].real * (-2 * (np.linspace(0, npts - 1, npts, dtype='int') % 2) + 1)
                iFid[k] = fid[k][1::2].real * (2 * (np.linspace(0, npts - 1, npts, dtype='int') % 2) - 1)

            fid = np.copy(rFid + 1j * iFid)

        if inc == 6:
            for k in range(len(fid)):
                rFid[k] = np.conj(fid[k][0::2] + fid[k][1::2])
                iFid[k] = np.conj(fid[k][0::2] - fid[k][1::2])
                iFid[k] = iFid[k].imag - 1j * iFid[k].real
                rFid[k] = rFid[k].imag + 1j * iFid[k].imag

            fid = np.copy(rFid)

        return fid
        # end quad_2d

    def read_pipe_2d(self, p_name, f_name):
        npd = nmrpipeData.NmrPipeData()
        npd.read_pipe(p_name, f_name)
        self.spc = npd.spc
        self.proc.sw_h[0] = npd.fdf2sw
        self.ref_shift[0] = npd.fdf2orig / self.acq.sfo1
        self.ref_point[0] = 0
        self.proc.sw_h[1] = npd.fdf1sw
        self.ref_shift[1] = npd.fdf1orig / self.acq.sfo2
        self.ref_point[1] = 0
        self.proc.phase_inversion = False
        self.proc.mult_factor = [1, 1]
        self.proc.n_points[0] = len(self.spc[0])
        self.proc.n_points[1] = len(self.spc)
        self.calc_ppm()
        # end read_pipe_2d

    def read_spc(self):
        self.acq.read(self.data_set_name + os.sep + self.data_set_number)
        self.proc.read(self.data_set_name + os.sep + self.data_set_number)
        # self.proc.sw_h = self.acq.sw_h
        self.orig_data_set = self.data_set_name + os.sep + self.data_set_number
        if self.acq.manufacturer == 'Bruker':
            title_file = self.data_set_name + os.sep + self.data_set_number + os.sep + 'pdata' + os.sep + '1' + os.sep + 'title'
            if (os.path.isfile(title_file)):
                fid = open(title_file, "r")
                self.title = fid.read()
                fid.close()

            acqu_file = self.data_set_name + os.sep + self.data_set_number + os.sep + 'acqu'
            if (os.path.isfile(acqu_file)):
                fid = open(acqu_file, "r")
                self.acqu = fid.read()
                fid.close()

            acqu2_file = self.data_set_name + os.sep + self.data_set_number + os.sep + 'acqu2'
            if (os.path.isfile(acqu2_file)):
                fid = open(acqu2_file, "r")
                self.acqu2 = fid.read()
                fid.close()

            acqu3_file = self.data_set_name + os.sep + self.data_set_number + os.sep + 'acqu3'
            if (os.path.isfile(acqu3_file)):
                fid = open(acqu3_file, "r")
                self.acqu3 = fid.read()
                fid.close()

            audita_file = self.data_set_name + os.sep + self.data_set_number + os.sep + 'audita.txt'
            if (os.path.isfile(audita_file)):
                fid = open(audita_file, "r")
                self.audita_txt = fid.read()
                fid.close()

            format_file = self.data_set_name + os.sep + self.data_set_number + os.sep + 'format.temp'
            if (os.path.isfile(format_file)):
                fid = open(format_file, "r")
                self.format_temp = fid.read()
                fid.close()

            fq1list_file = self.data_set_name + os.sep + self.data_set_number + os.sep + 'fq1list'
            if (os.path.isfile(fq1list_file)):
                fid = open(fq1list_file, "r")
                self.fq1list = fid.read()
                fid.close()

            scon2_file = self.data_set_name + os.sep + self.data_set_number + os.sep + 'scon2'
            if (os.path.isfile(scon2_file)):
                fid = open(scon2_file, "r")
                self.scon2 = fid.read()
                fid.close()

            shimvalues_file = self.data_set_name + os.sep + self.data_set_number + os.sep + 'shimvalues'
            if (os.path.isfile(shimvalues_file)):
                fid = open(shimvalues_file, "r")
                self.shimvalues = fid.read()
                fid.close()

            uxnmrpar_file = self.data_set_name + os.sep + self.data_set_number + os.sep + 'uxnmr.par'
            if (os.path.isfile(uxnmrpar_file)):
                fid = open(uxnmrpar_file, "r")
                self.uxnmr_par = fid.read()
                fid.close()

            proc_file = self.data_set_name + os.sep + self.data_set_number + os.sep + 'pdata' + os.sep + '1' + os.sep + 'proc'
            if (os.path.isfile(proc_file)):
                fid = open(proc_file, "r")
                self.proc1 = fid.read()
                fid.close()

            proc2_file = self.data_set_name + os.sep + self.data_set_number + os.sep + 'pdata' + os.sep + '1' + os.sep + 'proc2'
            if (os.path.isfile(proc2_file)):
                fid = open(proc2_file, "r")
                self.proc2 = fid.read()
                fid.close()

            proc3_file = self.data_set_name + os.sep + self.data_set_number + os.sep + 'pdata' + os.sep + '1' + os.sep + 'proc3'
            if (os.path.isfile(proc3_file)):
                fid = open(proc3_file, "r")
                self.proc3 = fid.read()
                fid.close()

            outd_file = self.data_set_name + os.sep + self.data_set_number + os.sep + 'pdata' + os.sep + '1' + os.sep + 'outd'
            if (os.path.isfile(outd_file)):
                fid = open(outd_file, "r")
                self.outd = fid.read()
                fid.close()

            self.acq.sfo2 = self.proc.sf[1]
            self.display.y_label = self.proc.axis_nucleus[1]

        else:
            title_file1 = self.data_set_name + os.sep + self.data_set_number + os.sep + 'text'
            title_file2 = self.data_set_name + os.sep + self.data_set_number + os.sep + 'TEXT'
            if os.path.isfile(title_file1):
                fid = open(title_file1, "r")
                self.title = fid.read()
                fid.close()

            if os.path.isfile(title_file2):
                fid = open(title_file2, "r")
                self.title = fid.read()
                fid.close()

        if self.acq.manufacturer == 'Bruker':
            pul_prog_file = self.data_set_name + os.sep + self.data_set_number + os.sep + 'pulseprogram'
        else:
            pul_prog_file = self.data_set_name + os.sep + self.data_set_number + os.sep + self.acq.pul_prog_name + ".c"

        if os.path.isfile(pul_prog_file):
            fid = open(pul_prog_file, "r")
            self.pulse_program = fid.read()
            fid.close()

        if self.acq.manufacturer == 'Bruker':
            fid_file1 = self.data_set_name + os.sep + self.data_set_number + os.sep + 'fid'
            spc_file1r = self.data_set_name + os.sep + self.data_set_number + os.sep + 'pdata' + os.sep + '1' + os.sep + '1r'
            spc_file1i = self.data_set_name + os.sep + self.data_set_number + os.sep + 'pdata' + os.sep + '1' + os.sep + '1i'
            fid_file2 = self.data_set_name + os.sep + self.data_set_number + os.sep + 'ser'
            if os.path.isfile(fid_file1):
                # read 1D FID file
                self.fid = np.resize(self.fid, (1, int(self.acq.n_data_points[0] / 2)))
                f = open(fid_file1, 'rb')
                if self.acq.byte_order == 1:
                    fid = np.fromfile(f, dtype='>i4')
                else:
                    if self.acq.data_type == 0:
                        fid = np.fromfile(f, dtype='<i4')
                    else:
                        if self.acq.data_type == 1:
                            fid = np.fromfile(f, dtype=np.float64)
                        else:
                            if self.acq.data_type == 2:
                                fid = np.fromfile(f, dtype=np.double)
                            else:
                                fid = np.fromfile(f, dtype=np.float32)

                f.close()
                self.fid[0].real = fid[0::2]
                self.fid[0].imag = -fid[1::2]
                self.dim = 1

            if os.path.isfile(spc_file1r):
                # read 1D spectrum file (real part)
                f = open(spc_file1r, 'rb')
                fid = np.fromfile(f, dtype=np.int32)
                self.spc = np.resize(self.spc, (1, int(len(fid))))
                self.spc[0].real = fid
                f.close()
                if os.path.isfile(spc_file1i):
                    # read 1D spectrum file (imaginary part)
                    f = open(spc_file1i, 'rb')
                    fid = np.fromfile(f, dtype=np.int32)
                    f.close()
                    self.spc[0].imag = fid


            elif (os.path.isfile(fid_file2)):
                # read 2D spectrum
                np1 = int(self.acq.n_data_points[0])
                np2 = int(self.acq.n_data_points[1])
                self.fid = np.resize(self.fid, (int(np2), int(np1 / 2)))
                f = open(fid_file2, 'rb')
                fid = np.fromfile(f, dtype=np.int32)
                f.close()
                fid = fid.reshape(int(np2), int(np1))
                for x in range(np2):
                    self.fid[x].real = fid[x][0::2]
                    self.fid[x].imag = -fid[x][1::2]

                self.dim = 2
                if self.acq.pul_prog_name.find("j_res"):
                    self.j_res = True
                    self.proc.tilt = True
                    self.proc.symj = True

            self.proc.sw_h = np.copy(self.acq.sw_h)
        elif self.acq.manufacturer == 'Varian':
            fid_file1 = self.data_set_name + os.sep + self.data_set_number + os.sep + 'fid'
            fid_file2 = self.data_set_name + os.sep + self.data_set_number + os.sep + 'FID'
            fid_file = ''
            if os.path.isfile(fid_file1):
                fid_file = fid_file1
            elif os.path.isfile(fid_file2):
                fid_file = fid_file2

            if len(fid_file) > 0:
                self.dim = self.acq.np + self.acq.np2 - 1
                if self.dim == 1 and self.acq.ni > 1:
                    self.dim = 2
                    self.j_res = True
                    self.proc.tilt = True
                    self.proc.symj = True
                    self.acq.fn_mode = 1

                f = open(fid_file, 'rb')
                dh1 = np.fromfile(f, dtype='>i4', count=6)
                dh2 = np.fromfile(f, dtype='>i2', count=2)
                dh3 = np.fromfile(f, dtype='>i4', count=1)
                status = np.binary_repr(dh2[1], 16)
                data_format = '>i2'
                if status[12] == '1':
                    data_format = '>f'
                elif status[13] == '1':
                    data_format = '>i4'

                n_blocks = dh1[0]
                n_traces = dh1[1]
                n_points = dh1[2]
                if self.dim == 1:
                    self.fid = np.resize(self.fid, (1, int(self.acq.n_data_points[0] / 2)))
                    fid = np.array([])
                    for k in range(n_blocks):
                        bh1 = np.fromfile(f, dtype='>i2', count=4)
                        bh2 = np.fromfile(f, dtype='>i4', count=1)
                        bh3 = np.fromfile(f, dtype='>f', count=4)
                        if status[10] == '1':
                            bh4 = np.fromfile(f, dtype='>i2', count=4)
                            bh5 = np.fromfile(f, dtype='>i4', count=1)
                            bh6 = np.fromfile(f, dtype='>f', count=4)

                        data = np.fromfile(f, dtype=data_format, count=n_traces * n_points)
                        fid = np.append(fid, data[0::2] + 1j * data[1::2])

                    self.fid[0] = fid

                else:
                    ni = int(self.acq.n_data_points[1] / 2)
                    td = int(self.acq.n_data_points[0] / 2)
                    self.fid = np.zeros((ni, td), dtype='complex')
                    for k in range(ni):
                        bh1 = np.fromfile(f, dtype='>i2', count=4)
                        bh2 = np.fromfile(f, dtype='>i4', count=1)
                        bh3 = np.fromfile(f, dtype='>f', count=4)
                        status = np.binary_repr(bh1[1], 8)
                        data_format = '>i2'
                        if status[4] == '1':
                            data_format = '>f'
                        elif status[5] == '1':
                            data_format = '>i4'

                        data = np.fromfile(f, dtype=data_format, count=n_traces * n_points)
                        self.fid[k] = data[0::2] + 1j * data[1::2]

                f.close()

        # end read_spc

    def set_ref(self, ref_shift, ref_point):
        for k in range(len(ref_shift)):
            self.ref_shift[k] = ref_shift[k]

        for k in range(len(ref_point)):
            self.ref_point[k] = ref_point[k]

        self.calc_ppm()
        # end set_ref    

    def set_window_function(self, dim, wdwf):
        try:
            self.proc.window_type[dim] = self.window_function[wdwf]

        except:
            self.proc.window_type[dim] = self.window_function['none']
            print('Unknown windowType, setting window_function to "none"')

        # end set_window_function

    def smo(self, fid):
        x = np.linspace(0, len(fid) - 1, len(fid))
        f0 = np.copy(fid)
        fid = np.roll(fid, math.floor(-self.acq.group_delay))
        pp = np.polyfit(x, fid, self.proc.poly_order)
        fid = fid - np.polyval(pp, x)
        fid = np.roll(fid, math.ceil(self.acq.group_delay))
        return fid
        # end smo

    def symj_res(self):
        for k in range(int(len(self.spc) / 2)):
            tmp = np.minimum(self.spc[k], self.spc[len(self.spc) - k - 1])
            self.spc[k] = tmp
            self.spc[len(self.spc) - k - 1] = tmp

        # end symj

    def tiltj_res(self):
        hz_per_point = self.acq.sw_h[0] / len(self.spc[0])
        sw2 = self.acq.sw_h[1]
        npts2 = len(self.spc)
        hz_vect = np.linspace(sw2 / 2.0, -sw2 / 2.0, npts2)
        # hz_vect = hz_vect - (hz_vect[1] - hz_vect[0]) / 2.0
        for k in range(npts2):
            npts = hz_vect[k] / hz_per_point
            fid1 = np.copy(ifft(self.spc[k]))
            fid1 = self.phase2(fid1, 0, -npts * 360.0)

            self.spc[k] = np.copy(fft(fid1))

        # end tiltj_res

    def water_supp(self, fid):
        if (self.proc.water_suppression == 1):
            fid = self.conv(fid)

        if (self.proc.water_suppression == 2):
            fid = self.smo(fid)

        return fid
        # end water_supp

    def zero_fill(self, fid, dim=0):
        fid1 = np.zeros(self.proc.n_points[dim], dtype='complex')
        fid1[:int(len(fid))] = fid
        return fid1
        # end zero_fill
