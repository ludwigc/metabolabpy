#!/usr/bin/env python


"""

test_nmrDataSet

author(s): C. Ludwig
origin: 26-07-2019


"""

import unittest
import metabolabpy.nmr.nmrDataSet as nmrDataSet
import os
import numpy as np
import shutil
import darkdetect
import pandas as pd


class NrDataSetTestCase(unittest.TestCase):

    def test_ft(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.ft()
        self.assertEqual(len(nd.nmrdat[0][0].spc[0]),
                         65536)  # check number of data points in Fourier transformed spectrum

    def test_ft_all(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.ft_all()
        self.assertEqual(len(nd.nmrdat[0][0].spc[0]),
                         65536)  # check number of data points in Fourier transformed spectrum

    def test_load(self):
        nd = nmrDataSet.NmrDataSet()
        data_set_path = os.path.join(os.path.dirname(__file__), "data", "loadData.mlpy")
        nd.load(data_set_path)
        self.assertEqual(nd.nmrdat[0][0].ref_shift[0], 0.0)
        self.assertEqual(nd.nmrdat[0][0].ref_point[0], 10342)

    def test_set_zero_fill(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.set_zero_fill([131072])
        nd.ft()
        self.assertEqual(len(nd.nmrdat[0][0].spc[0]),
                         131072)  # check number of data points in Fourier transformed spectrum

    def test_autobaseline1d(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.autobaseline1d()

    def test_autobaseline1d_all(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.autobaseline1d_all()

    def test_autophase1d(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.nmrdat[0][0].exclude_water = True
        nd.nmrdat[0][0].autophase1d1()
        nd.nmrdat[0][0].exclude_water = True
        nd.autophase1d()
        nd.nmrdat[0][0].set_peak(np.array([0.01]), np.array([-0.01]), np.array(['TMSP']))
        self.assertAlmostEqual(nd.nmrdat[0][0].peak_max_ppm[0] / 10.0, 0.0, places=1)

    def test_autophase1d_all(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.autophase1d_all()

    def test_auto_ref(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.auto_ref()
        self.assertEqual(nd.nmrdat[0][0].ref_point[0], 10342)
        self.assertEqual(nd.nmrdat[0][0].ref_shift[0], 0.0)
        nd = [[]]
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "3"  # 2D HSQC NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.ft()
        nd.auto_ref()
        self.assertEqual(nd.nmrdat[0][0].ref_point[0], 254)
        self.assertEqual(nd.nmrdat[0][0].ref_shift[0], 0.0)
        self.assertEqual(nd.nmrdat[0][0].ref_point[1], 2048)
        self.assertAlmostEqual(nd.nmrdat[0][0].ref_shift[1], 80.0, 4)

    def test_clear(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.clear()
        self.assertEqual(len(nd.nmrdat), 1)
        self.assertEqual(len(nd.nmrdat[0]), 0)

    def test_exclude_region(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.pp.exclude_start = [8.9]
        nd.pp.exclude_end = [11.0]
        nd.pp.flag_exclude_region = True
        nd.data_pre_processing()
        pts = len(nd.nmrdat[0][0].spc[0]) - nd.nmrdat[0][0].ppm2points([nd.pp.exclude_start, nd.pp.exclude_end], 0)
        pts2 = range(int(min(pts)[0]), int(max(pts)[0]))
        ssum = np.sum(nd.nmrdat[0][0].spc[0][pts2].real)
        self.assertEqual(ssum, 0.0)
        nd = [[]]
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "2"  # 2D Jres NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.ft()
        nd.auto_ref()
        nd.pjres(2, 'skyline')
        nd.pp.exclude_start = [8.9]
        nd.pp.exclude_end = [10.5]
        nd.pp.flag_exclude_region = True
        nd.s = 1
        nd.e = 0
        nd.data_pre_processing()
        pts = len(nd.nmrdat[1][0].spc[0]) - nd.nmrdat[1][0].ppm2points([nd.pp.exclude_start, nd.pp.exclude_end], 0)
        pts2 = range(int(min(pts)), int(max(pts)))
        ssum = np.sum(nd.nmrdat[1][0].spc[0][pts2].real)
        self.assertEqual(ssum, 0.0)

    def test_export_data_set(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.pp.exclude_start = [8.9]
        nd.pp.exclude_end = [11.0]
        nd.pp.flag_exclude_region = True
        nd.pp.bucket_points = 64
        nd.pp.flag_bucket_spectra = True
        nd.data_pre_processing()
        nd.pp.export_excel_path = os.path.expanduser("~")
        nd.pp.export_path_name = os.path.expanduser("~")
        nd.pp.export_metabo_analyst_path = os.path.expanduser("~")
        nd.pp.export_r_dolphin_path = os.path.expanduser("~")
        nd.pp.export_batman_path = os.path.expanduser("~")
        nd.pp.export_bruker_path = os.path.expanduser("~")
        nd.pp.export_file_name = "testExport.csv"
        nd.pp.export_samples_in_rows_cols = 0
        nd.pp.class_select = ["1"]
        nd.pp.export_method = 1
        nd.export_data_set()
        self.assertEqual(os.path.isfile(nd.pp.export_path_name + os.sep + nd.pp.export_file_name), True)
        os.remove(nd.pp.export_path_name + os.sep + nd.pp.export_file_name)
        nd.pp.export_file_name = "testExport.csv"
        nd.pp.export_samples_in_rows_cols = 1
        nd.pp.export_method = 1
        nd.export_data_set()
        self.assertEqual(os.path.isfile(nd.pp.export_path_name + os.sep + nd.pp.export_file_name), True)
        os.remove(nd.pp.export_path_name + os.sep + nd.pp.export_file_name)
        nd.pp.export_excel = "testExport.xlsx"
        nd.pp.export_samples_in_rows_cols = 0
        nd.pp.export_method = 0
        nd.export_data_set('init')
        nd.export_data_set('variable_meta')
        nd.export_data_set('finish')
        self.assertEqual(os.path.isfile(nd.pp.export_excel_path + os.sep + nd.pp.export_excel), True)
        nd.pp.export_samples_in_rows_cols = 1
        nd.pp.export_method = 0
        nd.export_data_set('init')
        nd.export_data_set('variable_meta')
        nd.export_data_set('finish')
        self.assertEqual(os.path.isfile(nd.pp.export_excel_path + os.sep + nd.pp.export_excel), True)
        os.remove(nd.pp.export_excel_path + os.sep + nd.pp.export_excel)
        nd.pp.export_excel = "testExport.xlsx"
        nd.pp.export_samples_in_rows_cols = 1
        nd.pp.export_method = 0
        nd.export_data_set('init')
        nd.export_data_set('finish')
        self.assertEqual(os.path.isfile(nd.pp.export_excel_path + os.sep + nd.pp.export_excel), True)
        os.remove(nd.pp.export_excel_path + os.sep + nd.pp.export_excel)
        nd.pp.export_metabo_analyst = "metaboAnalystExport.csv"
        nd.pp.export_samples_in_rows_cols = 0
        nd.pp.export_method = 2
        nd.export_data_set()
        self.assertEqual(os.path.isfile(nd.pp.export_metabo_analyst_path + os.sep + nd.pp.export_metabo_analyst), True)
        os.remove(nd.pp.export_metabo_analyst_path + os.sep + nd.pp.export_metabo_analyst)
        nd.pp.export_r_dolphin = "rDolphinExport.csv"
        nd.pp.export_samples_in_rows_cols = 0
        nd.pp.export_method = 3
        nd.export_data_set()
        self.assertEqual(os.path.isfile(nd.pp.export_r_dolphin_path + os.sep + nd.pp.export_r_dolphin), True)
        os.remove(nd.pp.export_r_dolphin_path + os.sep + nd.pp.export_r_dolphin)
        nd.pp.export_batman = "batmanExport.txt"
        nd.pp.export_samples_in_rows_cols = 0
        nd.pp.export_method = 4
        nd.export_data_set()
        self.assertEqual(os.path.isfile(nd.pp.export_batman_path + os.sep + nd.pp.export_batman), True)
        os.remove(nd.pp.export_batman_path + os.sep + nd.pp.export_batman)
        nd.pp.export_bruker = "brukerExport"
        nd.pp.export_samples_in_rows_cols = 0
        nd.pp.export_method = 5
        nd.export_data_set()
        self.assertEqual(
            os.path.isfile(nd.pp.export_bruker_path + os.sep + nd.pp.export_bruker + os.sep + '1' + os.sep + 'acqus'),
            True)
        shutil.rmtree(nd.pp.export_bruker_path + os.sep + nd.pp.export_bruker)

    def test_noise_filtering(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.pp.noise_start = 10.0
        nd.pp.noise_end = 10.5
        nd.pp.noise_threshold = 8.0
        nd.pp.th_line_width = 2.0
        nd.pp.flag_noise_filtering = True
        nd.data_pre_processing()
        pts = len(nd.nmrdat[0][0].spc[0]) - nd.nmrdat[0][0].ppm2points([9.5, 11.5], 0)
        pts2 = range(int(min(pts)), int(max(pts)))
        ssum = np.sum(nd.nmrdat[0][0].spc[0][pts2].real)
        self.assertEqual(ssum, 0.0)

    def test_bucket_spectra(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.pp.bucket_points = 64
        nd.pp.flag_bucket_spectra = True
        nd.data_pre_processing()
        pts = len(nd.nmrdat[0][0].spc[0])
        self.assertEqual(pts, 2850)

    def test_compress_buckets(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.nmrdat[0][0].proc.ph0[0] = -18.49143139318872
        nd.nmrdat[0][0].proc.ph1[0] = -0.7082407011136772
        nd.ft()
        nd.pp.compress_start = [6.5448]
        nd.pp.compress_end = [7.4705]
        nd.pp.flag_compress_buckets = True
        nd.data_pre_processing()
        pts = nd.nmrdat[0][0].ppm2points([nd.pp.compress_end[0], nd.pp.compress_start[0]], 0)
        npts = len(nd.nmrdat[0][0].spc[0])
        spc = nd.nmrdat[0][0].spc[0][npts - pts[0]:npts - pts[1]].real
        self.assertEqual(spc.max(), spc.sum())

    def test_avoid_negative_values(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.nmrdat[0][0].proc.ph0[0] = -18.49143139318872
        nd.nmrdat[0][0].proc.ph1[0] = -0.7082407011136772
        nd.ft()
        nd.pp.avoid_negative_values = True
        nd.data_pre_processing()
        self.assertEqual(nd.nmrdat[0][0].spc[0].real.min(), 0.0)

    def test_segmental_alignment(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "10"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.nmrdat[0][0].proc.ph0[0] = -85.94556
        nd.nmrdat[0][0].proc.ph1[0] = 6.64201
        nd.nmrdat[0][0].proc.lb[0] = 2.0
        nd.ft()
        nd.nmrdat[0][0].autobaseline1d()
        nd.nmrdat[0][0].auto_ref()
        e_name = "11"  # 1D NMR data in exp 1
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.nmrdat[0][1].proc.ph0[0] = -79.29513383337031
        nd.nmrdat[0][1].proc.ph1[0] = 7.790498486770286
        nd.nmrdat[0][1].proc.lb[0] = 2.0
        nd.e = 1
        nd.ft()
        nd.nmrdat[0][0].autobaseline1d()
        nd.nmrdat[0][1].autobaseline1d()
        nd.nmrdat[0][0].auto_ref()
        nd.nmrdat[0][1].auto_ref()
        start_ppm = 8.1091
        end_ppm = 8.1825
        nd.pp.seg_start = [start_ppm]
        nd.pp.seg_end = [end_ppm]
        nd.pp.flag_segmental_alignment = True
        nd.pp.seg_align_ref_spc = 1
        nd.data_pre_processing()
        pts1 = len(nd.nmrdat[0][0].spc[0]) - nd.nmrdat[0][0].ppm2points([end_ppm, start_ppm], 0)
        spc1 = nd.nmrdat[0][0].spc[0][pts1[0]:pts1[1]].real
        spc2 = nd.nmrdat[0][1].spc[0][pts1[0]:pts1[1]].real
        self.assertEqual(np.argmax(spc1), np.argmax(spc2))

    def test_scale_spectra(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "10"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.nmrdat[0][0].proc.ph0[0] = -85.94556
        nd.nmrdat[0][0].proc.ph1[0] = 6.64201
        nd.nmrdat[0][0].proc.lb[0] = 0.3
        nd.ft()
        nd.nmrdat[0][0].autobaseline1d()
        nd.nmrdat[0][0].auto_ref()
        e_name = "11"  # 1D NMR data in exp 1
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.nmrdat[0][1].proc.ph0[0] = -79.29513383337031
        nd.nmrdat[0][1].proc.ph1[0] = 7.790498486770286
        nd.nmrdat[0][1].proc.lb[0] = 0.3
        nd.e = 1
        nd.ft()
        nd.pp.flag_scale_spectra = True
        nd.pp.scale_spectra_ref_spc = 1
        self.assertAlmostEqual(nd.nmrdat[0][1].spc[0].real.max() / 5746500922.356481, 1.0, 1)
        nd.pp.scale_pqn = True
        nd.data_pre_processing()
        self.assertAlmostEqual(nd.nmrdat[0][1].spc[0].real.max() / 4250520058.933967, 1.0, 1)
        nd.reset_data_pre_processing()
        self.assertAlmostEqual(nd.nmrdat[0][0].spc[0].real.max() / 5673105467.238647, 1.0, 1)
        nd.pp.scale_pqn = False
        nd.pp.flag_scale_spectra = True
        nd.pp.preserve_overall_scale = True
        nd.data_pre_processing()
        self.assertAlmostEqual(nd.nmrdat[0][0].spc[0].real.max() / 5673105467.2386465, 1.0, 1)
        nd.reset_data_pre_processing()
        nd.pp.preserve_overall_scale = False
        nd.data_pre_processing()
        self.assertAlmostEqual(nd.nmrdat[0][0].spc[0].real.max() / 0.006644593036069566, 1.0, 3)

    def test_pjres(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "2"  # 2D JresNMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.ft()
        nd.pjres(2, 'skyline')
        self.assertEqual(len(nd.nmrdat[1][0].spc[0]), 8192)

    def test_pre_proc_init(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.pp.init(1)
        self.assertEqual(nd.pp.class_select, "1")
        # end test_preProcInit

    def test_pre_proc_init_plot_colours(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.pp.cf.read_config()
        nd.pp.init_plot_colours()
        if nd.pp.cf.mode == 'dark' or (nd.pp.cf.mode == 'system' and darkdetect.isDark()):
            self.assertEqual(nd.pp.plot_colours[0], (1.0, 1.0, 0.0))
            self.assertEqual(nd.pp.plot_colours[-1], (0.6, 0.3, 0.6))
        else:
            self.assertEqual(nd.pp.plot_colours[0], (0.0, 0.0, 0.4))
            self.assertEqual(nd.pp.plot_colours[-1], (0.8, 0.5, 0.8))

        # end test_pre_proc_init_plot_colours

    def test_pre_proc_set_alpha(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        dd = nd.pp.set_alpha(0.5)
        self.assertEqual(nd.pp.alpha, 0.5)
        # end test_pre_proc_set_alpha

    def test_read_spc(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc()  # check if Bruker data can be read
        self.assertEqual(len(nd.nmrdat[0]), 0)
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.data_set_name = p_name
        nd.read_spc()  # check if Bruker data can be read
        self.assertEqual(len(nd.nmrdat[0]), 0)
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()  # check if Bruker data can be read
        self.assertEqual(len(nd.nmrdat[0][0].fid[0]), 32768)  # check number of data points in fid
        self.assertEqual(len(nd.nmrdat[0][0].spc[0]), 65536)  # check number of data points in 1r
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        self.assertEqual(len(nd.nmrdat[0][0].fid[0]), 32768)  # check number of data points in fid
        self.assertEqual(len(nd.nmrdat[0][0].spc[0]), 65536)  # check number of data points in 1r
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.data_set_name = p_name
        nd.data_set_number = e_name
        nd.read_spc()  # check if Bruker data can be read
        self.assertEqual(len(nd.nmrdat[0][0].fid[0]), 32768)  # check number of data points in fid
        self.assertEqual(len(nd.nmrdat[0][0].spc[0]), 65536)  # check number of data points in 1r

    def test_read_spcs(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = ["1", "1"]  # 1D data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spcs([p_name], e_name)  # check if Bruker data can be read
        self.assertEqual(len(nd.nmrdat[0][0].spc[0]), 65536)
        self.assertEqual(len(nd.nmrdat[0][1].spc[0]), 65536)

    def test_read_nmrpipe_spc(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "5"  # 2D HSQC NMRPipe data in exp 5
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_nmrpipe_spc(p_name, e_name)
        self.assertEqual(len(nd.nmrdat[0][0].spc), 16384)
        self.assertEqual(len(nd.nmrdat[0][0].spc[0]), 922)

    def test_read_nmrpipe_spcs(self):
        p_name = [os.path.join(os.path.dirname(__file__), "data", "nmrData")]  # directory of test data set
        e_name = ["5"]  # 2D HSQC NMRPipe data in exp 5
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_nmrpipe_spcs(p_name, e_name)
        self.assertEqual(len(nd.nmrdat[0][0].spc), 16384)
        self.assertEqual(len(nd.nmrdat[0][0].spc[0]), 922)

    def test_reset_data_pre_processing(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.pp.exclude_start = [8.9]
        nd.pp.exclude_end = [11.0]
        nd.pp.flag_exclude_region = True
        nd.data_pre_processing()
        nd.reset_data_pre_processing()
        self.assertEqual(len(np.where(nd.nmrdat[0][0].spc[0].real == 0)[0]), 0)

    def test_save(self):
        nd = nmrDataSet.NmrDataSet()
        data_set_path = os.path.join(os.path.dirname(__file__), "data", "loadData.mlpy")
        nd.load(data_set_path)
        data_set_path = os.path.join(os.path.dirname(__file__), "data", "saveData.mlpy")
        nd.save(data_set_path)
        self.assertEqual(os.path.exists(data_set_path), True)
        shutil.rmtree(data_set_path)

    def test_set_gb(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.set_gb([15.0])
        self.assertEqual(nd.nmrdat[0][0].proc.gb[0], 15.0)

    def test_set_lb(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.set_lb([15.0])
        self.assertEqual(nd.nmrdat[0][0].proc.lb[0], 15.0)

    def test_set_ph_from_exp(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.set_ph_from_exp()
        self.assertEqual(nd.nmrdat[0][0].proc.ph0[0], 22.0)
        self.assertEqual(nd.nmrdat[0][0].proc.ph1[0], 22.0)

    def test_ph0(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.set_ph0([15.0])
        self.assertEqual(nd.nmrdat[0][0].proc.ph0[0], 15.0)

    def test_ph1(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.set_ph1([15.0])
        self.assertEqual(nd.nmrdat[0][0].proc.ph1[0], 15.0)

    def test_set_ssb(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.set_ssb([1.0])
        self.assertEqual(nd.nmrdat[0][0].proc.ssb[0], 1.0)

    def test_set_window_type(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.set_window_type([2])
        self.assertEqual(nd.nmrdat[0][0].proc.window_type[0], 2)

    def test_shift_ref(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
        nd.shift_ref()
        self.assertEqual(nd.nmrdat[0][0].ref_shift[0], 0.0)

    def test_add_peak(self):
        f_name = os.path.join(os.path.dirname(__file__), "data", "loadData.mlpy")  # directory of test data set
        nd = nmrDataSet.NmrDataSet()
        nd.load(f_name)
        nd.add_peak(np.array([0.01, -0.01]), ['TMSP'])
        self.assertAlmostEqual(nd.nmrdat[0][0].peak_max_ppm, 0.0)
    # end

    def test_set_peak(self):
        f_name = os.path.join(os.path.dirname(__file__), "data", "loadData.mlpy")  # directory of test data set
        nd = nmrDataSet.NmrDataSet()
        nd.load(f_name)
        nd.set_peak(np.array([0.01]), np.array([-0.01]), np.array(['TMSP']), n_protons=[9])
        self.assertAlmostEqual(nd.nmrdat[0][0].peak_max_ppm[0], 0.0)
    # end

    def test_clear_peak(self):
        f_name = os.path.join(os.path.dirname(__file__), "data", "loadData.mlpy")  # directory of test data set
        nd = nmrDataSet.NmrDataSet()
        nd.load(f_name)
        nd.set_peak(np.array([0.01]), np.array([-0.01]), np.array(['TMSP']), n_protons=[9])
        nd.clear_peak()
        self.assertEqual(len(nd.nmrdat[0][0].peak_max_ppm), 0)
    # end

    def test_create_titles(self):
        excel_name = os.path.join(os.path.dirname(__file__), "data", "sampleTitleSpreadSheet.xlsx")
        xls = pd.read_excel(excel_name).fillna('')
        f_name = os.path.join(os.path.dirname(__file__), "data", "loadData.mlpy")  # directory of test data set
        nd = nmrDataSet.NmrDataSet()
        nd.load(f_name)
        rack_label = 'Rack'
        pos_label = 'Position'
        dataset_label = 'dataset'
        replace_title = True
        nd.create_titles(xls, dataset_label, pos_label, rack_label, replace_title, excel_name)
        nd.nmrdat[0][0].title.index('sample : control')
        # end test_create_title

    def test_read_spc(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()
        nd.read_spc(p_name, e_name, 1)
        self.assertEqual(len(nd.nmrdat[0][0].fid[0]), 32768)  # check number of data points in fid
        self.assertEqual(len(nd.nmrdat[0][0].spc[0]), 65536)  # check number of data points in 1r
        # end test_read_spc

    def test_read_spcs(self):
        p_name = [os.path.join(os.path.dirname(__file__), "data", "nmrData")]  # directory of test data set
        e_name = ["1", "6"]  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()
        nd.read_spcs(p_name, e_name, 1)
        self.assertEqual(len(nd.nmrdat[0][0].fid[0]), 32768)  # check number of data points in fid
        self.assertEqual(len(nd.nmrdat[0][0].spc[0]), 65536)  # check number of data points in 1r
        self.assertEqual(len(nd.nmrdat[0]), 2)
        # end test_read_spc

    def test_pre_proc_init(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()
        nd.read_spc(p_name, e_name, 1)
        self.assertEqual(len(nd.pp.plot_select), 0)
        nd.pre_proc_init()
        self.assertEqual(len(nd.pp.plot_select), 1)
        # end test_pre_proc_init

    def test_read_title_file_information_excel(self):
        excel_name = os.path.join(os.path.dirname(__file__), "data", "sampleTitleSpreadSheet.xlsx")
        nd = nmrDataSet.NmrDataSet()
        xls = nd.read_title_file_information_excel(excel_name)
        self.assertEqual(xls['Rack'][0], 1)
        # end

    def test_reference1d_all(self):
        p_name = [os.path.join(os.path.dirname(__file__), "data", "nmrData")]  # directory of test data set
        e_name = ["1", "6"]  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()
        nd.read_spcs(p_name, e_name, 1)
        nd.auto_ref_all()
        ref_point1 = nd.nmrdat[0][0].ref_point[0]
        ref_point2 = nd.nmrdat[0][1].ref_point[0]
        nd.reference1d_all(0.0, 10.0)
        self.assertEqual(nd.nmrdat[0][0].ref_shift[0], 10.0)  # check reference shift in spectrum 1
        self.assertEqual(nd.nmrdat[0][1].ref_shift[0], 10.0)  # check reference shift in spectrum 2
        self.assertEqual(nd.nmrdat[0][0].ref_point[0], ref_point1)  # check reference point in spectrum 1
        self.assertEqual(nd.nmrdat[0][1].ref_point[0], ref_point2)  # check reference point in spectrum 2
        # end test_reference1d_all

    def test_set_standard_plot_colours(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        e_name = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()
        nd.read_spc(p_name, e_name, 1)
        nd.cf.mode = 'light'
        nd.set_standard_plot_colours()
        self.assertEqual(nd.nmrdat[0][0].display.pos_col_rgb, (0.0, 0.0, 1.0))
        nd.cf.mode = 'dark'
        nd.set_standard_plot_colours()
        self.assertEqual(nd.nmrdat[0][0].display.pos_col_rgb, (0.8, 0.8, 1.0))
        # end test_set_standard_plot_colours

    def test_spline_correct(self):
        f_name = os.path.join(os.path.dirname(__file__), "data", "loadData.mlpy")  # directory of test data set
        nd = nmrDataSet.NmrDataSet()
        nd.load(f_name)
        nd.nmrdat[0][0].spline_baseline.baseline_points = [11.787481124553866, -2.2087244721275536]
        nd.nmrdat[0][0].add_baseline_points()
        baseline = nd.nmrdat[0][0].calc_spline_baseline()
        nd.spline_correct()
        # end test_calc_spline_baseline

    def test_variance_stabilisation(self):
        ds_name = os.path.join(os.path.dirname(__file__), "data", "pre_proc_test.mlpy")  # directory of test data set
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.load(ds_name)  # check if Bruker data can be read
        nd.pp.flag_variance_stabilisation = True
        nd.pp.auto_scaling = False
        nd.pp.pareto_scaling = True
        nd.pp.g_log_transform = False
        max1 = np.max(nd.nmrdat[0][0].spc[0].real)
        nd.data_pre_processing()
        max2 = np.max(nd.nmrdat[0][0].spc[0].real)
        self.assertAlmostEqual(max1, 10185302564.706278, places=1)
        self.assertAlmostEqual(max2, 521808.8253127394, places=1)
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.load(ds_name)  # check if Bruker data can be read
        nd.pp.flag_variance_stabilisation = True
        nd.pp.auto_scaling = True
        nd.pp.pareto_scaling = False
        nd.pp.g_log_transform = False
        max1 = np.max(nd.nmrdat[0][0].spc[0].real)
        nd.data_pre_processing()
        max2 = np.max(nd.nmrdat[0][0].spc[0].real)
        self.assertAlmostEqual(max1, 10185302564.706278, places=1)
        self.assertAlmostEqual(max2, 1.7320373446154178, places=1)
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.load(ds_name)  # check if Bruker data can be read
        nd.pp.flag_variance_stabilisation = True
        nd.pp.auto_scaling = False
        nd.pp.pareto_scaling = False
        nd.pp.g_log_transform = True
        max1 = np.max(nd.nmrdat[0][0].spc[0].real)
        nd.data_pre_processing()
        max2 = np.max(nd.nmrdat[0][0].spc[0].real)
        self.assertAlmostEqual(max1, 10185302564.706278, places=1)
        self.assertAlmostEqual(max2, 9.727370430885156, places=1)

    def test_export_hsqc_data(self):
        p_name = os.path.join(os.path.dirname(__file__), "data", "nmrData")
        e_name = "5"
        export_excel_path = os.path.expanduser("~")
        export_excel = "testHsqcExport.xlsx"
        nd = nmrDataSet.NmrDataSet()
        nd.read_nmrpipe_spc(p_name, e_name, "test.dat")
        nd.nmrdat[0][0].proc.ph0 = [153.05550133, -1.55127223, 0.0]
        nd.nmrdat[0][0].proc.ph1 = [165.12883873, 6.06083814, 0.0]
        nd.nmrdat[0][0].ref_shift = [1.314, 22.8972, 0.0]
        nd.nmrdat[0][0].ref_point = [143, 2355, 0]
        nd.nmrdat[0][0].calc_ppm()
        nd.nmrdat[0][0].phase2d(nd.nmrdat[0][0].proc.ph0[0], nd.nmrdat[0][0].proc.ph1[0], 0)
        nd.nmrdat[0][0].phase2d(nd.nmrdat[0][0].proc.ph0[1], nd.nmrdat[0][0].proc.ph1[1], 1)
        nd.nmrdat[0][0].autobaseline2d()
        metabolite_name = 'L-AsparticAcid'
        cur_peak = 2
        nd.nmrdat[0][0].hsqc.read_metabolite_information(metabolite_name)
        nd.nmrdat[0][0].hsqc.set_metabolite_information(metabolite_name, nd.nmrdat[0][0].hsqc.metabolite_information)
        nd.nmrdat[0][0].hsqc.cur_metabolite = metabolite_name
        nd.nmrdat[0][0].hsqc.cur_peak = cur_peak
        nd.nmrdat[0][0].hsqc.set_peak_information()
        nd.nmrdat[0][0].autopick_hsqc()
        nd.nmrdat[0][0].autofit_hsqc()
        nd.export_hsqc_data(os.path.join(export_excel_path, export_excel))
        self.assertEqual(os.path.isfile(os.path.join(export_excel_path, export_excel)), True)
        os.remove(os.path.join(export_excel_path, export_excel))
        # end

if __name__ == "__main__":
    unittest.main()
