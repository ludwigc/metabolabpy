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
        nd.autophase1d()

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
        self.assertEqual(nd.nmrdat[0][0].ref_point[0], 152)
        self.assertEqual(nd.nmrdat[0][0].ref_shift[0], 0.0)
        self.assertEqual(nd.nmrdat[0][0].ref_point[1], 1792)
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
        pts2 = range(int(min(pts)), int(max(pts)))
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
        self.assertEqual(pts, 1024)

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
        self.assertAlmostEqual(nd.nmrdat[0][1].spc[0].real.max() / 10331062121.45446, 1.0, 1)
        nd.pp.scale_pqn = True
        nd.data_pre_processing()
        self.assertAlmostEqual(nd.nmrdat[0][1].spc[0].real.max() / 7644944702.898952, 1.0, 1)
        nd.reset_data_pre_processing()
        self.assertAlmostEqual(nd.nmrdat[0][0].spc[0].real.max() / 9406699558.794441, 1.0, 1)
        nd.pp.scale_pqn = False
        nd.pp.flag_scale_spectra = True
        nd.pp.preserve_overall_scale = True
        nd.data_pre_processing()
        self.assertAlmostEqual(nd.nmrdat[0][0].spc[0].real.max() / 11501776512.517786, 1.0, 1)
        nd.reset_data_pre_processing()
        nd.pp.preserve_overall_scale = False
        nd.data_pre_processing()
        self.assertAlmostEqual(nd.nmrdat[0][0].spc[0].real.max() / 0.011113980938806603, 1.0, 1)

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
        if nd.pp.cf.mode == 'dark':
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
        nd.read_spc(p_name, e_name)  # check if Bruker data can be read
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


if __name__ == "__main__":
    unittest.main()
