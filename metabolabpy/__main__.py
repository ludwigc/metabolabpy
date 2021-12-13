#!/usr/bin/env python
import sys  # pragma: no cover
import matplotlib  # pragma: no cover

matplotlib.use("Agg")
matplotlib.rcParams['agg.path.chunksize'] = 64000 #64_000_000_000

try:
    from PySide2.QtUiTools import QUiLoader  # pragma: no cover
    from PySide2.QtCore import QFile  # pragma: no cover
    from PySide2.QtCore import QCoreApplication  # pragma: no cover
    from PySide2.QtWidgets import *  # pragma: no cover
    from PySide2 import QtWidgets  # pragma: no cover
    from PySide2.QtGui import *  # pragma: no cover
    from PySide2 import QtGui  # pragma: no cover
    from PySide2 import QtCore  # pragma: no cover
    from PySide2.QtWidgets import QFileDialog  # pragma: no cover
    from PySide2.QtCore import SIGNAL  # pragma: no cover
    from PySide2.QtWebEngineWidgets import QWebEngineView #, QWebEngineProfile, QWebEnginePage, \
    #    QWebEngineSettings  # pragma: no cover
    from PySide2.QtCore import QUrl, Qt  # pragma: no cover
    from PySide2.QtWebEngineCore import QWebEngineUrlSchemeHandler  # pragma: no cover
    import PySide2  # pragma: no cover
    import qtmodern.styles  # pragma: no cover
except:
    pass

import matplotlib.pyplot as pl  # pragma: no cover

if "linux" in sys.platform:  # pragma: no cover
    gui_env = ['TkAgg', 'GTKAgg', 'Qt5Agg', 'WXAgg']  # pragma: no cover
elif sys.platform == "darwin":  # pragma: no cover
    try:  # pragma: no cover
        gui_env = ['Qt5Agg']  # pragma: no cover
    except ImportError:  # pragma: no cover
        gui_env = ['TkAgg', 'GTKAgg', 'Qt5Agg', 'WXAgg']  # pragma: no cover

    # matplotlib.use("Qt5Agg")
else:  # pragma: no cover
    gui_env = ['TkAgg', 'GTKAgg', 'Qt5Agg', 'WXAgg']  # pragma: no cover

if sys.platform != "win32":  # pragma: no cover
    for gui in gui_env:  # pragma: no cover
        try:  # pragma: no cover
            matplotlib.use(gui, warn=False, force=True)  # pragma: no cover
            break  # pragma: no cover
        except:  # pragma: no cover
            continue  # pragma: no cover

try:
    from matplotlib.backends.backend_qt5agg import (FigureCanvasQTAgg as FigureCanvas,
                                                    NavigationToolbar2QT as NavigationToolbar)  # pragma: no cover
except:
    pass

from matplotlib.figure import Figure  # pragma: no cover
import argparse  # pragma: no cover
# from PySide2.QtUiTools import QUiLoader  # pragma: no cover
# from PySide2.QtCore import QFile  # pragma: no cover
# from PySide2.QtCore import QCoreApplication  # pragma: no cover
# from PySide2.QtWidgets import *  # pragma: no cover
# from PySide2 import QtWidgets  # pragma: no cover
# from PySide2.QtGui import *  # pragma: no cover
# from PySide2 import QtGui  # pragma: no cover
# from PySide2 import QtCore  # pragma: no cover
# from PySide2.QtWidgets import QFileDialog # pragma: no cover
# from PySide2.QtCore import SIGNAL  # pragma: no cover
from time import sleep  # pragma: no cover

try:  # pragma: no cover
    import pyautogui  # pragma: no cover
except:  # pragma: no cover
    pass  # pragma: no cover

import numpy as np  # pragma: no cover
import io  # pragma: no cover
from metabolabpy.nmr import nmrDataSet  # pragma: no cover
from metabolabpy.GUI import phCorr  # pragma: no cover
import time  # pragma: no cover
##import platform  # pragma: no cover
import math  # pragma: no cover
from metabolabpy.nmr import nmrConfig  # pragma: no cover
import os  # pragma: no cover
import traceback  # pragma: no cover
import shutil  # pragma: no cover
import scipy.io  # pragma: no cover
##import inspect
from io import StringIO
import contextlib
import zipfile
# from notebook import notebookapp
# import multiprocess
import subprocess
# import jupyterthemes
import itertools
import xlsxwriter
from string import ascii_uppercase

## import pandas as pd                       # pragma: no cover

try:
    # ------------------ MplWidget ------------------
    class MplWidget(QWidget):  # pragma: no cover

        def __init__(self, parent=None):
            QWidget.__init__(self, parent)
            fig = Figure()
            self.canvas = FigureCanvas(fig)
            vertical_layout = QVBoxLayout()
            vertical_layout.addWidget(self.canvas)
            home = NavigationToolbar.home
            def new_home(self, *args, **kwargs):
                self.canvas.axes.autoscale()
                self.canvas.draw()
                self.canvas.toolbar.update()
                home(self, *args, **kwargs)

            NavigationToolbar.home = new_home
            self.toolbar = NavigationToolbar(self.canvas, self)
            #self.toolbar._actions['back'].setEnabled(False)
            vertical_layout.addWidget(self.toolbar)

            self.canvas.axes = self.canvas.figure.add_subplot(111)
            self.setLayout(vertical_layout)
            self.ph_corr = phCorr.PhCorr()
            # end __init__


    # ------------------ MplWidget ------------------

    # ------------------ MplWidget ------------------
    class QWebEngineView2(QWebEngineView):
        def print_cmd(self):
            print("QWebEngineView2")

    # ------------------ MplWidget2 ------------------
except:
    pass


class main_w(object):  # pragma: no cover
    def __init__(self):
        self.exited_peak_picking = False
        self.__version__ = '0.6.51'
        self.zoom_was_on = True
        self.pan_was_on = False
        self.std_pos_col1 = (0.0, 0.0, 1.0)
        self.std_neg_col1 = (1.0, 0.0, 0.0)
        self.std_pos_col2 = (0.8, 0.8, 1.0)
        self.std_neg_col2 = (1.0, 0.8, 0.8)
        self.n_clicks = 1
        self.cur_clicks = 0
        self.xy = [[]]
        self.xdata = []
        self.ydata = []
        self.nd = nmrDataSet.NmrDataSet()
        self.ph_corr = phCorr.PhCorr()
        # load ui; create w
        f_name = os.path.join(os.path.dirname(__file__), "ui", "metabolabpy_mainwindow.ui")
        self.file = QFile(f_name)
        self.file.open(QFile.ReadOnly)
        self.loader = QUiLoader()
        self.loader.registerCustomWidget(QWebEngineView2)
        self.loader.registerCustomWidget(MplWidget)
        # self.loader.registerCustomWidget(MplWidget2)
        self.w = self.loader.load(self.file)
        self.zoom = False

        self.hide_pre_processing()
        self.hide_peak_picking()
        self.w.preprocessing.setVisible(False)
        self.w.peakPicking.setVisible(False)
        self.w.preProcPeak.setVisible(False)
        self.w.hsqcAnalysis.setVisible(False)
        self.w.multipletAnalysis.setVisible(False)
        self.w.isotopomerAnalysis.setVisible(False)
        self.w.nmrSpectrum.setTabEnabled(1, False)
        self.w.nmrSpectrum.setTabEnabled(2, False)
        self.w.nmrSpectrum.setTabEnabled(3, False)
        self.w.nmrSpectrum.setTabEnabled(4, False)
        self.w.nmrSpectrum.setStyleSheet(
            "QTabBar::tab::disabled {width: 0; height: 0; margin: 0; padding: 0; border: none;} ")
        # connections
        # self.w.rDolphinExport.clicked.connect(self.setrDolphinExport)
        self.w.exportPath.textChanged.connect(self.set_export_path)
        self.w.invertMatrix_1.stateChanged.connect(self.set_invert)
        self.w.invertMatrix_2.stateChanged.connect(self.set_invert)
        self.w.exportFileName.textChanged.connect(self.set_export_file_name)
        self.w.exportDelimiterTab.toggled.connect(self.set_export_delimiter_tab)
        self.w.exportCharacter.textChanged.connect(self.set_export_character)
        self.w.samplesInComboBox.currentIndexChanged.connect(self.set_samples_in_combo_box)
        self.w.runPreProcessingButton.clicked.connect(self.data_pre_processing)
        self.w.resetPreProcessingButton.clicked.connect(self.reset_data_pre_processing)
        self.w.avoidNegValues.stateChanged.connect(self.set_avoid_neg_values)
        self.w.excludeRegion.stateChanged.connect(self.set_exclude_region)
        self.w.segmentalAlignment.stateChanged.connect(self.set_segmental_alignment)
        self.w.compressBuckets.stateChanged.connect(self.set_compress_buckets)
        self.w.noiseFiltering.stateChanged.connect(self.set_noise_filtering)
        self.w.bucketSpectra.stateChanged.connect(self.set_bucket_spectra)
        self.w.scaleSpectraRefSpc.valueChanged.connect(self.change_scale_spectra_ref_spc)
        self.w.segAlignRefSpc.valueChanged.connect(self.change_seg_align_ref_spc)
        self.w.scaleSpectra.stateChanged.connect(self.set_scale_spectra)
        self.w.pqnButton.clicked.connect(self.set_pqn_tsa_scaling)
        self.w.tsaButton.clicked.connect(self.set_pqn_tsa_scaling)
        self.w.autoScaling.clicked.connect(self.set_variance_stabilisation_options)
        self.w.paretoScaling.clicked.connect(self.set_variance_stabilisation_options)
        self.w.gLogTransform.clicked.connect(self.set_variance_stabilisation_options)
        self.w.varianceStabilisation.stateChanged.connect(self.set_variance_stabilisation)
        self.w.exportDataSet.stateChanged.connect(self.set_export_data_set)
        self.w.excludeRegionTW.cellChanged.connect(self.set_exclude_pre_proc)
        self.w.segAlignTW.cellChanged.connect(self.set_seg_align_pre_proc)
        self.w.selectClassTW.itemSelectionChanged.connect(self.set_plot_pre_proc)
        self.w.selectClassTW.cellChanged.connect(self.set_change_pre_proc)
        self.w.excludeClearButton.clicked.connect(self.select_clear_exclude_pre_proc)
        self.w.segAlignClearButton.clicked.connect(self.select_clear_seg_align_pre_proc)
        self.w.compressClearButton.clicked.connect(self.select_clear_compress_pre_proc)
        self.w.excludeAddButton.clicked.connect(self.select_add_exclude_pre_proc)
        self.w.segAlignAddButton.clicked.connect(self.select_add_seg_align_pre_proc)
        self.w.compressAddButton.clicked.connect(self.select_add_compress_pre_proc)
        self.w.selectAllButton.clicked.connect(self.select_all_pre_proc)
        self.w.selectEvenButton.clicked.connect(self.select_even_pre_proc)
        self.w.selectOddButton.clicked.connect(self.select_odd_pre_proc)
        self.w.selectClassButton.clicked.connect(self.select_class_pre_proc)
        self.w.selectClassLE.textChanged.connect(self.select_class_pre_proc)
        self.w.cmdLine.returnPressed.connect(self.exec_cmd)
        self.w.noiseThresholdLE.textChanged.connect(self.set_noise_reg_pre_proc)
        self.w.noiseRegionStartLE.textChanged.connect(self.set_noise_reg_pre_proc)
        self.w.noiseRegionEndLE.textChanged.connect(self.set_noise_reg_pre_proc)
        self.w.thLineWidthLE.textChanged.connect(self.set_noise_reg_pre_proc)
        self.w.bucketPpmLE.textChanged.connect(self.set_bucket_ppm_pre_proc)
        self.w.bucketDataPointsLE.textChanged.connect(self.set_bucket_points_pre_proc)
        self.w.actionVertical_AutoScale.triggered.connect(self.vertical_auto_scale)
        self.w.actionZoom.triggered.connect(self.set_zoom)
        self.w.actionPan.triggered.connect(self.set_pan)
        self.w.actionShow_Next_Tab.triggered.connect(self.next_tab)
        self.w.actionShow_Previous_Tab.triggered.connect(self.previous_tab)
        self.w.actionPlot_spc.triggered.connect(self.plot_spc)
        self.w.actionSave.triggered.connect(self.save_button)
        self.w.actionLoad.triggered.connect(self.load_button)
        self.w.actionOpen_NMRPipe.triggered.connect(self.read_nmrpipe_spc)
        self.w.actionActivate_Command_Line.triggered.connect(self.activate_command_line)
        self.w.actionPrevious_command.triggered.connect(self.previous_command)
        self.w.actionNext_command.triggered.connect(self.next_command)
        self.w.actionCorrect_Phase.triggered.connect(self.start_stop_ph_corr)
        # self.w.actionZoomCorrect_Phase.triggered.connect(self.zoom_ph_corr)
        self.w.zoomPhCorr1d.clicked.connect(self.zoom_ph_corr)
        self.w.exitZoomPhCorr1d.clicked.connect(self.zoom_ph_corr)
        self.w.exitPhCorr1d.clicked.connect(self.start_stop_ph_corr)
        self.w.actionClear.triggered.connect(self.clear)
        self.w.lambdaLE.textChanged.connect(self.set_var_lambda)
        self.w.y0LE.textChanged.connect(self.set_var_y0)
        self.w.actionRead_NMR_Spectrum.triggered.connect(self.read_nmr_spc)
        self.w.preprocessing.stateChanged.connect(self.set_pre_processing)
        self.w.peakPicking.stateChanged.connect(self.set_peak_picking)
        self.w.peakAddButton.clicked.connect(self.add_peak)
        self.w.peakClearButton.clicked.connect(self.clear_peak)
        self.w.peakWidget.cellChanged.connect(self.set_add_peak)
        self.w.peakExportButton.clicked.connect(self.export_peak)
        self.w.intAllDS.stateChanged.connect(self.set_datasets_exps)
        self.w.intAllExps.stateChanged.connect(self.set_datasets_exps)
        self.w.exportFormatCB.currentIndexChanged.connect(self.set_datasets_exps)
        self.w.hsqcAnalysis.stateChanged.connect(self.set_hsqc_analysis)
        self.w.multipletAnalysis.stateChanged.connect(self.set_multiplet_analysis)
        self.w.isotopomerAnalysis.stateChanged.connect(self.set_isotopomer_analysis)
        self.w.preserveOverallScale.stateChanged.connect(self.set_preserve_overall_scale)
        self.w.actionReset.triggered.connect(self.reset_plot)
        self.w.actionShow_NMR_Spectrum.triggered.connect(self.show_nmr_spectrum)
        self.w.actionSetup_Processing_Parameters.triggered.connect(self.setup_processing_parameters)
        self.w.actionShow_Display_Parameters.triggered.connect(self.show_display_parameters)
        self.w.actionShow_Acquisition_Parameters.triggered.connect(self.show_acquisition_parameters)
        self.w.actionShow_Title_File_Information.triggered.connect(self.show_title_file_information)
        self.w.actionShow_pulseProgram.triggered.connect(self.show_pulse_program)
        self.w.actionFourier_Transform.triggered.connect(self.ft)
        self.w.actionScript_Editor.triggered.connect(self.script_editor)
        self.w.actionChange_to_next_Exp.triggered.connect(self.change_to_next_exp)
        self.w.actionChange_to_previous_Exp.triggered.connect(self.change_to_previous_exp)
        self.w.actionChange_to_next_DS.triggered.connect(self.change_to_next_ds)
        self.w.actionChange_to_previous_DS.triggered.connect(self.change_to_previous_ds)
        self.w.exampleScripts.view().pressed.connect(self.load_example_script)
        self.w.actionAutomatic_Phase_Correction.triggered.connect(self.autophase1d)
        self.w.actionAutomatic_Baseline_Correction.triggered.connect(self.autobaseline1d)
        self.w.actionScale_2D_Spectrum_Up.triggered.connect(self.scale_2d_spectrum_up)
        self.w.actionScale_2D_Spectrum_Down.triggered.connect(self.scale_2d_spectrum_down)
        self.w.actionScale_all_2D_Spectra_Up.triggered.connect(self.scale_all_2d_spectra_up)
        self.w.actionScale_all_2D_Spectra_Down.triggered.connect(self.scale_all_2d_spectra_down)
        self.w.actionSelect_All.triggered.connect(self.select_plot_all)
        self.w.actionClear_All.triggered.connect(self.select_plot_clear)
        self.w.actionConsole.triggered.connect(self.show_console)
        self.w.actionShow_SplashScreen.triggered.connect(self.splash)
        self.w.actionHelp.triggered.connect(self.show_help)
        self.w.actionToggle_FullScreen.triggered.connect(self.show_main_window)
        self.w.setBox.valueChanged.connect(self.change_data_set_exp)
        self.w.expBox.valueChanged.connect(self.change_data_set_exp)
        self.w.posCol.currentIndexChanged.connect(self.get_disp_pars1)
        self.w.negCol.currentIndexChanged.connect(self.get_disp_pars2)
        self.w.posColR.textChanged.connect(self.get_disp_pars3)
        self.w.posColG.textChanged.connect(self.get_disp_pars3)
        self.w.posColB.textChanged.connect(self.get_disp_pars3)
        self.w.negColR.textChanged.connect(self.get_disp_pars3)
        self.w.negColG.textChanged.connect(self.get_disp_pars3)
        self.w.negColB.textChanged.connect(self.get_disp_pars3)
        self.w.nLevels.textChanged.connect(self.get_disp_pars4)
        self.w.minLevel.textChanged.connect(self.get_disp_pars5)
        self.w.maxLevel.textChanged.connect(self.get_disp_pars6)
        self.w.axisType1.currentIndexChanged.connect(self.get_disp_pars7)
        self.w.axisType2.currentIndexChanged.connect(self.get_disp_pars8)
        self.w.displaySpc.currentIndexChanged.connect(self.get_disp_pars9)
        self.w.baselineCorrection.currentIndexChanged.connect(self.check_baseline_correction)
        self.w.baselineOrder.currentIndexChanged.connect(self.check_baseline_order)
        self.w.spcOffset.textChanged.connect(self.get_disp_pars10)
        self.w.spcScale.textChanged.connect(self.get_disp_pars11)
        self.w.fontSize.valueChanged.connect(self.set_font_size)
        self.w.xLabel.textChanged.connect(self.get_disp_pars12)
        self.w.yLabel.textChanged.connect(self.get_disp_pars13)
        self.w.spcLabel.textChanged.connect(self.get_disp_pars14)
        self.w.preProcessingSelect.currentIndexChanged.connect(self.set_pre_processing_options)
        self.w.exportMethod.currentIndexChanged.connect(self.set_export_method_options)
        self.w.tilt.currentIndexChanged.connect(self.set_tilt)
        self.w.symJ.currentIndexChanged.connect(self.set_sym_j)
        self.w.windowFunction.currentIndexChanged.connect(self.get_proc_pars1)
        self.w.windowFunction_2.currentIndexChanged.connect(self.get_proc_pars2)
        self.w.phaseCorrection.currentIndexChanged.connect(self.get_proc_pars3)
        self.w.phaseCorrection_2.currentIndexChanged.connect(self.get_proc_pars4)
        self.w.waterSuppression.currentIndexChanged.connect(self.get_proc_pars5)
        self.w.winType.currentIndexChanged.connect(self.get_proc_pars6)
        self.w.gibbs.currentIndexChanged.connect(self.get_proc_pars7)
        self.w.gibbs_2.currentIndexChanged.connect(self.get_proc_pars8)
        self.w.zeroFilling.textChanged.connect(self.get_proc_pars9)
        self.w.zeroFilling_2.textChanged.connect(self.get_proc_pars10)
        self.w.lb.textChanged.connect(self.get_proc_pars11)
        self.w.gb.textChanged.connect(self.get_proc_pars12)
        self.w.ssb.textChanged.connect(self.get_proc_pars13)
        self.w.lb_2.textChanged.connect(self.get_proc_pars14)
        self.w.gb_2.textChanged.connect(self.get_proc_pars15)
        self.w.ssb_2.textChanged.connect(self.get_proc_pars16)
        self.w.ph0.textChanged.connect(self.get_proc_pars17)
        self.w.ph1.textChanged.connect(self.get_proc_pars18)
        self.w.ph0_2.textChanged.connect(self.get_proc_pars19)
        self.w.ph1_2.textChanged.connect(self.get_proc_pars20)
        self.w.polyOrder.textChanged.connect(self.get_proc_pars21)
        self.w.extrapolationSize.textChanged.connect(self.get_proc_pars22)
        self.w.windowSize.textChanged.connect(self.get_proc_pars23)
        self.w.fidOffsetCorrection.textChanged.connect(self.get_proc_pars24)
        self.w.stripTransformStart.textChanged.connect(self.get_proc_pars25)
        self.w.stripTransformEnd.textChanged.connect(self.get_proc_pars26)
        self.w.phRefDS.valueChanged.connect(self.change_data_set_exp_ph_ref)
        self.w.phRefExp.valueChanged.connect(self.change_data_set_exp_ph_ref)
        self.w.phRefColour.currentIndexChanged.connect(self.get_disp_pars15)
        self.w.fourierTransformButton.clicked.connect(self.ft)
        self.w.executeScript.clicked.connect(self.exec_script)
        self.w.openScript.clicked.connect(self.open_script)
        self.w.saveScript.clicked.connect(self.save_script)
        self.w.actionOpen_Script.triggered.connect(self.open_script)
        self.w.actionSave_Script.triggered.connect(self.save_script)
        self.w.actionExecute_Script.triggered.connect(self.exec_script)
        # self.w.helpComboBox.currentIndexChanged.connect(self.set_help)
        self.w.helpComboBox.activated.connect(self.set_help)
        # Quit Button
        self.w.quitButton.clicked.connect(self.quit_app)
        self.w.saveButton.clicked.connect(self.save_button)
        self.w.loadButton.clicked.connect(self.load_button)
        self.w.exportPathSelectButton.clicked.connect(self.set_export_table)
        self.w.actionQuit.triggered.connect(self.quit_app)
        self.w.dispPlotButton.clicked.connect(self.plot_spc_disp)
        self.show_version()
        self.keep_zoom = False
        self.keep_x_zoom = False
        self.ph_corr_active = False
        self.set_font_size()
        self.cf = nmrConfig.NmrConfig()
        self.cf.read_config()
        self.w.autoPlot.setChecked(self.cf.auto_plot)
        self.w.keepZoom.setChecked(self.cf.keep_zoom)
        self.w.fontSize.setValue(self.cf.font_size)
        self.std_pos_col1 = (self.cf.pos_col10, self.cf.pos_col11, self.cf.pos_col12)
        self.std_neg_col1 = (self.cf.neg_col10, self.cf.neg_col11, self.cf.neg_col12)
        self.std_pos_col2 = (self.cf.pos_col20, self.cf.pos_col21, self.cf.pos_col22)
        self.std_neg_col2 = (self.cf.neg_col20, self.cf.neg_col21, self.cf.neg_col22)
        self.w.actionSave_as_Default.triggered.connect(self.save_config)
        self.w.actionLoad_Default.triggered.connect(self.load_config)
        self.w.actionReset_Config.triggered.connect(self.reset_config)
        self.w.rSpc_p0.textChanged.connect(self.get_r_spc_p0)
        self.w.rSpc_p1.textChanged.connect(self.get_r_spc_p1)
        self.w.rSpc_p2.textChanged.connect(self.get_r_spc_p2)
        self.w.rSpc_p3.textChanged.connect(self.get_r_spc_p3)
        self.w.rSpc_p4.textChanged.connect(self.get_r_spc_p4)
        self.w.rSpc_p5.textChanged.connect(self.get_r_spc_p5)
        self.w.rSpc_p6.textChanged.connect(self.get_r_spc_p6)
        self.w.iSpc_p0.textChanged.connect(self.get_i_spc_p0)
        self.w.iSpc_p1.textChanged.connect(self.get_i_spc_p1)
        self.w.iSpc_p2.textChanged.connect(self.get_i_spc_p2)
        self.w.iSpc_p3.textChanged.connect(self.get_i_spc_p3)
        self.w.iSpc_p4.textChanged.connect(self.get_i_spc_p4)
        self.w.iSpc_p5.textChanged.connect(self.get_i_spc_p5)
        self.w.iSpc_p6.textChanged.connect(self.get_i_spc_p6)
        self.set_font_size()
        self.w.MplWidget.toolbar.setVisible(False)
        self.w.MplWidget2.toolbar.setVisible(False)
        self.w.MplWidget3.toolbar.setVisible(False)
        self.w.MplWidget.setFocus()
        self.set_zoom()
        self.w.pickRowColPhCorr2d.clicked.connect(self.pick_col_row)
        self.w.emptyRowColPhCorr2d.clicked.connect(self.empty_col_row)
        self.w.removeRowColPhCorr2d.clicked.connect(self.remove_last_col_row)
        self.w.horzPhCorr2d.clicked.connect(self.horz_ph_corr_2d)
        self.w.vertPhCorr2d.clicked.connect(self.vert_ph_corr_2d)
        self.w.exitPhCorr2d.clicked.connect(self.start_stop_ph_corr)
        self.w.applyPhCorr2d.clicked.connect(self.apply_2d_ph_corr)
        self.w.cancelPhCorr2d.clicked.connect(self.cancel_2d_ph_corr)
        self.w.zoomPhCorr2d.clicked.connect(self.zoom_ph_corr)
        self.w.exitZoomPhCorr2d.clicked.connect(self.zoom_ph_corr)
        self.w.exitPhCorr1d.setVisible(False)
        self.w.zoomPhCorr1d.setVisible(False)
        self.w.exitZoomPhCorr1d.setVisible(False)
        self.w.pickRowColPhCorr2d.setVisible(False)
        self.w.emptyRowColPhCorr2d.setVisible(False)
        self.w.removeRowColPhCorr2d.setVisible(False)
        self.w.horzPhCorr2d.setVisible(False)
        self.w.vertPhCorr2d.setVisible(False)
        self.w.zoomPhCorr2d.setVisible(False)
        self.w.applyPhCorr2d.setVisible(False)
        self.w.cancelPhCorr2d.setVisible(False)
        self.w.exitPhCorr2d.setVisible(False)
        self.w.exitZoomPhCorr2d.setVisible(False)
        self.w.actionSet_light_mode_requires_restart.triggered.connect(self.set_light_mode)
        self.w.actionSet_dark_mode_requires_restart.triggered.connect(self.set_dark_mode)
        self.w.MplWidget.canvas.draw()
        self.w.setStyleSheet("font-size: " + str(self.cf.font_size) + "pt")
        self.w.actionreInitialise_pre_processing_plot_colours.triggered.connect(self.nd.pp.init_plot_colours)
        self.w.actionreInitialise_plot_colours.triggered.connect(self.set_standard_plot_colours)
        # print(sys.platform)
        if sys.platform == 'darwin':
            self.w.actionCreate.setText('Create Launchpad Icon')
            self.w.actionCreate.triggered.connect(self.create_icon_mac)
        elif sys.platform == 'win' or sys.platform == 'win32' or sys.platform == 'win64':
            self.w.actionCreate.setText('Create Desktop Icon')
            self.w.actionCreate.triggered.connect(self.create_icon_win)
        else:
            # print(sys.platform)
            self.w.actionCreate.setText('Create Desktop Starter')
            # print('doing stuffs3....')

            self.w.actionCreate.setVisible(False)

        self.emp_ref_shift = 0.0
        self.p = []
        if self.cf.mode == 'dark':
            self.load_dark_mode()
        else:
            self.load_light_mode()
        #
        self.w.helpView.page().profile().downloadRequested.connect(self._download_requested)
        self.w.peakWidget.setColumnWidth(2, 182)
        # end __init__

    def activate_command_line(self):
        if (self.w.cmdLine.hasFocus() == True):
            self.w.cmdLine.clearFocus()
        else:
            self.w.cmdLine.setFocus()

        # end activate_command_line

    def add_peak(self):
        self.ginput_add_peak(2)
        # end add_peak

    def clear_peak(self):
        self.nd.clear_peak()
        self.w.peakWidget.setRowCount(0)
        self.plot_spc()
        # end add_peak

    def apply_2d_ph_corr(self):
        s = self.nd.s
        e = self.nd.e
        if self.nd.nmrdat[s][e].proc.phase_inversion:
            self.ph_corr.ph0_2d[self.ph_corr.dim] *= -1
            self.ph_corr.ph1_2d[self.ph_corr.dim] *= -1

        cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.on_ph_corr_click_2d)
        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.on_ph_corr_release_2d)
        cid = self.w.MplWidget.canvas.mpl_disconnect(cid)
        cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
        # self.w.actionApplyPhCorr.triggered.disconnect()
        # self.w.actionCancelPhCorr.triggered.disconnect()
        self.w.pickRowColPhCorr2d.setVisible(True)
        self.w.emptyRowColPhCorr2d.setVisible(True)
        self.w.removeRowColPhCorr2d.setVisible(True)
        self.w.horzPhCorr2d.setVisible(True)
        self.w.vertPhCorr2d.setVisible(True)
        self.w.zoomPhCorr2d.setVisible(False)
        self.w.applyPhCorr2d.setVisible(False)
        self.w.cancelPhCorr2d.setVisible(False)
        self.w.exitPhCorr2d.setVisible(True)
        self.w.exitZoomPhCorr2d.setVisible(False)
        ph0 = ((self.ph_corr.ph0_2d[self.ph_corr.dim] + 180.0) % 360.0) - 180.0
        ph1 = self.ph_corr.ph1_2d[self.ph_corr.dim]
        if self.nd.nmrdat[s][e].proc.phase_inversion is False:
            self.nd.nmrdat[s][e].phase2a(ph0, ph1, self.ph_corr.dim)
        else:
            self.nd.nmrdat[s][e].phase2a(-ph0, -ph1, self.ph_corr.dim)

        ph0 = ((ph0 + self.nd.nmrdat[s][e].proc.ph0[self.ph_corr.dim] + 180.0) % 360.0) - 180.0
        ph1 = ph1 + self.nd.nmrdat[s][e].proc.ph1[self.ph_corr.dim]

        self.nd.nmrdat[s][e].proc.ph0[self.ph_corr.dim] = ph0
        self.nd.nmrdat[s][e].proc.ph1[self.ph_corr.dim] = ph1
        self.ph_corr.ph0_2d[self.ph_corr.dim] = 0
        self.ph_corr.ph1_2d[self.ph_corr.dim] = 0
        self.ph_corr.spc = np.array([[]], dtype='complex')
        self.ph_corr.spc2 = np.array([[]], dtype='complex')
        zoom_status = self.w.keepZoom.isChecked()
        self.w.keepZoom.setChecked(False)
        self.plot_spc()
        self.w.keepZoom.setChecked(zoom_status)
        self.plot_2d_col_row()
        if (self.zoom_was_on == True):
            self.set_zoom_off()
            self.set_zoom()

        if (self.pan_was_on == True):
            self.set_pan()

        self.show_ph_corr2d()
        self.set_proc_pars()
        self.show_acquisition_parameters()
        self.show_nmr_spectrum()
        # end apply_2d_ph_corr

    def autobaseline1d(self):
        code_out = io.StringIO()
        code_err = io.StringIO()
        sys.stdout = code_out
        sys.stderr = code_err
        self.show_auto_baseline()
        self.nd.ft()
        self.nd.auto_ref()
        self.nd.autobaseline1d()
        self.w.baselineCorrection.setCurrentIndex(1)
        self.nd.ft()
        self.nd.baseline1d()
        # self.w.baselineCorrection.setCurrentIndex(1)
        self.set_proc_pars()
        self.show_version()
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.change_data_set_exp()
        self.plot_spc()
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        if self.cf.mode == 'dark':
            txt_col = QColor.fromRgbF(1.0, 1.0, 1.0, 1.0)
            err_col = QColor.fromRgbF(1.0, 0.5, 0.5, 1.0)
        else:
            txt_col = QColor.fromRgbF(0.0, 0.0, 0.0, 1.0)
            err_col = QColor.fromRgbF(1.0, 0.0, 0.0, 1.0)

        self.w.console.setTextColor(txt_col)
        self.w.console.append(code_out.getvalue())
        self.w.console.setTextColor(err_col)
        self.w.console.append(code_err.getvalue())
        code_out.close()
        code_err.close()
        # end autobaseline1d

    def autobaseline1d_all(self):
        code_out = io.StringIO()
        code_err = io.StringIO()
        sys.stdout = code_out
        sys.stderr = code_err
        self.show_auto_baseline()
        self.nd.ft()
        self.nd.auto_ref()
        self.nd.autobaseline1d_all()
        self.w.baselineCorrection.setCurrentIndex(1)
        self.nd.ft()
        self.nd.baseline1d()
        # self.w.baselineCorrection.setCurrentIndex(1)
        self.set_proc_pars()
        self.show_version()
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.change_data_set_exp()
        self.plot_spc()
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        if self.cf.mode == 'dark':
            txt_col = QColor.fromRgbF(1.0, 1.0, 1.0, 1.0)
            err_col = QColor.fromRgbF(1.0, 0.5, 0.5, 1.0)
        else:
            txt_col = QColor.fromRgbF(0.0, 0.0, 0.0, 1.0)
            err_col = QColor.fromRgbF(1.0, 0.0, 0.0, 1.0)

        self.w.console.setTextColor(txt_col)
        self.w.console.append(code_out.getvalue())
        self.w.console.setTextColor(err_col)
        self.w.console.append(code_err.getvalue())
        code_out.close()
        code_err.close()
        # end autobaseline1d_all

    def autophase1d(self):
        code_out = io.StringIO()
        code_err = io.StringIO()
        sys.stdout = code_out
        sys.stderr = code_err
        self.show_auto_phase()
        self.nd.ft()
        self.nd.auto_ref()
        self.nd.autophase1d()
        self.w.baselineCorrection.setCurrentIndex(1)
        self.nd.ft()
        self.nd.baseline1d()
        self.plot_spc()
        self.set_proc_pars()
        self.show_version()
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.change_data_set_exp()
        self.plot_spc()
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        if self.cf.mode == 'dark':
            txt_col = QColor.fromRgbF(1.0, 1.0, 1.0, 1.0)
            err_col = QColor.fromRgbF(1.0, 0.5, 0.5, 1.0)
        else:
            txt_col = QColor.fromRgbF(0.0, 0.0, 0.0, 1.0)
            err_col = QColor.fromRgbF(1.0, 0.0, 0.0, 1.0)

        self.w.console.setTextColor(txt_col)
        self.w.console.append(code_out.getvalue())
        self.w.console.setTextColor(err_col)
        self.w.console.append(code_err.getvalue())
        code_out.close()
        code_err.close()
        # end autophase1d

    def autophase1d_all(self):
        code_out = io.StringIO()
        code_err = io.StringIO()
        sys.stdout = code_out
        sys.stderr = code_err
        self.show_auto_phase()
        self.nd.ft()
        self.nd.auto_ref()
        self.nd.autophase1d_all()
        self.w.baselineCorrection.setCurrentIndex(1)
        self.nd.ft()
        self.nd.baseline1d()
        self.plot_spc()
        self.set_proc_pars()
        self.show_version()
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.change_data_set_exp()
        self.plot_spc()
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        if self.cf.mode == 'dark':
            txt_col = QColor.fromRgbF(1.0, 1.0, 1.0, 1.0)
            err_col = QColor.fromRgbF(1.0, 0.5, 0.5, 1.0)
        else:
            txt_col = QColor.fromRgbF(0.0, 0.0, 0.0, 1.0)
            err_col = QColor.fromRgbF(1.0, 0.0, 0.0, 1.0)

        self.w.console.setTextColor(txt_col)
        self.w.console.append(code_out.getvalue())
        self.w.console.setTextColor(err_col)
        self.w.console.append(code_err.getvalue())
        code_out.close()
        code_err.close()
        # end autophase1d_all

    def autoref(self, tmsp=True):
        self.nd.auto_ref(tmsp)
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.change_data_set_exp()
        self.plot_spc()
        return "autoref"
        # end autoref

    def baseline1d(self):
        self.nd.baseline1d()
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.change_data_set_exp()
        self.plot_spc()
        # end baseline1d

    def cancel_2d_ph_corr(self):
        cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.on_ph_corr_click_2d)
        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.on_ph_corr_release_2d)
        cid = self.w.MplWidget.canvas.mpl_disconnect(cid)
        cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
        # self.w.actionApplyPhCorr.triggered.disconnect()
        # self.w.actionCancelPhCorr.triggered.disconnect()
        self.w.pickRowColPhCorr2d.setVisible(True)
        self.w.emptyRowColPhCorr2d.setVisible(True)
        self.w.removeRowColPhCorr2d.setVisible(True)
        self.w.horzPhCorr2d.setVisible(True)
        self.w.vertPhCorr2d.setVisible(True)
        self.w.zoomPhCorr2d.setVisible(False)
        self.w.applyPhCorr2d.setVisible(False)
        self.w.cancelPhCorr2d.setVisible(False)
        self.w.exitPhCorr2d.setVisible(True)
        self.w.exitZoomPhCorr2d.setVisible(False)
        self.ph_corr.ph0_2d[self.ph_corr.dim] = 0
        self.ph_corr.ph1_2d[self.ph_corr.dim] = 0
        zoomStatus = self.w.keepZoom.isChecked()
        self.w.keepZoom.setChecked(False)
        self.plot_spc()
        self.w.keepZoom.setChecked(zoomStatus)
        self.plot_2d_col_row()
        if (self.zoom_was_on == True):
            self.set_zoom_off()
            self.set_zoom()

        if (self.pan_was_on == True):
            self.set_pan()

        self.show_ph_corr2d()
        self.show_acquisition_parameters()
        self.show_nmr_spectrum()
        # end cancel2dPhCorr

    def change_data_set_exp(self):
        cidx = self.w.nmrSpectrum.currentIndex()
        if (len(self.nd.nmrdat) > 0):
            if (len(self.nd.nmrdat[self.nd.s]) > 0):
                self.keep_zoom = self.w.keepZoom.isChecked()
                old_set = self.nd.s
                old_exp = self.nd.e
                if (self.w.setBox.value() < 1):
                    self.w.setBox.setValue(1)

                if (self.w.expBox.value() < 1):
                    self.w.expBox.setValue(1)

                if (self.w.setBox.value() > len(self.nd.nmrdat)):
                    self.w.setBox.setValue(len(self.nd.nmrdat))

                self.nd.s = self.w.setBox.value() - 1
                if (self.w.expBox.value() > len(self.nd.nmrdat[self.nd.s])):
                    self.w.expBox.setValue(len(self.nd.nmrdat[self.nd.s]))

                self.nd.e = self.w.expBox.value() - 1
                keep_zoom = self.w.keepZoom.isChecked()
                if not ((old_set == self.nd.s) and (old_exp == self.nd.e)):
                    if (self.nd.nmrdat[old_set][old_exp].dim != self.nd.nmrdat[self.nd.s][self.nd.e].dim):
                        self.keep_x_zoom = True
                        self.keep_zoom = False
                        self.w.keepZoom.setChecked(False)

                    self.set_disp_pars()
                    self.set_proc_pars()
                    self.set_acq_pars()
                    self.set_title_file()
                    self.set_pulse_program()
                    if (self.ph_corr_active == False):
                        if (self.w.autoPlot.isChecked()):
                            self.plot_spc()
                        elif (self.w.nmrSpectrum.currentIndex() == 0):
                            self.plot_spc()

                    else:
                        if self.nd.nmrdat[self.nd.s][self.nd.e].dim == 1:
                            self.ph_corr.spc = self.nd.nmrdat[self.nd.s][self.nd.e].spc
                            self.ph_corr_plot_spc()
                        else:
                            self.plot_spc()
                            self.plot_2d_col_row()

                    self.keep_zoom = keep_zoom
                    self.w.keepZoom.setChecked(keep_zoom)
                # else:
                #    if (self.ph_corr_active == False):
                #        if (self.w.autoPlot.isChecked()):
                #            self.plot_spc()
                #        elif (self.w.nmrSpectrum.currentIndex() == 0):
                #            self.plot_spc()
                #
                #    else:
                #        self.ph_corr.spc = self.nd.nmrdat[self.nd.s][self.nd.e].spc
                #        self.ph_corr_plot_spc()

                self.keep_zoom = False

            else:
                self.w.setBox.valueChanged.disconnect()
                self.w.expBox.valueChanged.disconnect()
                self.w.expBox.setValue(0)
                self.w.setBox.setValue(0)
                self.w.setBox.valueChanged.connect(lambda: self.change_data_set_exp())
                self.w.expBox.valueChanged.connect(lambda: self.change_data_set_exp())

            self.update_gui()
        else:
            self.w.setBox.valueChanged.disconnect()
            self.w.expBox.valueChanged.disconnect()
            self.w.expBox.setValue(0)
            self.w.setBox.setValue(0)
            self.w.setBox.valueChanged.connect(lambda: self.change_data_set_exp())
            self.w.expBox.valueChanged.connect(lambda: self.change_data_set_exp())

        if (self.w.autoPlot.isChecked() is False):
            self.w.nmrSpectrum.setCurrentIndex(cidx)

        # end change_data_set_exp

    def change_data_set_exp_ph_ref(self):
        if (len(self.nd.nmrdat) > 0):
            s = self.nd.s
            e = self.nd.e
            if (len(self.nd.nmrdat[self.nd.s]) > 0):
                if (self.w.phRefDS.value() < 0):
                    self.w.phRefDS.setValue(0)

                if (self.w.phRefExp.value() < 0):
                    self.w.phRefExp.setValue(0)

                if (self.w.phRefDS.value() > len(self.nd.nmrdat)):
                    self.w.phRefExp.setValue(len(self.nd.nmrdat))

                if (self.w.expBox.value() > len(self.nd.nmrdat[self.nd.s])):
                    self.w.expBox.setValue(len(self.nd.nmrdat[self.nd.s]))

                for k in range(len(self.nd.nmrdat)):
                    for l in range(len(self.nd.nmrdat[k])):
                        self.nd.nmrdat[k][l].display.ph_ref_ds = self.w.phRefDS.value()
                        self.nd.nmrdat[k][l].display.ph_ref_exp = self.w.phRefExp.value()

        # end change_data_set_exp_ph_ref

    def change_scale_spectra_ref_spc(self):
        self.nd.pp.scale_spectra_ref_spc = self.w.scaleSpectraRefSpc.value()
        # end change_scale_spectra_ref_spc

    def change_seg_align_ref_spc(self):
        self.nd.pp.seg_align_ref_spc = self.w.segAlignRefSpc.value()
        # end change_seg_align_ref_spc

    def change_standard_colours(self, pos_col1=(), neg_col1=(), pos_col2=(), neg_col2=()):
        if len(pos_col1) != 3:
            pos_col1 = self.std_pos_col1

        if len(neg_col1) != 3:
            neg_col1 = self.std_neg_col1

        if len(pos_col2) != 3:
            pos_col2 = self.std_pos_col2

        if len(neg_col2) != 3:
            neg_col2 = self.std_neg_col2

        self.std_pos_col1 = pos_col1
        self.std_neg_col1 = neg_col1
        self.std_pos_col2 = pos_col2
        self.std_neg_col2 = neg_col2
        self.set_standard_colours()

    # end change_standard_colours

    def change_to_next_ds(self):
        self.w.setBox.setValue(self.w.setBox.value() + 1)
        # end change_to_next_ds

    def change_to_next_exp(self):
        self.w.expBox.setValue(self.w.expBox.value() + 1)
        # end change_to_next_exp

    def change_to_previous_ds(self):
        self.w.setBox.setValue(self.w.setBox.value() - 1)
        # end change_to_previous_ds

    def change_to_previous_exp(self):
        self.w.expBox.setValue(self.w.expBox.value() - 1)
        # end change_to_previous_exp

    def check_baseline_correction(self):
        cbl = self.w.baselineCorrection.currentIndex()
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.correct_baseline = cbl
        if (cbl == 1):
            self.w.baselineOrder.setEnabled(True)
        else:
            self.w.baselineOrder.setEnabled(False)

        self.check_baseline_order()
        # end check_baseline_correction

    def check_baseline_order(self):
        blo = self.w.baselineOrder.currentIndex()
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.n_order = blo
        self.w.rSpc_p0.setEnabled(False)
        self.w.rSpc_p1.setEnabled(False)
        self.w.rSpc_p2.setEnabled(False)
        self.w.rSpc_p3.setEnabled(False)
        self.w.rSpc_p4.setEnabled(False)
        self.w.rSpc_p5.setEnabled(False)
        self.w.rSpc_p6.setEnabled(False)
        self.w.iSpc_p0.setEnabled(False)
        self.w.iSpc_p1.setEnabled(False)
        self.w.iSpc_p2.setEnabled(False)
        self.w.iSpc_p3.setEnabled(False)
        self.w.iSpc_p4.setEnabled(False)
        self.w.iSpc_p5.setEnabled(False)
        self.w.iSpc_p6.setEnabled(False)
        if (self.w.baselineOrder.isEnabled() == True):
            self.w.rSpc_p0.setEnabled(True)
            self.w.iSpc_p0.setEnabled(True)
            if (blo > 0):
                self.w.rSpc_p1.setEnabled(True)
                self.w.iSpc_p1.setEnabled(True)

            if (blo > 1):
                self.w.rSpc_p2.setEnabled(True)
                self.w.iSpc_p2.setEnabled(True)

            if (blo > 2):
                self.w.rSpc_p3.setEnabled(True)
                self.w.iSpc_p3.setEnabled(True)

            if (blo > 3):
                self.w.rSpc_p4.setEnabled(True)
                self.w.iSpc_p4.setEnabled(True)

            if (blo > 4):
                self.w.rSpc_p5.setEnabled(True)
                self.w.iSpc_p5.setEnabled(True)

            if (blo > 5):
                self.w.rSpc_p6.setEnabled(True)
                self.w.iSpc_p6.setEnabled(True)

        # end check_baseline_order

    def clear(self):
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        self.w.MplWidget.canvas.axes.clear()
        self.w.MplWidget.canvas.draw()
        self.zero_disp_pars()
        self.zero_proc_pars()
        self.zero_acq_pars()
        self.zero_console()
        self.zero_title_file()
        self.zero_pulse_program()
        self.nd.nmrdat = [[]]
        self.nd.s = 0
        self.nd.e = -1
        self.w.setBox.valueChanged.disconnect()
        self.w.expBox.valueChanged.disconnect()
        self.w.expBox.setValue(0)
        self.w.setBox.setValue(0)
        self.w.setBox.valueChanged.connect(lambda: self.change_data_set_exp())
        self.w.expBox.valueChanged.connect(lambda: self.change_data_set_exp())
        code_out = io.StringIO()
        code_err = io.StringIO()
        sys.stdout = code_out
        sys.stderr = code_err
        return "Workspace cleared"
        # end clear

    def create_icon_mac(self):
        home_dir = os.path.expanduser('~')
        app_dir2 = os.path.join(home_dir, 'Applications')
        app_dir = os.path.join(app_dir2, 'MetaboLabPy.app')
        app_dir1 = 'MetaboLabPy'
        try:
            shutil.rmtree(app_dir)
        except:
            pass

        appify = os.path.join(app_dir2, 'appify')
        f = open(appify, 'w')
        f.write('#!/usr/bin/env bash\n\n')
        f.write('APp_name=${2:-$(basename "$1" ".sh")}\n')
        f.write('DIR="$APp_name/$3.app/Contents/MacOS"\n\n')
        f.write('mkdir -p "$DIR"\n')
        f.write('cp "$1" "$DIR/$3"\n')
        f.write('chmod +x "$DIR/$3"\n')
        f.close()
        os.chmod(appify, 0o777)
        base_dir = os.path.dirname(__file__)
        ml_starter = os.path.join(base_dir, 'mlStarter')
        contents = os.path.join(ml_starter, 'Contents')
        icon = os.path.join(ml_starter, 'Icon')
        starter = os.path.join(app_dir2, 'createStarter')
        f = open(starter, 'w')
        f.write('#!/usr/bin/env bash\n\n')
        f.write(appify + ' $(which metabolabpy) ' + app_dir2 + ' ' + app_dir1 + '\n')
        f.write('cp -r ' + contents + ' ' + app_dir + '\n')
        f.write("cp " + icon + " " + app_dir + "/Icon$'\\r'\n")
        f.close()
        os.chmod(starter, 0o777)
        subprocess.os.system(starter)
        os.remove(appify)
        os.remove(starter)
        # end create_icon_mac

    def create_icon_linux(self):
        print("Linux!")
        # end create_icon_linux

    def create_icon_win(self):
        base_dir = os.path.dirname(__file__)
        icon_file = os.path.join(base_dir, 'icon', 'icon.ico')
        user_dir = os.environ.get('USERPROFILE')
        desktop_dir = os.path.join(user_dir, 'Desktop')
        link_file = os.path.join(desktop_dir, 'MetaboLabPy.lnk')
        ml_bat = os.path.join(base_dir, 'ml.bat')
        ml_exec_bat = os.path.join(base_dir, 'ml_exec.bat')
        f = open(ml_bat, 'w')
        f.write('start /min ' + ml_exec_bat)
        f.close()
        f = open(ml_exec_bat, 'w')
        venv = sys.prefix.find('env')
        cnda = subprocess.check_output('where conda').decode()
        if venv == -1:
            f.write('metabolabpy && exit')
        else:
            idx = sys.prefix.rfind('\\') + 1
            env = sys.prefix[idx:]
            f.write(cnda[:-2] + ' activate ' + env + ' && metabolabpy && exit')

        f.close()
        subprocess.os.system('pip install pylnk3')
        subprocess.os.system('pylnk3 create ' + ml_bat + ' ' + link_file + ' -m Minimized --icon ' + icon_file)
        subprocess.os.system('pip uninstall pylnk3 --yes')
        # end create_icon_win

    def data_pre_processing(self):
        self.reset_data_pre_processing()
        self.nd.data_pre_processing()
        self.plot_spc_pre_proc()
        self.vertical_auto_scale()
        self.w.MplWidget.canvas.flush_events()
        self.w.MplWidget.canvas.draw()
        # end data_pre_processing

    def empty_col_row(self):
        while len(self.w.MplWidget.canvas.axes.lines) > 0:
            self.w.MplWidget.canvas.axes.lines[0].remove()

        self.w.MplWidget.canvas.draw()
        self.ph_corr.spc_row = []
        self.ph_corr.spc_col = []
        self.ph_corr.spc_row_pts = []
        self.ph_corr.spc_col_pts = []
        self.show_acquisition_parameters()
        self.show_nmr_spectrum()
        # end empty_col_row

    def enable_baseline(self):
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            self.nd.nmrdat[self.nd.s][k].apc.correct_baseline = 1

        self.w.baselineOrder.setCurrentIndex(self.nd.nmrdat[self.nd.s][0].apc.n_order)
        self.w.baselineCorrection.setCurrentIndex(1)
        return "baselineCorrection enabled"
        # end enableBaseline

    def exec_cmd(self):
        if self.cf.mode == 'dark':
            txt_col = QColor.fromRgbF(1.0, 1.0, 1.0, 1.0)
            err_col = QColor.fromRgbF(1.0, 0.5, 0.5, 1.0)
        else:
            txt_col = QColor.fromRgbF(0.0, 0.0, 0.0, 1.0)
            err_col = QColor.fromRgbF(1.0, 0.0, 0.0, 1.0)

        cmd_text = self.w.cmdLine.text()
        if (len(cmd_text) > 0):
            self.w.nmrSpectrum.setCurrentIndex(11)
            self.w.cmdLine.setText("")
            self.nd.cmd_buffer = np.append(self.nd.cmd_buffer, cmd_text)
            self.nd.cmd_idx = len(self.nd.cmd_buffer)
            code_out = io.StringIO()
            code_err = io.StringIO()
            sys.stdout = code_out
            sys.stderr = code_err
            print(">>> " + cmd_text)
            try:
                output = eval(cmd_text)
                print(output)
                self.w.console.setTextColor(txt_col)
                self.w.console.append(code_out.getvalue())
            except:  # (SyntaxError, NameError, TypeError, ZeroDivisionError, AttributeError, ArithmeticError, BufferError, LookupError):
                cmd_text2 = "self." + cmd_text
                try:
                    output = eval(cmd_text2)
                    print(output)
                    self.w.console.setTextColor(txt_col)
                    self.w.console.append(code_out.getvalue())
                except:
                    traceback.print_exc()
                    self.w.console.setTextColor(err_col)
                    self.w.console.append(code_out.getvalue())
                    self.w.console.append(code_err.getvalue())

            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__
            code_out.close()
            code_err.close()
            self.w.console.verticalScrollBar().setValue(self.w.console.verticalScrollBar().maximum())

        # end execCmd

    def exec_script(self):
        if self.cf.mode == 'dark':
            txt_col = QColor.fromRgbF(1.0, 1.0, 1.0, 1.0)
            err_col = QColor.fromRgbF(1.0, 0.5, 0.5, 1.0)
            scr_col = QColor.fromRgbF(0.5, 0.5, 1.0, 1.0)
            scr_col2 = QColor.fromRgbF(0.4, 0.4, 1.0, 1.0)
            adm_col = QColor.fromRgbF(1.0, 1.0, 0.5, 1.0)
        else:
            txt_col = QColor.fromRgbF(0.0, 0.0, 0.0, 1.0)
            err_col = QColor.fromRgbF(1.0, 0.0, 0.0, 1.0)
            scr_col = QColor.fromRgbF(0.0, 0.0, 1.0, 1.0)
            scr_col2 = QColor.fromRgbF(0.0, 0.0, 0.6, 1.0)
            adm_col = QColor.fromRgbF(0.4, 0.4, 0.4, 1.0)

        zoom_checked = self.w.keepZoom.isChecked()
        self.w.keepZoom.setChecked(False)
        code_out = io.StringIO()
        code_err = io.StringIO()
        sys.stdout = code_out
        sys.stderr = code_err
        code = self.w.script.toPlainText()
        code = code.replace('\\', '\\' * 2)
        try:
            exec(code)
            code = self.w.script.toPlainText()
            code = code.replace(' interactive ', ' abcint2 ')
            code = code.replace('''interactive''', '' + self.nd.nmrdat[self.nd.s][self.nd.e].data_set_name + '')
            code = code.replace(' abcint2 ', ' interactive ')
            self.w.script.setText(code)


        except:  # (SyntaxError, NameError, TypeError, ZeroDivisionError, AttributeError):
            self.w.nmrSpectrum.setCurrentIndex(11)
            traceback.print_exc()

        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        self.w.console.setTextColor(adm_col)
        self.w.console.append('--- Script_start -------------------------\n')
        self.w.console.setTextColor(scr_col2)
        self.w.console.append('Executing script...\n')
        self.w.console.setTextColor(scr_col)
        codeSplit = code.split('\n')
        for k in range(len(codeSplit)):
            self.w.console.append(str(k + 1) + ': ' + str(codeSplit[k]))

        self.w.console.setTextColor(adm_col)
        self.w.console.append('\n--- ScriptOutput ------------------------\n')
        self.w.console.setTextColor(txt_col)
        self.w.console.append(code_out.getvalue())
        self.w.console.setTextColor(err_col)
        self.w.console.append(code_err.getvalue())
        self.w.console.setTextColor(adm_col)
        self.w.console.append('--- Script_end ---------------------------\n')
        self.w.console.setTextColor(txt_col)
        code_out.close()
        code_err.close()
        if (len(self.nd.nmrdat[0]) > 0):
            self.update_gui()

        self.w.setBox.valueChanged.disconnect()
        self.w.expBox.valueChanged.disconnect()
        self.w.expBox.setValue(self.nd.e + 1)
        self.w.setBox.setValue(self.nd.s + 1)
        self.w.setBox.valueChanged.connect(lambda: self.change_data_set_exp())
        self.w.expBox.valueChanged.connect(lambda: self.change_data_set_exp())
        if (zoom_checked == True):
            self.w.keepZoom.setChecked(True)

        # end exec_script

    def export_peak(self):
        f_name = QFileDialog.getSaveFileName(None, "Save Excel file", "", "*.xlsx", "*.xlsx")
        f_name = f_name[0]
        if len(f_name) == 0:
            return

        workbook = xlsxwriter.Workbook(f_name)
        if self.nd.int_all_data_sets == True:
            ds = range(len(self.nd.nmrdat))
        else:
            ds = [self.nd.s]

        for k in range(len(ds)):
            worksheet = workbook.add_worksheet('Dataset ' + str(ds[k] + 1))
            worksheet.set_column(1, 4, 15)
            worksheet.set_column(5, 255, 30)
            worksheet.set_row(2, 30)
            # worksheet.set_default_row(64)
            # my_format = workbook.add_format()
            # my_format.set_align('vtop')
            # my_format.set_align('hleft')
            # worksheet.set_column('A:XFD', None, my_format)
            if self.nd.int_all_exps == True:
                exps = range(len(self.nd.nmrdat[ds[k]]))
            else:
                exps = [self.nd.e]

            abc_string = []
            for s in itertools.islice(self.iter_all_strings(), 5 + len(exps)):
                abc_string.append(s)

            worksheet.write('B1', 'peak_label')
            worksheet.write('C1', 'peak_max_ppm')
            worksheet.write('D1', 'start_peakPPM')
            worksheet.write('E1', 'stopPeakPPM')
            worksheet.write('A2', 'Sample #')
            worksheet.write('A3', 'Sample ID')
            for m in range(len(exps)):
                worksheet.write(abc_string[m + 5] + '1', 'peak_int(exp ' + str(exps[m] + 1) + ')')
                worksheet.write(abc_string[m + 5] + '2', str(exps[m] + 1))
                worksheet.write(abc_string[m + 5] + '3', self.nd.nmrdat[k][exps[m]].title)

            for l in range(len(self.nd.nmrdat[self.nd.s][self.nd.e].start_peak)):
                worksheet.write('B' + str(l + 4), self.nd.nmrdat[self.nd.s][self.nd.e].peak_label[l])
                worksheet.write('C' + str(l + 4), self.nd.nmrdat[self.nd.s][self.nd.e].peak_max_ppm[l])
                worksheet.write('D' + str(l + 4), self.nd.nmrdat[self.nd.s][self.nd.e].start_peak[l])
                worksheet.write('E' + str(l + 4), self.nd.nmrdat[self.nd.s][self.nd.e].end_peak[l])
                for m in range(len(exps)):
                    worksheet.write(abc_string[m + 5] + str(l + 4), self.nd.nmrdat[k][exps[m]].peak_int[l])

        workbook.close()
        # end export_peak

    def fill_peak_numbers(self):
        self.nd.peak_fill = True
        s = self.nd.s
        e = self.nd.e
        n_peaks = len(self.nd.nmrdat[s][e].start_peak)
        start_peak = self.nd.nmrdat[s][e].start_peak
        end_peak = self.nd.nmrdat[s][e].end_peak
        peak_label = self.nd.nmrdat[s][e].peak_label
        self.w.peakWidget.setRowCount(n_peaks)
        for k in range(n_peaks):
            # peakNumber = QTableWidgetItem(str(k))
            # peakNumber.setTextAlignment(QtCore.Qt.AlignHCenter)
            # self.w.peakWidget.setItem(k, 0, peakNumber)
            # print('start: {}. end: {}, label: {}'.format(start_peak[k], end_peak[k], peak_label[k]))
            start_peak_tw = QTableWidgetItem(str(start_peak[k]))
            end_peak_tw = QTableWidgetItem(str(end_peak[k]))
            peak_label_tw = QTableWidgetItem(str(peak_label[k]))
            start_peak_tw.setTextAlignment(QtCore.Qt.AlignHCenter)
            end_peak_tw.setTextAlignment(QtCore.Qt.AlignHCenter)
            peak_label_tw.setTextAlignment(QtCore.Qt.AlignHCenter)
            self.w.peakWidget.setItem(k, 0, start_peak_tw)
            self.w.peakWidget.setItem(k, 1, end_peak_tw)
            self.w.peakWidget.setItem(k, 2, peak_label_tw)

        self.nd.peak_fill = False
        # end fill_pre_processing_numbers

    def fill_pre_processing_numbers(self):
        self.nd.pp.pre_proc_fill = True
        n_spc = len(self.nd.pp.class_select)
        self.w.selectClassTW.setRowCount(n_spc)
        for k in range(n_spc):
            spc_number = QTableWidgetItem(str(k))
            spc_number.setTextAlignment(QtCore.Qt.AlignHCenter)
            self.w.selectClassTW.setItem(k, 0, spc_number)
            # self.w.selectClassTW.setItemSelected(spc_number, False)
            class_number = QTableWidgetItem(self.nd.pp.class_select[k])
            class_number.setTextAlignment(QtCore.Qt.AlignHCenter)
            self.w.selectClassTW.setItem(k, 1, class_number)

        self.w.selectClassTW.selectAll()
        sel_it = self.w.selectClassTW.selectedItems()
        for k in np.arange(len(sel_it) - 1, -1, -1):
            if (self.w.selectClassTW.selectedItems()[k].column() == 1):
                self.w.selectClassTW.selectedItems()[k].setSelected(False)

        sel_it = self.w.selectClassTW.selectedItems()
        for k in np.arange(len(sel_it) - 1, -1, -1):
            if (np.isin(self.w.selectClassTW.selectedItems()[k].row(), self.nd.pp.plot_select)):
                self.w.selectClassTW.selectedItems()[k].setSelected(True)
            else:
                self.w.selectClassTW.selectedItems()[k].setSelected(False)

        for k in range(len(self.nd.pp.exclude_start)):
            excl_number1 = QTableWidgetItem(str(2 * k))
            excl_number1.setTextAlignment(QtCore.Qt.AlignHCenter)
            excl_number2 = QTableWidgetItem(str(2 * k + 1))
            excl_number2.setTextAlignment(QtCore.Qt.AlignHCenter)
            self.w.excludeRegionTW.setItem(k, 0, excl_number1)
            self.w.excludeRegionTW.setItem(k, 1, excl_number2)
            self.w.excludeRegionTW.item(k, 0).setText(str(self.nd.pp.exclude_start[k]))
            self.w.excludeRegionTW.item(k, 1).setText(str(self.nd.pp.exclude_end[k]))

        for k in range(len(self.nd.pp.seg_start)):
            seg_number1 = QTableWidgetItem(str(2 * k))
            seg_number1.setTextAlignment(QtCore.Qt.AlignHCenter)
            seg_number2 = QTableWidgetItem(str(2 * k + 1))
            seg_number2.setTextAlignment(QtCore.Qt.AlignHCenter)
            self.w.segAlignTW.setItem(k, 0, seg_number1)
            self.w.segAlignTW.setItem(k, 1, seg_number2)
            self.w.segAlignTW.item(k, 0).setText(str(self.nd.pp.seg_start[k]))
            self.w.segAlignTW.item(k, 1).setText(str(self.nd.pp.seg_end[k]))

        for k in range(len(self.nd.pp.compress_start)):
            comp_number1 = QTableWidgetItem(str(2 * k))
            comp_number1.setTextAlignment(QtCore.Qt.AlignHCenter)
            comp_number2 = QTableWidgetItem(str(2 * k + 1))
            comp_number2.setTextAlignment(QtCore.Qt.AlignHCenter)
            self.w.compressBucketsTW.setItem(k, 0, comp_number1)
            self.w.compressBucketsTW.setItem(k, 1, comp_number2)
            self.w.compressBucketsTW.item(k, 0).setText(str(self.nd.pp.compress_start[k]))
            self.w.compressBucketsTW.item(k, 1).setText(str(self.nd.pp.compress_end[k]))

        self.w.noiseThresholdLE.setText(str(self.nd.pp.noise_threshold))
        self.w.noiseRegionStartLE.setText(str(self.nd.pp.noise_start))
        self.w.noiseRegionEndLE.setText(str(self.nd.pp.noise_end))
        self.w.thLineWidthLE.setText(str(self.nd.pp.th_line_width))
        self.w.bucketPpmLE.setText(str(self.nd.pp.bucket_ppm))
        self.set_bucket_ppm_pre_proc()
        self.w.excludeRegion.setChecked(self.nd.pp.flag_exclude_region)
        self.w.segmentalAlignment.setChecked(self.nd.pp.flag_segmental_alignment)
        self.w.noiseFiltering.setChecked(self.nd.pp.flag_noise_filtering)
        self.w.bucketSpectra.setChecked(self.nd.pp.flag_bucket_spectra)
        self.w.compressBuckets.setChecked(self.nd.pp.flag_compress_buckets)
        self.w.scaleSpectra.setChecked(self.nd.pp.flag_scale_spectra)
        if self.nd.pp.scale_pqn is True:
            self.w.pqnButton.setChecked(True)
        else:
            self.w.tsaButton.setChecked(True)

        self.w.varianceStabilisation.setChecked(self.nd.pp.flag_variance_stabilisation)
        self.w.exportDataSet.setChecked(self.nd.pp.flag_export_data_set)
        self.w.exportDelimiterTab.setChecked(self.nd.pp.export_delimiter_tab)
        self.w.exportDelimiterCharacter.setChecked(not self.nd.pp.export_delimiter_tab)
        self.w.exportCharacter.setText(self.nd.pp.export_character)
        self.w.exportMethod.setCurrentIndex(self.nd.pp.export_method)
        self.w.samplesInComboBox.setCurrentIndex(self.nd.pp.export_samples_in_rows_cols)
        self.w.segAlignRefSpc.setMaximum(len(self.nd.nmrdat[self.nd.s]))
        self.w.scaleSpectraRefSpc.setMaximum(len(self.nd.nmrdat[self.nd.s]))
        self.w.segAlignRefSpc.setValue(self.nd.pp.seg_align_ref_spc)
        self.w.scaleSpectraRefSpc.setValue(self.nd.pp.scale_spectra_ref_spc)
        self.w.preserveOverallScale.setChecked(self.nd.pp.preserve_overall_scale)
        self.w.preserveOverallScale.setDisabled(self.nd.pp.scale_pqn)
        self.w.autoScaling.setChecked(self.nd.pp.auto_scaling)
        self.w.paretoScaling.setChecked(self.nd.pp.pareto_scaling)
        self.w.gLogTransform.setChecked(self.nd.pp.g_log_transform)
        self.w.lambdaText.setEnabled(self.nd.pp.g_log_transform)
        self.w.y0Text.setEnabled(self.nd.pp.g_log_transform)
        self.w.lambdaLE.setEnabled(self.nd.pp.g_log_transform)
        self.w.y0LE.setEnabled(self.nd.pp.g_log_transform)
        self.w.lambdaLE.setText(str(self.nd.pp.var_lambda))
        self.w.y0LE.setText(str(self.nd.pp.var_y0))
        self.nd.pp.pre_proc_fill = False
        # end fill_pre_processing_numbers

    def ft(self):
        self.nd.ft()
        if (self.w.baselineCorrection.currentIndex() > 0):
            self.baseline1d()

        self.w.nmrSpectrum.setCurrentIndex(0)
        self.change_data_set_exp()
        self.plot_spc()
        # end ft

    def ft_all(self):
        self.nd.ft_all()
        if (self.w.baselineCorrection.currentIndex() > 0):
            self.baseline1d_all()

        self.w.nmrSpectrum.setCurrentIndex(0)
        self.change_data_set_exp()
        self.plot_spc()
        # end ft_all

    def get_disp_pars1(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.pos_col = d.colours.get(self.w.posCol.currentIndex())
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        # end get_disp_pars1

    def get_disp_pars2(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.neg_col = d.colours.get(self.w.negCol.currentIndex())
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        # end get_disp_pars2

    def get_disp_pars3(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        pos_r = float(self.w.posColR.text())
        pos_g = float(self.w.posColG.text())
        pos_b = float(self.w.posColB.text())
        neg_r = float(self.w.negColR.text())
        neg_g = float(self.w.negColG.text())
        neg_b = float(self.w.negColB.text())
        d.pos_col_rgb = (pos_r, pos_g, pos_b)
        d.neg_col_rgb = (neg_r, neg_g, neg_b)
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        # end get_disp_pars3

    def get_disp_pars4(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.n_levels = round(float(self.w.nLevels.text()))
        # end get_disp_pars4

        self.nd.nmrdat[self.nd.s][self.nd.e].display = d

    def get_disp_pars5(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.min_level = float(self.w.minLevel.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        # end get_disp_pars5

    def get_disp_pars6(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.max_level = float(self.w.maxLevel.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        # end get_disp_pars6

    def get_disp_pars7(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.axis_type1 = d.axes.get(self.w.axisType1.currentIndex())
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        self.nd.nmrdat[self.nd.s][self.nd.e].calc_ppm()
        # end get_disp_pars7

    def get_disp_pars8(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.axis_type2 = d.axes.get(self.w.axisType2.currentIndex())
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        self.nd.nmrdat[self.nd.s][self.nd.e].calc_ppm()
        # end get_disp_pars8

    def get_disp_pars9(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.display_spc = d.false_true.get(self.w.displaySpc.currentIndex())
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        # end get_disp_pars9

    def get_disp_pars10(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.spc_offset = float(self.w.spcOffset.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        # end get_disp_pars10

    def get_disp_pars11(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.spc_scale = float(self.w.spcScale.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        # end get_disp_pars11

    def get_disp_pars12(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.x_label = self.w.xLabel.text()
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        # end get_disp_pars12

    def get_disp_pars13(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.y_label = self.w.yLabel.text()
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        # end get_disp_pars13

    def get_disp_pars14(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.spc_label = self.w.spcLabel.text()
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        # end get_disp_pars14

    def get_disp_pars15(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.ph_ref_col = d.colours2.get(self.w.phRefColour.currentIndex())
        for k in range(len(self.nd.nmrdat)):
            for l in range(len(self.nd.nmrdat[k])):
                self.nd.nmrdat[k][l].display.ph_ref_col = d.ph_ref_col

        # end get_disp_pars15

    def get_proc_pars1(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.window_type[0] = self.w.windowFunction.currentIndex()
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars1

    def get_proc_pars2(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.window_type[1] = self.w.windowFunction_2.currentIndex()
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars2

    def get_proc_pars3(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.ph_corr[0] = self.w.phaseCorrection.currentIndex()
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars3

    def get_proc_pars4(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.ph_corr[1] = self.w.phaseCorrection_2.currentIndex()
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars4

    def get_proc_pars5(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.water_suppression = self.w.waterSuppression.currentIndex()
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars5

    def get_proc_pars6(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.conv_window_type[0] = self.w.winType.currentIndex()
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars6

    def get_proc_pars7(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.gibbs[0] = p.gibbs_p.get(self.w.gibbs.currentIndex())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars7

    def get_proc_pars8(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.gibbs[1] = p.gibbs_p.get(self.w.gibbs_2.currentIndex())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars8

    def get_proc_pars9(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.n_points[0] = int(self.w.zeroFilling.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars9

    def get_proc_pars10(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.n_points[1] = int(self.w.zeroFilling_2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars10

    def get_proc_pars11(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.lb[0] = float(self.w.lb.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars11

    def get_proc_pars12(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.gb[0] = float(self.w.gb.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars12

    def get_proc_pars13(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.ssb[0] = float(self.w.ssb.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars13

    def get_proc_pars14(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.lb[1] = float(self.w.lb_2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars14

    def get_proc_pars15(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.gb[1] = float(self.w.gb_2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars15

    def get_proc_pars16(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.ssb[1] = float(self.w.ssb_2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars16

    def get_proc_pars17(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.ph0[0] = float(self.w.ph0.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars17

    def get_proc_pars18(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.ph1[0] = float(self.w.ph1.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars18

    def get_proc_pars19(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.ph0[1] = float(self.w.ph0_2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars19

    def get_proc_pars20(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.ph1[1] = float(self.w.ph1_2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars20

    def get_proc_pars21(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.poly_order = int(self.w.polyOrder.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars21

    def get_proc_pars22(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.conv_extrapolation_size[0] = int(self.w.extrapolationSize.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars22

    def get_proc_pars23(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.conv_window_size[0] = int(self.w.windowSize.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars23

    def get_proc_pars24(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.fid_offset_correction = int(self.w.fidOffsetCorrection.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars24

    def get_proc_pars25(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.strip_start = int(self.w.stripTransformStart.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars25

    def get_proc_pars26(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.strip_end = int(self.w.stripTransformEnd.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end get_proc_pars25

    def get_bottom_top(self, line):
        margin = 0.1
        xd = line.get_xdata()
        yd = line.get_ydata()
        lo, hi = self.w.MplWidget.canvas.axes.get_xlim()
        y_displayed = yd[((xd > min(lo, hi)) & (xd < max(lo, hi)))]
        h = np.max(y_displayed) - np.min(y_displayed)
        bot = np.min(y_displayed) - margin * h
        top = np.max(y_displayed) + margin * h
        return bot, top

    def get_r_spc_p0(self):
        r = self.nd.nmrdat[self.nd.s][self.nd.e].apc.r_spc
        r[0] = float(self.w.rSpc_p0.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.r_spc = r
        # end get_r_spc_p0

    def get_r_spc_p1(self):
        r = self.nd.nmrdat[self.nd.s][self.nd.e].apc.r_spc
        r[1] = float(self.w.rSpc_p1.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.r_spc = r
        # end get_r_spc_p1

    def get_r_spc_p2(self):
        r = self.nd.nmrdat[self.nd.s][self.nd.e].apc.r_spc
        r[2] = float(self.w.rSpc_p2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.r_spc = r
        # end get_r_spc_p2

    def get_r_spc_p3(self):
        r = self.nd.nmrdat[self.nd.s][self.nd.e].apc.r_spc
        r[3] = float(self.w.rSpc_p3.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.r_spc = r
        # end get_r_spc_p3

    def get_r_spc_p4(self):
        r = self.nd.nmrdat[self.nd.s][self.nd.e].apc.r_spc
        r[4] = float(self.w.rSpc_p4.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.r_spc = r
        # end get_r_spc_p4

    def get_r_spc_p5(self):
        r = self.nd.nmrdat[self.nd.s][self.nd.e].apc.r_spc
        r[5] = float(self.w.rSpc_p5.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.r_spc = r
        # end get_r_spc_p5

    def get_r_spc_p6(self):
        r = self.nd.nmrdat[self.nd.s][self.nd.e].apc.r_spc
        r[6] = float(self.w.rSpc_p6.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.r_spc = r
        # end get_r_spc_p6

    def get_i_spc_p0(self):
        i = self.nd.nmrdat[self.nd.s][self.nd.e].apc.i_spc
        i[0] = float(self.w.iSpc_p0.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.r_spc = i
        # end get_i_spc_p0

    def get_i_spc_p1(self):
        i = self.nd.nmrdat[self.nd.s][self.nd.e].apc.i_spc
        i[1] = float(self.w.iSpc_p1.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.r_spc = i
        # end get_i_spc_p1

    def get_i_spc_p2(self):
        i = self.nd.nmrdat[self.nd.s][self.nd.e].apc.i_spc
        i[2] = float(self.w.iSpc_p2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.r_spc = i
        # end get_i_spc_p2

    def get_i_spc_p3(self):
        i = self.nd.nmrdat[self.nd.s][self.nd.e].apc.i_spc
        i[3] = float(self.w.iSpc_p3.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.r_spc = i
        # end get_i_spc_p3

    def get_i_spc_p4(self):
        i = self.nd.nmrdat[self.nd.s][self.nd.e].apc.i_spc
        i[4] = float(self.w.iSpc_p4.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.r_spc = i
        # end get_i_spc_p4

    def get_i_spc_p5(self):
        i = self.nd.nmrdat[self.nd.s][self.nd.e].apc.i_spc
        i[5] = float(self.w.iSpc_p5.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.r_spc = i
        # end get_i_spc_p5

    def get_i_spc_p6(self):
        i = self.nd.nmrdat[self.nd.s][self.nd.e].apc.i_spc
        i[6] = float(self.w.iSpc_p6.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.r_spc = i
        # end get_i_spc_p6

    def on_g_input_click(self, event):
        self.cur_clicks += 1
        if self.cur_clicks < self.n_clicks:
            self.xdata.append(event.xdata)
            self.ydata.append(event.ydata)
        else:
            self.xdata.append(event.xdata)
            self.ydata.append(event.ydata)
            self.n_clicks = 1
            self.cur_clicks = 0
            cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.on_g_input_click)
            cid = self.w.MplWidget.canvas.mpl_disconnect(cid)
            cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.on_g_input_click)
            cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
            code_out = io.StringIO()
            code_err = io.StringIO()
            sys.stdout = code_out
            sys.stderr = code_err
            print("x-values: {} / xDiff [ppm]: {} / xDiff [Hz]: {}".format(self.xdata, np.abs(np.diff(self.xdata)),
                                                                           np.abs(np.diff(self.xdata)) *
                                                                           self.nd.nmrdat[self.nd.s][
                                                                               self.nd.e].acq.sfo1))
            if self.nd.nmrdat[self.nd.s][self.nd.e].dim == 1:
                print("y-values: {} / yDiff: {}".format(self.ydata, -np.diff(self.ydata)))
            else:
                print("y-values: {} / yDiff: {} / yDiff [Hz]: {}".format(self.ydata, np.abs(np.diff(self.ydata)),
                                                                         np.abs(np.diff(self.ydata)) *
                                                                         self.nd.nmrdat[self.nd.s][self.nd.e].acq.sfo2))

            self.w.console.append(code_out.getvalue())
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__
            code_out.close()
            code_err.close()
            self.xdata = []
            self.ydata = []
            self.show_console()

    def on_g_input_click_add_peak(self, event):
        self.cur_clicks += 1
        if self.cur_clicks < self.n_clicks:
            self.xdata.append(event.xdata)
            self.ydata.append(event.ydata)
        else:
            self.xdata.append(event.xdata)
            self.ydata.append(event.ydata)
            self.n_clicks = 1
            self.cur_clicks = 0
            cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.on_g_input_click_add_peak)
            cid = self.w.MplWidget.canvas.mpl_disconnect(cid)
            cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.on_g_input_click_add_peak)
            cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
            xy = []
            for k in range(len(self.xdata)):
                xy.append((self.xdata[k], self.ydata[k]))

            self.xdata = []
            self.ydata = []
            t = np.round(1e4 * np.array([xy[0][0], xy[1][0]])) / 1e4
            self.nd.add_peak(t, '')
            # self.nd.pp.exclude_start = np.append(self.nd.pp.exclude_start, min(t))
            # self.nd.pp.exclude_end = np.append(self.nd.pp.exclude_end, max(t))
            self.fill_peak_numbers()
            # self.w.excludeRegionTW.setFocus()
            # self.set_plot_pre_proc()
            # self.w.excludeRegionTW.setFocus()
            self.plot_spc()
            # self.set_exclude_pre_proc()

    def on_g_input_click_ref_1d(self, event):
        self.cur_clicks += 1
        if self.cur_clicks < self.n_clicks:
            self.xdata.append(event.xdata)
            self.ydata.append(event.ydata)
        else:
            self.xdata.append(event.xdata)
            self.ydata.append(event.ydata)
            self.n_clicks = 1
            self.cur_clicks = 0
            cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.on_g_input_click_ref_1d)
            cid = self.w.MplWidget.canvas.mpl_disconnect(cid)
            cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.on_g_input_click_ref_1d)
            cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
            xy = []
            for k in range(len(self.xdata)):
                xy.append((self.xdata[k], self.ydata[k]))

            self.xdata = []
            self.ydata = []
            self.nd.nmrdat[self.nd.s][self.nd.e].ref_point[0] = self.nd.nmrdat[self.nd.s][self.nd.e].ppm2points(
                xy[0][0],
                0)
            self.nd.nmrdat[self.nd.s][self.nd.e].ref_shift[0] = self.temp_ref_shift
            self.nd.nmrdat[self.nd.s][self.nd.e].calc_ppm()
            self.reset_plot()

    def on_g_input_click_compress(self, event):
        self.cur_clicks += 1
        if self.cur_clicks < self.n_clicks:
            self.xdata.append(event.xdata)
            self.ydata.append(event.ydata)
        else:
            self.xdata.append(event.xdata)
            self.ydata.append(event.ydata)
            self.n_clicks = 1
            self.cur_clicks = 0
            cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.on_g_input_click_compress)
            cid = self.w.MplWidget.canvas.mpl_disconnect(cid)
            cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.on_g_input_click_compress)
            cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
            xy = []
            for k in range(len(self.xdata)):
                xy.append((self.xdata[k], self.ydata[k]))

            self.xdata = []
            self.ydata = []
            t = np.round(1e4 * np.array([xy[0][0], xy[1][0]])) / 1e4
            self.nd.pp.compress_start = np.append(self.nd.pp.compress_start, min(t))
            self.nd.pp.compress_end = np.append(self.nd.pp.compress_end, max(t))
            self.fill_pre_processing_numbers()
            self.w.compressBucketsTW.setFocus()
            self.set_plot_pre_proc()
            self.w.compressBucketsTW.setFocus()
            self.plot_spc_pre_proc()
            self.set_compress_pre_proc()

    def on_g_input_click_exclude(self, event):
        self.cur_clicks += 1
        if self.cur_clicks < self.n_clicks:
            self.xdata.append(event.xdata)
            self.ydata.append(event.ydata)
        else:
            self.xdata.append(event.xdata)
            self.ydata.append(event.ydata)
            self.n_clicks = 1
            self.cur_clicks = 0
            cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.on_g_input_click_exclude)
            cid = self.w.MplWidget.canvas.mpl_disconnect(cid)
            cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.on_g_input_click_exclude)
            cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
            xy = []
            for k in range(len(self.xdata)):
                xy.append((self.xdata[k], self.ydata[k]))

            self.xdata = []
            self.ydata = []
            t = np.round(1e4 * np.array([xy[0][0], xy[1][0]])) / 1e4
            self.nd.pp.exclude_start = np.append(self.nd.pp.exclude_start, min(t))
            self.nd.pp.exclude_end = np.append(self.nd.pp.exclude_end, max(t))
            self.fill_pre_processing_numbers()
            self.w.excludeRegionTW.setFocus()
            self.set_plot_pre_proc()
            self.w.excludeRegionTW.setFocus()
            self.plot_spc_pre_proc()
            self.set_exclude_pre_proc()

    def on_g_input_click_seg_align(self, event):
        self.cur_clicks += 1
        if self.cur_clicks < self.n_clicks:
            self.xdata.append(event.xdata)
            self.ydata.append(event.ydata)
        else:
            self.xdata.append(event.xdata)
            self.ydata.append(event.ydata)
            self.n_clicks = 1
            self.cur_clicks = 0
            cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.on_g_input_click_seg_align)
            cid = self.w.MplWidget.canvas.mpl_disconnect(cid)
            cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.on_g_input_click_seg_align)
            cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
            xy = []
            for k in range(len(self.xdata)):
                xy.append((self.xdata[k], self.ydata[k]))

            self.xdata = []
            self.ydata = []
            t = np.round(1e4 * np.array([xy[0][0], xy[1][0]])) / 1e4
            self.nd.pp.seg_start = np.append(self.nd.pp.seg_start, min(t))
            self.nd.pp.seg_end = np.append(self.nd.pp.seg_end, max(t))
            self.fill_pre_processing_numbers()
            self.w.segAlignTW.setFocus()
            self.set_plot_pre_proc()
            self.w.segAlignTW.setFocus()
            self.plot_spc_pre_proc()
            self.set_seg_align_pre_proc()

    def on_g_input_click_ref_2d(self, event):
        self.cur_clicks += 1
        if self.cur_clicks < self.n_clicks:
            self.xdata.append(event.xdata)
            self.ydata.append(event.ydata)
        else:
            self.xdata.append(event.xdata)
            self.ydata.append(event.ydata)
            self.n_clicks = 1
            self.cur_clicks = 0
            cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.on_g_input_click_ref_2d)
            cid = self.w.MplWidget.canvas.mpl_disconnect(cid)
            cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.on_g_input_click_ref_2d)
            cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
            xy = []
            for k in range(len(self.xdata)):
                xy.append((self.xdata[k], self.ydata[k]))

            self.xdata = []
            self.ydata = []
            self.nd.nmrdat[self.nd.s][self.nd.e].ref_point[0] = self.nd.nmrdat[self.nd.s][self.nd.e].ppm2points(
                xy[0][0],
                0)
            self.nd.nmrdat[self.nd.s][self.nd.e].ref_shift[0] = self.temp_ref_shift[0]
            self.nd.nmrdat[self.nd.s][self.nd.e].ref_point[1] = self.nd.nmrdat[self.nd.s][self.nd.e].ppm2points(
                xy[0][1],
                1)
            self.nd.nmrdat[self.nd.s][self.nd.e].ref_shift[1] = self.temp_ref_shift[1]
            self.nd.nmrdat[self.nd.s][self.nd.e].proc.ref_point[0] = self.nd.nmrdat[self.nd.s][self.nd.e].ref_point[0] * \
                                                                     self.nd.nmrdat[self.nd.s][self.nd.e].proc.n_points[
                                                                         0] / (len(
                self.nd.nmrdat[self.nd.s][self.nd.e].fid[0]) * self.nd.nmrdat[self.nd.s][self.nd.e].proc.mult_factor[0])
            self.nd.nmrdat[self.nd.s][self.nd.e].proc.ref_point[1] = self.nd.nmrdat[self.nd.s][self.nd.e].ref_point[1] * \
                                                                     self.nd.nmrdat[self.nd.s][self.nd.e].proc.n_points[
                                                                         1] / (len(
                self.nd.nmrdat[self.nd.s][self.nd.e].fid) * self.nd.nmrdat[self.nd.s][self.nd.e].proc.mult_factor[1])
            self.nd.nmrdat[self.nd.s][self.nd.e].calc_ppm()
            self.reset_plot()

    def on_g_input_2d_click(self, event):
        self.cur_clicks += 1
        if self.cur_clicks < self.n_clicks:
            self.xdata.append(event.xdata)
            self.ydata.append(event.ydata)
        else:
            self.xdata.append(event.xdata)
            self.ydata.append(event.ydata)
            self.n_clicks = 1
            self.cur_clicks = 0
            cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.on_g_input_2d_click)
            cid = self.w.MplWidget.canvas.mpl_disconnect(cid)
            cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.on_g_input_2d_click)
            cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
            self.xy = []
            self.xy = np.resize(self.xy, (1, 2))
            self.xy[0][0] = self.xdata[0]
            self.xy[0][1] = self.ydata[0]
            self.xdata = []
            self.ydata = []
            xy = self.xy
            self.show_ph_corr2d()
            xyPts = []
            xy2 = []
            xyPts.append(self.nd.nmrdat[self.nd.s][self.nd.e].ppm2points(xy[0][0], 0))
            xyPts.append(self.nd.nmrdat[self.nd.s][self.nd.e].ppm2points(xy[0][1], 1))
            self.ph_corr.spc_row_pts.append(xyPts[1])
            self.ph_corr.spc_col_pts.append(xyPts[0])
            xy2.append(self.nd.nmrdat[self.nd.s][self.nd.e].points2ppm(xyPts[0], 0))
            xy2.append(self.nd.nmrdat[self.nd.s][self.nd.e].points2ppm(xyPts[1], 1))
            self.ph_corr.spc_row.append(xy2[1])
            self.ph_corr.spc_col.append(xy2[0])
            self.plot_2d_col_row()

    def ginput(self, n_clicks=1):
        self.w.MplWidget.canvas.setFocus()
        self.show_nmr_spectrum()
        self.n_clicks = n_clicks
        cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.on_g_input_click)
        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.on_g_input_click)
        cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
        # end ginput

    def ginput_add_peak(self, n_clicks=2):
        self.w.MplWidget.canvas.setFocus()
        self.show_nmr_spectrum()
        self.n_clicks = n_clicks
        cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.on_g_input_click_add_peak)
        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.on_g_input_click_add_peak)
        cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
        # end ginput

    def ginput_ref_1d(self, n_clicks=1):
        self.w.MplWidget.canvas.setFocus()
        self.show_nmr_spectrum()
        self.n_clicks = n_clicks
        # print("ginput_ref_1d")
        cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.on_g_input_click_ref_1d)
        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.on_g_input_click_ref_1d)
        cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
        # end ginput

    def ginput_ref_2d(self, n_clicks=1):
        self.w.MplWidget.canvas.setFocus()
        self.show_nmr_spectrum()
        self.n_clicks = n_clicks
        cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.on_g_input_click_ref_2d)
        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.on_g_input_click_ref_2d)
        cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
        # end ginput

    def ginput_compress(self, n_clicks=1):
        self.w.MplWidget.canvas.setFocus()
        self.show_nmr_spectrum()
        self.n_clicks = n_clicks
        cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.on_g_input_click_compress)
        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.on_g_input_click_compress)
        cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
        # end ginput

    def ginput_exclude(self, n_clicks=2):
        self.w.MplWidget.canvas.setFocus()
        self.show_nmr_spectrum()
        self.n_clicks = n_clicks
        cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.on_g_input_click_exclude)
        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.on_g_input_click_exclude)
        cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
        # end ginput

    def ginput_seg_align(self, n_clicks=1):
        self.w.MplWidget.canvas.setFocus()
        self.show_nmr_spectrum()
        self.n_clicks = n_clicks
        cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.on_g_input_click_seg_align)
        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.on_g_input_click_seg_align)
        cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
        # end ginput

    def ginput2d(self):
        self.w.MplWidget.canvas.setFocus()
        self.show_nmr_spectrum()
        self.n_clicks = 1
        cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.on_g_input_2d_click)
        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.on_g_input_2d_click)
        cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
        # end ginput2d

    def h(self):
        print("Command history: ")
        print(">>><<<")
        for k in range(len(self.nd.cmd_buffer)):
            print(self.nd.cmd_buffer[k])

        return (">>><<<")
        # end h

    def hide_pre_processing(self):
        self.w.preProcessingTab.setHidden(True)
        self.w.preProcPeak.setTabEnabled(0, False)
        self.w.preProcessingGroupBox.setHidden(True)
        self.w.preProcessingSelect.setHidden(True)
        self.w.preProcessingWidget.setHidden(True)
        self.w.runPreProcessingButton.setHidden(True)
        self.w.resetPreProcessingButton.setHidden(True)
        self.w.writeScriptButton.setHidden(True)
        self.plot_spc(True)
        # end hide_pre_processing

    def hide_peak_picking(self):
        self.w.peakPickingTab.setHidden(True)
        self.w.preProcPeak.setTabEnabled(1, False)
        self.w.peakWidget.setHidden(True)
        self.w.peakAddButton.setHidden(True)
        self.w.peakClearButton.setHidden(True)
        self.w.peakExportButton.setHidden(True)
        self.w.intAllExps.setHidden(True)
        self.w.intAllDS.setHidden(True)
        self.w.exportFormatCB.setHidden(True)

    def show_peak_picking(self):
        self.w.preProcPeak.setHidden(False)
        self.w.peakPickingTab.setHidden(False)
        self.w.preProcPeak.setTabEnabled(0, False)
        self.w.preProcPeak.setTabEnabled(1, True)
        self.w.peakWidget.setHidden(False)
        self.w.peakAddButton.setHidden(False)
        self.w.peakClearButton.setHidden(False)
        self.w.peakExportButton.setHidden(False)
        self.w.intAllExps.setHidden(False)
        self.w.intAllDS.setHidden(False)
        self.w.exportFormatCB.setHidden(False)

    def hilbert(self, mat):
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
        return mat
        # end hilbert

    def horz_ph_corr_2d(self):
        s = self.nd.s
        e = self.nd.e
        self.ph_corr.n_dims = 2
        self.ph_corr.dim = 0
        n_lines = len(self.ph_corr.spc_row_pts)
        if n_lines > 0:
            npts0 = len(self.nd.nmrdat[s][e].spc)
            npts = len(self.nd.nmrdat[s][e].spc[0])
            self.ph_corr.spc = np.zeros((n_lines, npts), dtype='complex')
            spc1 = np.copy(self.nd.nmrdat[s][e].spc)
            for k in range(n_lines):
                spc = np.array([spc1[npts0 - self.ph_corr.spc_row_pts[k]]])
                spc = self.hilbert(spc)
                self.ph_corr.spc[k] = spc[0]

            self.ph_corr.ppm = self.nd.nmrdat[s][e].ppm1
            if self.ph_corr.pivot_points2d[0] < 0:
                self.ph_corr.pivot_points2d[0] = int(len(self.ph_corr.ppm) / 2)
                self.ph_corr.pivot2d[0] = self.nd.nmrdat[s][e].points2ppm(self.ph_corr.pivot_points2d[0], 0)

        self.show_ph_corr2d_1d(self.ph_corr.dim)
        self.ph_corr.spc_max = np.max(np.max(np.abs(self.ph_corr.spc)))
        try:
            zwo = True
            self.w.MplWidget.canvas.figure.canvas.toolbar.zoom()
        except:
            pass

        self.set_zoom_off()
        self.ph_corr.max_ph0 = 90.0
        self.ph_corr.max_ph1 = 90.0
        cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.on_ph_corr_click_2d)
        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.on_ph_corr_release_2d)
        # self.w.actionApplyPhCorr.triggered.connect(self.apply_2d_ph_corr)
        # self.w.actionCancelPhCorr.triggered.connect(self.cancel2dPhCorr)
        self.w.pickRowColPhCorr2d.setVisible(False)
        self.w.emptyRowColPhCorr2d.setVisible(False)
        self.w.removeRowColPhCorr2d.setVisible(False)
        self.w.horzPhCorr2d.setVisible(False)
        self.w.vertPhCorr2d.setVisible(False)
        self.w.zoomPhCorr2d.setVisible(True)
        self.w.applyPhCorr2d.setVisible(True)
        self.w.cancelPhCorr2d.setVisible(True)
        self.w.exitPhCorr2d.setVisible(False)
        self.w.exitZoomPhCorr2d.setVisible(False)
        self.ph_corr_plot_spc_2d(False)
        self.show_acquisition_parameters()
        self.show_nmr_spectrum()
        # end horz_ph_corr_2d

    def iter_all_strings(self):
        for size in itertools.count(1):
            for s in itertools.product(ascii_uppercase, repeat=size):
                yield "".join(s)

        # end iter_all_strings+

    def load_button(self):
        selectedDirectory = QFileDialog.getExistingDirectory()
        if (len(selectedDirectory) > 0):
            self.clear()
            self.zero_script()
        else:
            return

        self.load_file(selectedDirectory)
        # end save_button

    def load_config(self):
        self.cf.read_config()
        self.w.phRefColour.setCurrentIndex(self.nd.nmrdat[0][0].display.colours2.get(self.cf.phase_reference_colour))
        self.w.autoPlot.setChecked(self.cf.auto_plot)
        self.w.keepZoom.setChecked(self.cf.keep_zoom)
        self.w.fontSize.setValue(self.cf.font_size)
        self.std_pos_col1 = (self.cf.pos_col10, self.cf.pos_col11, self.cf.pos_col12)
        self.std_neg_col1 = (self.cf.neg_col10, self.cf.neg_col11, self.cf.neg_col12)
        self.std_pos_col2 = (self.cf.pos_col20, self.cf.pos_col21, self.cf.pos_col22)
        self.std_neg_col2 = (self.cf.neg_col20, self.cf.neg_col21, self.cf.neg_col22)
        self.set_standard_colours()
        # end load_config

    def load_example_script(self):
        idx = self.w.exampleScripts.view().selectedIndexes()[0].row()
        self.w.exampleScripts.setCurrentIndex(idx)
        if (idx == 0):
            f_name = os.path.join(os.path.dirname(__file__), "exampleScripts", "example1DScript.py")

        if (idx == 1):
            f_name = os.path.join(os.path.dirname(__file__), "exampleScripts", "exampleAutoPhaseScript.py")

        if (idx == 2):
            f_name = os.path.join(os.path.dirname(__file__), "exampleScripts", "example2DJresScript.py")

        if (idx == 3):
            f_name = os.path.join(os.path.dirname(__file__), "exampleScripts", "examplePreprocessingScript.py")

        if (idx == 4):
            f_name = os.path.join(os.path.dirname(__file__), "exampleScripts", "example2DNMRPipeScript.py")

        f = open(f_name, 'r')
        script_text = f.read()
        self.w.script.setText(script_text)
        # end load_example_script

    def load_dark_mode(self):
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            self.nd.nmrdat[self.nd.s][k].display.pos_col_rgb = self.std_pos_col2
            self.nd.nmrdat[self.nd.s][k].display.neg_col_rgb = self.std_neg_col2

        idx = self.w.helpComboBox.currentIndex()
        url = []
        f_name = os.path.join(os.path.dirname(__file__), "nmr", "web", "introductionDark", "index.html")
        url.append("file:///" + f_name.replace('\\', '/'))
        url.append("https://www.hmdb.ca")
        url.append("https://www.smpdb.ca")
        url.append("https://bmrb.io/metabolomics/")
        url.append("https://www.genome.jp/kegg/pathway.html#metabolism")
        url.append("https://nmrshiftdb.nmr.uni-koeln.de")
        url.append("https://sdbs.db.aist.go.jp/sdbs/cgi-bin/cre_index.cgi")
        url.append("http://dmar.riken.jp/spincouple/")
        self.w.helpView.setUrl(url[idx])
        bg = (42 / 255, 42 / 255, 42 / 255)
        fg = (255 / 255, 255 / 255, 255 / 255)
        self.w.MplWidget.canvas.figure.set_facecolor(bg)
        self.w.MplWidget.canvas.axes.set_facecolor(bg)
        self.w.MplWidget.canvas.axes.xaxis.label.set_color(fg)
        self.w.MplWidget.canvas.axes.yaxis.label.set_color(fg)
        self.w.MplWidget.canvas.axes.tick_params(axis='x', colors=fg)
        self.w.MplWidget.canvas.axes.tick_params(axis='y', colors=fg)
        self.w.MplWidget.canvas.axes.spines['bottom'].set_color(fg)
        self.w.MplWidget.canvas.axes.spines['top'].set_color(fg)
        self.w.MplWidget.canvas.axes.spines['left'].set_color(fg)
        self.w.MplWidget.canvas.axes.spines['right'].set_color(fg)
        self.w.MplWidget2.canvas.figure.set_facecolor(bg)
        self.w.MplWidget2.canvas.axes.set_facecolor(bg)
        self.w.MplWidget2.canvas.axes.xaxis.label.set_color(fg)
        self.w.MplWidget2.canvas.axes.yaxis.label.set_color(fg)
        self.w.MplWidget2.canvas.axes.tick_params(axis='x', colors=fg)
        self.w.MplWidget2.canvas.axes.tick_params(axis='y', colors=fg)
        self.w.MplWidget2.canvas.axes.spines['bottom'].set_color(fg)
        self.w.MplWidget2.canvas.axes.spines['top'].set_color(fg)
        self.w.MplWidget2.canvas.axes.spines['left'].set_color(fg)
        self.w.MplWidget2.canvas.axes.spines['right'].set_color(fg)
        self.w.MplWidget3.canvas.figure.set_facecolor(bg)
        self.w.MplWidget3.canvas.axes.set_facecolor(bg)
        self.w.MplWidget3.canvas.axes.xaxis.label.set_color(fg)
        self.w.MplWidget3.canvas.axes.yaxis.label.set_color(fg)
        self.w.MplWidget3.canvas.axes.tick_params(axis='x', colors=fg)
        self.w.MplWidget3.canvas.axes.tick_params(axis='y', colors=fg)
        self.w.MplWidget3.canvas.axes.spines['bottom'].set_color(fg)
        self.w.MplWidget3.canvas.axes.spines['top'].set_color(fg)
        self.w.MplWidget3.canvas.axes.spines['left'].set_color(fg)
        self.w.MplWidget3.canvas.axes.spines['right'].set_color(fg)
        # end load_dark_mode

    def set_standard_plot_colours(self):
        self.cf.read_config()
        self.std_pos_col1 = (self.cf.pos_col10, self.cf.pos_col11, self.cf.pos_col12)
        self.std_neg_col1 = (self.cf.neg_col10, self.cf.neg_col11, self.cf.neg_col12)
        self.std_pos_col2 = (self.cf.pos_col20, self.cf.pos_col21, self.cf.pos_col22)
        self.std_neg_col2 = (self.cf.neg_col20, self.cf.neg_col21, self.cf.neg_col22)
        self.set_standard_colours()

    def load_file(self, fileName):
        self.nd.load(fileName)
        self.w.script.insertHtml(self.nd.script)
        self.w.console.insertHtml(self.nd.console)
        self.reset_plot()
        self.update_gui()
        self.w.console.verticalScrollBar().setValue(self.w.console.verticalScrollBar().maximum())
        self.show_title_file_information()
        self.show_acquisition_parameters()
        self.show_nmr_spectrum()
        # end load_file

    def load_light_mode(self):
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            self.nd.nmrdat[self.nd.s][k].display.pos_col_rgb = self.std_pos_col1
            self.nd.nmrdat[self.nd.s][k].display.neg_col_rgb = self.std_neg_col1

        idx = self.w.helpComboBox.currentIndex()
        url = []
        f_name = os.path.join(os.path.dirname(__file__), "nmr", "web", "introduction", "index.html")
        url.append("file:///" + f_name.replace('\\', '/'))
        url.append("https://www.hmdb.ca")
        url.append("https://www.smpdb.ca")
        url.append("https://bmrb.io/metabolomics/")
        url.append("https://www.genome.jp/kegg/pathway.html#metabolism")
        url.append("https://nmrshiftdb.nmr.uni-koeln.de")
        url.append("https://sdbs.db.aist.go.jp/sdbs/cgi-bin/cre_index.cgi")
        url.append("http://dmar.riken.jp/spincouple/")
        self.w.helpView.setUrl(url[idx])
        bg = (255 / 255, 255 / 255, 255 / 255)
        fg = (0 / 255, 0 / 255, 0 / 255)
        self.w.MplWidget.canvas.figure.set_facecolor(bg)
        self.w.MplWidget.canvas.axes.set_facecolor(bg)
        self.w.MplWidget.canvas.axes.xaxis.label.set_color(fg)
        self.w.MplWidget.canvas.axes.yaxis.label.set_color(fg)
        self.w.MplWidget.canvas.axes.tick_params(axis='x', colors=fg)
        self.w.MplWidget.canvas.axes.tick_params(axis='y', colors=fg)
        self.w.MplWidget.canvas.axes.spines['bottom'].set_color(fg)
        self.w.MplWidget.canvas.axes.spines['top'].set_color(fg)
        self.w.MplWidget.canvas.axes.spines['left'].set_color(fg)
        self.w.MplWidget.canvas.axes.spines['right'].set_color(fg)
        self.w.MplWidget2.canvas.figure.set_facecolor(bg)
        self.w.MplWidget2.canvas.axes.set_facecolor(bg)
        self.w.MplWidget2.canvas.axes.xaxis.label.set_color(fg)
        self.w.MplWidget2.canvas.axes.yaxis.label.set_color(fg)
        self.w.MplWidget2.canvas.axes.tick_params(axis='x', colors=fg)
        self.w.MplWidget2.canvas.axes.tick_params(axis='y', colors=fg)
        self.w.MplWidget2.canvas.axes.spines['bottom'].set_color(fg)
        self.w.MplWidget2.canvas.axes.spines['top'].set_color(fg)
        self.w.MplWidget2.canvas.axes.spines['left'].set_color(fg)
        self.w.MplWidget2.canvas.axes.spines['right'].set_color(fg)
        self.w.MplWidget3.canvas.figure.set_facecolor(bg)
        self.w.MplWidget3.canvas.axes.set_facecolor(bg)
        self.w.MplWidget3.canvas.axes.xaxis.label.set_color(fg)
        self.w.MplWidget3.canvas.axes.yaxis.label.set_color(fg)
        self.w.MplWidget3.canvas.axes.tick_params(axis='x', colors=fg)
        self.w.MplWidget3.canvas.axes.tick_params(axis='y', colors=fg)
        self.w.MplWidget3.canvas.axes.spines['bottom'].set_color(fg)
        self.w.MplWidget3.canvas.axes.spines['top'].set_color(fg)
        self.w.MplWidget3.canvas.axes.spines['left'].set_color(fg)
        self.w.MplWidget3.canvas.axes.spines['right'].set_color(fg)
        # end load_light_mode

    def next_command(self):
        if (self.w.cmdLine.hasFocus() == True):
            if (self.nd.cmd_idx < len(self.nd.cmd_buffer)):
                self.nd.cmd_idx += 1
                if (self.nd.cmd_idx == len(self.nd.cmd_buffer)):
                    self.w.cmdLine.setText("")
                else:
                    self.w.cmdLine.setText(self.nd.cmd_buffer[self.nd.cmd_idx])

        # end next_command

    def next_tab(self):
        cidx = self.w.nmrSpectrum.currentIndex()
        while self.w.nmrSpectrum.isTabEnabled(cidx + 1) is False and cidx < 11:
            cidx += 1

        if cidx < 12:
            self.w.nmrSpectrum.setCurrentIndex(cidx + 1)
        else:
            self.w.nmrSpectrum.setCurrentIndex(0)

        # end next_tab

    def previous_tab(self):
        cidx = self.w.nmrSpectrum.currentIndex()
        while self.w.nmrSpectrum.isTabEnabled(cidx - 1) is False and cidx > 1:
            cidx -= 1

        if cidx > 0:
            self.w.nmrSpectrum.setCurrentIndex(cidx - 1)
            self.w.nmrSpectrum.setFocus()
        else:
            self.w.nmrSpectrum.setCurrentIndex(12)

        # end previous_tab

    def on_ph_corr_click(self, event):
        s = self.nd.s
        e = self.nd.e
        if (self.zoom == False):
            self.ph_corr.spc = self.nd.nmrdat[s][e].spc
            self.ph_corr.spc_max = max(max(abs(self.ph_corr.spc)))
            # self.w.MplWidget.canvas.toolbar._zoom_mode.__init__()
            if (event.button == 1):
                mods = QApplication.queryKeyboardModifiers()
                if (mods == QtCore.Qt.ControlModifier):
                    # set pivot for phase correction
                    self.ph_corr.start = event.xdata
                    self.ph_corr.pivot = event.xdata
                    self.ph_corr.piv_points = self.nd.nmrdat[s][e].ppm2points(self.ph_corr.pivot, 0)

                if (mods == QtCore.Qt.ShiftModifier):
                    # first order phase correction
                    self.ph_corr.start = event.ydata

                if (mods == QtCore.Qt.NoModifier):
                    # zero order phase correction
                    self.ph_corr.start = event.ydata

                if (mods == QtCore.Qt.AltModifier):
                    self.w.MplWidget.canvas.manager.toolbar.zoom()

            else:
                if (event.button == 2):
                    # set pivot for phase correction
                    self.ph_corr.start = event.xdata
                    self.ph_corr.pivot = event.xdata
                    self.ph_corr.piv_points = self.nd.nmrdat[s][e].ppm2points(self.ph_corr.pivot, 0)
                else:
                    # first order phase correction
                    self.ph_corr.start = event.ydata

            cid3 = self.w.MplWidget.canvas.mpl_connect('motion_notify_event', self.on_ph_corr_draw)

        # end on_ph_corr_click

    def on_ph_corr_click_2d(self, event):
        s = self.nd.s
        e = self.nd.e
        if (self.zoom == False):
            self.ph_corr.spc2 = np.copy(self.ph_corr.spc)
            # self.ph_corr.spc = self.nd.nmrdat[s][e].spc
            # self.ph_corr.spc_max = max(max(abs(self.ph_corr.spc)))
            if (event.button == 1):
                mods = QApplication.queryKeyboardModifiers()
                if (mods == QtCore.Qt.ControlModifier):
                    # set pivot for phase correction
                    self.ph_corr.start = event.xdata
                    self.ph_corr.pivot = event.xdata
                    self.ph_corr.pivot_points2d[self.ph_corr.dim] = self.nd.nmrdat[s][e].ppm2points(
                        self.ph_corr.pivot2d[self.ph_corr.dim], self.ph_corr.dim)

                if (mods == QtCore.Qt.ShiftModifier):
                    # first order phase correction
                    self.ph_corr.start = event.ydata

                if (mods == QtCore.Qt.NoModifier):
                    # zero order phase correction
                    self.ph_corr.start = event.ydata

                if (mods == QtCore.Qt.AltModifier):
                    self.w.MplWidget.canvas.manager.toolbar.zoom()

            else:
                if (event.button == 2):
                    # set pivot for phase correction
                    self.ph_corr.start = event.xdata
                    self.ph_corr.pivot2d[self.ph_corr.dim] = event.xdata
                    self.ph_corr.pivot_points2d[self.ph_corr.dim] = self.nd.nmrdat[s][e].ppm2points(
                        self.ph_corr.pivot2d[self.ph_corr.dim], self.ph_corr.dim)
                else:
                    # first order phase correction
                    self.ph_corr.start = event.ydata

            cid3 = self.w.MplWidget.canvas.mpl_connect('motion_notify_event', self.on_ph_corr_draw2d)

        # end on_ph_corr_click_2d

    def on_ph_corr_draw(self, event):
        if (self.zoom == False):
            s = self.nd.s
            e = self.nd.e
            if ((event.xdata != None) & (event.ydata != None)):
                self.ph_corr.x_data = event.xdata
                self.ph_corr.y_data = event.ydata
                if (event.button == 1):
                    mods = QApplication.queryKeyboardModifiers()
                    if (mods == QtCore.Qt.ControlModifier):
                        # set pivot for phase correction
                        self.ph_corr.pivot = event.xdata
                        self.ph_corr.piv_points = self.nd.nmrdat[s][e].ppm2points(self.ph_corr.pivot, 0)

                    if (mods == QtCore.Qt.ShiftModifier):
                        # first order phase correction
                        ph0 = 0
                        ph1 = self.ph_corr.max_ph1 * (event.ydata - self.ph_corr.start) / self.ph_corr.spc_max
                        self.ph_corr.spc = self.phase1d(self.nd.nmrdat[s][e].spc, ph0, ph1, self.ph_corr.piv_points)

                    if (mods == QtCore.Qt.NoModifier):
                        # zero order phase correction
                        ph0 = self.ph_corr.max_ph0 * (event.ydata - self.ph_corr.start) / self.ph_corr.spc_max
                        ph1 = 0
                        self.ph_corr.spc = self.phase1d(self.nd.nmrdat[s][e].spc, ph0, ph1, self.ph_corr.piv_points)

                else:
                    if (event.button == 2):
                        # set pivot for phase correction
                        self.ph_corr.x_data = event.xdata
                        self.ph_corr.y_data = event.ydata
                        self.ph_corr.pivot = event.xdata
                        self.ph_corr.piv_points = self.nd.nmrdat[s][e].ppm2points(self.ph_corr.pivot, 0)
                    else:
                        # first order phase correction
                        self.ph_corr.x_data = event.xdata
                        self.ph_corr.y_data = event.ydata
                        ph0 = 0
                        ph1 = self.ph_corr.max_ph1 * (event.ydata - self.ph_corr.start) / self.ph_corr.spc_max
                        self.ph_corr.spc = self.phase1d(self.nd.nmrdat[s][e].spc, ph0, ph1, self.ph_corr.piv_points)

            self.ph_corr_plot_spc()

        # end on_ph_corr_draw

    def on_ph_corr_draw2d(self, event):
        if (self.zoom == False):
            s = self.nd.s
            e = self.nd.e
            if ((event.xdata != None) & (event.ydata != None)):
                self.ph_corr.x_data = event.xdata
                self.ph_corr.y_data = event.ydata
                if (event.button == 1):
                    mods = QApplication.queryKeyboardModifiers()
                    if (mods == QtCore.Qt.ControlModifier):
                        # set pivot for phase correction
                        self.ph_corr.pivot2d[self.ph_corr.dim] = event.xdata
                        self.ph_corr.pivot_points2d[self.ph_corr.dim] = self.nd.nmrdat[s][e].ppm2points(
                            self.ph_corr.pivot2d[self.ph_corr.dim], self.ph_corr.dim)

                    if (mods == QtCore.Qt.ShiftModifier):
                        # first order phase correction
                        ph0 = 0
                        ph1 = self.ph_corr.max_ph1 * (event.ydata - self.ph_corr.start) / self.ph_corr.spc_max
                        self.ph_corr.spc = self.phase1d(self.ph_corr.spc2, ph0, ph1,
                                                        self.ph_corr.pivot_points2d[self.ph_corr.dim])

                    if (mods == QtCore.Qt.NoModifier):
                        # zero order phase correction
                        ph0 = self.ph_corr.max_ph0 * (event.ydata - self.ph_corr.start) / self.ph_corr.spc_max
                        ph1 = 0
                        self.ph_corr.spc = self.phase1d(self.ph_corr.spc2, ph0, ph1,
                                                        self.ph_corr.pivot_points2d[self.ph_corr.dim])

                else:
                    if (event.button == 2):
                        # set pivot for phase correction
                        self.ph_corr.x_data = event.xdata
                        self.ph_corr.y_data = event.ydata
                        self.ph_corr.pivot2d[self.ph_corr.dim] = event.xdata
                        self.ph_corr.pivot_points2d[self.ph_corr.dim] = self.nd.nmrdat[s][e].ppm2points(
                            self.ph_corr.pivot2d[self.ph_corr.dim], self.ph_corr.dim)
                    else:
                        # first order phase correction
                        self.ph_corr.x_data = event.xdata
                        self.ph_corr.y_data = event.ydata
                        ph0 = 0
                        ph1 = self.ph_corr.max_ph1 * (event.ydata - self.ph_corr.start) / self.ph_corr.spc_max
                        self.ph_corr.spc = self.phase1d(self.ph_corr.spc2, ph0, ph1,
                                                        self.ph_corr.pivot_points2d[self.ph_corr.dim])

            self.ph_corr_plot_spc_2d()

        # end on_ph_corr_drawd

    def on_ph_corr_release(self, event):
        s = self.nd.s
        e = self.nd.e
        if ((event.xdata != None) & (event.ydata != None)):
            xdata = event.xdata
            ydata = event.ydata
        else:
            xdata = self.ph_corr.x_data
            ydata = self.ph_corr.y_data

        if (self.zoom == False):
            if (event.button == 1):
                mods = QApplication.queryKeyboardModifiers()
                if (mods == QtCore.Qt.ControlModifier):
                    # set pivot for phase correction
                    self.ph_corr.pivot = xdata
                    self.ph_corr.piv_points = self.nd.nmrdat[s][e].ppm2points(self.ph_corr.pivot, 0)

                if (mods == QtCore.Qt.ShiftModifier):
                    # first order phase correction
                    ph1 = (self.ph_corr.max_ph1 * (ydata - self.ph_corr.start) / self.ph_corr.spc_max)
                    ph = self.phases_remove_pivot(0.0, ph1, self.ph_corr.piv_points, len(self.ph_corr.spc[0]))
                    ph0 = ((self.nd.nmrdat[s][e].proc.ph0[0] + ph[0] + 180.0) % 360.0) - 180.0
                    ph1 = self.nd.nmrdat[s][e].proc.ph1[0] + ph[1]
                    self.nd.nmrdat[s][e].proc.ph0[0] = ph0
                    self.nd.nmrdat[s][e].proc.ph1[0] = ph1

                if (mods == QtCore.Qt.NoModifier):
                    # zero order phase correction
                    ph0a = (self.ph_corr.max_ph0 * (ydata - self.ph_corr.start) / self.ph_corr.spc_max) % 360.0
                    ph1a = 0.0
                    ph = self.phases_remove_pivot(ph0a, ph1a, self.ph_corr.piv_points, len(self.ph_corr.spc[0]))
                    ph0 = ((self.nd.nmrdat[s][e].proc.ph0[0] + ph[0] + 180.0) % 360.0) - 180.0
                    ph1 = self.nd.nmrdat[s][e].proc.ph1[0] + ph[1]
                    self.nd.nmrdat[s][e].proc.ph0[0] = ph0
                    self.nd.nmrdat[s][e].proc.ph1[0] = ph1

            else:
                if (event.button == 2):
                    # set pivot for phase correction
                    self.ph_corr.pivot = xdata
                    self.ph_corr.piv_points = self.nd.nmrdat[s][e].ppm2points(self.ph_corr.pivot, 0)
                else:
                    # first order phase correction
                    ph1 = (self.ph_corr.max_ph1 * (ydata - self.ph_corr.start) / self.ph_corr.spc_max)
                    ph = self.phases_remove_pivot(0.0, ph1, self.ph_corr.piv_points, len(self.ph_corr.spc[0]))
                    ph0 = ((self.nd.nmrdat[s][e].proc.ph0[0] + ph[0] + 180.0) % 360.0) - 180.0
                    ph1 = self.nd.nmrdat[s][e].proc.ph1[0] + ph[1]
                    self.nd.nmrdat[s][e].proc.ph0[0] = ph0
                    self.nd.nmrdat[s][e].proc.ph1[0] = ph1

            cid3 = self.w.MplWidget.canvas.mpl_connect('motion_notify_event', self.on_ph_corr_draw)
            cid3 = self.w.MplWidget.canvas.mpl_disconnect(cid3)
            self.nd.nmrdat[s][e].spc = self.ph_corr.spc
            self.set_proc_pars()
            self.nd.ft()
            self.ph_corr_plot_spc()
        else:
            # zoom mode activated
            if (event.button > 1):
                # Right MB click will unzoom the plot
                try:
                    self.w.MplWidget.canvas.figure.canvas.toolbar.home()
                except:
                    pass

        # end on_ph_corr_release

    def on_ph_corr_release_2d(self, event):
        s = self.nd.s
        e = self.nd.e
        if ((event.xdata != None) & (event.ydata != None)):
            xdata = event.xdata
            ydata = event.ydata
        else:
            xdata = self.ph_corr.x_data
            ydata = self.ph_corr.y_data

        if (self.zoom == False):
            if event.button == 1:
                mods = QApplication.queryKeyboardModifiers()
                if mods == QtCore.Qt.ControlModifier:
                    # set pivot for phase correction
                    self.ph_corr.pivot2d[self.ph_corr.dim] = xdata
                    self.ph_corr.pivot_points2d[self.ph_corr.dim] = self.nd.nmrdat[s][e].ppm2points(
                        self.ph_corr.pivot2d[self.ph_corr.dim], self.ph_corr.dim)

                if mods == QtCore.Qt.ShiftModifier:
                    # first order phase correction
                    ph1 = (self.ph_corr.max_ph1 * (ydata - self.ph_corr.start) / self.ph_corr.spc_max)
                    ph = self.phases_remove_pivot(0.0, ph1, self.ph_corr.pivot_points2d[self.ph_corr.dim],
                                                  len(self.ph_corr.spc[0]))
                    ph0 = ((self.ph_corr.ph0_2d[self.ph_corr.dim] + ph[0] + 180.0) % 360.0) - 180.0
                    ph1 = self.ph_corr.ph1_2d[self.ph_corr.dim] + ph[1]
                    self.ph_corr.ph0_2d[self.ph_corr.dim] = ph0
                    self.ph_corr.ph1_2d[self.ph_corr.dim] = ph1

                if mods == QtCore.Qt.NoModifier:
                    # zero order phase correction
                    ph0a = (self.ph_corr.max_ph0 * (ydata - self.ph_corr.start) / self.ph_corr.spc_max) % 360.0
                    ph1a = 0.0
                    ph = self.phases_remove_pivot(ph0a, ph1a, self.ph_corr.pivot_points2d[self.ph_corr.dim],
                                                  len(self.ph_corr.spc[0]))
                    ph0 = ((self.ph_corr.ph0_2d[self.ph_corr.dim] + ph[0] + 180.0) % 360.0) - 180.0
                    ph1 = self.ph_corr.ph1_2d[self.ph_corr.dim] + ph[1]
                    self.ph_corr.ph0_2d[self.ph_corr.dim] = ph0
                    self.ph_corr.ph1_2d[self.ph_corr.dim] = ph1

            else:
                if event.button == 2:
                    # set pivot for phase correction
                    self.ph_corr.pivot2d[self.ph_corr.dim] = xdata
                    self.ph_corr.pivot_points2d[self.ph_corr.dim] = self.nd.nmrdat[s][e].ppm2points(
                        self.ph_corr.pivot2d[self.ph_corr.dim], self.ph_corr.dim)

                else:
                    # first order phase correction
                    ph1 = (self.ph_corr.max_ph1 * (ydata - self.ph_corr.start) / self.ph_corr.spc_max)
                    ph = self.phases_remove_pivot(0.0, ph1, self.ph_corr.pivot_points2d[self.ph_corr.dim],
                                                  len(self.ph_corr.spc[0]))
                    ph0 = ((self.ph_corr.ph0_2d[self.ph_corr.dim] + ph[0] + 180.0) % 360.0) - 180.0
                    ph1 = self.ph_corr.ph1_2d[self.ph_corr.dim] + ph[1]
                    self.ph_corr.ph0_2d[self.ph_corr.dim] = ph0
                    self.ph_corr.ph1_2d[self.ph_corr.dim] = ph1

            cid3 = self.w.MplWidget.canvas.mpl_connect('motion_notify_event', self.on_ph_corr_draw2d)
            cid3 = self.w.MplWidget.canvas.mpl_disconnect(cid3)
            self.ph_corr.spc2 = np.copy(self.ph_corr.spc)
            self.ph_corr_plot_spc_2d()
        else:
            # zoom mode activated
            if (event.button > 1):
                # Right MB click will unzoom the plot
                try:
                    self.w.MplWidget.canvas.figure.canvas.toolbar.home()
                except:
                    pass

        # end on_ph_corr_release_2d

    def open_script(self, f_name=""):
        if (f_name == False):
            f_name = ""

        if (len(f_name) == 0):
            f_name = QFileDialog.getOpenFileName(None, 'Open Script File', '', 'Python files (*.py)')
            f_name = f_name[0]

        if (len(f_name) > 0):
            f = open(f_name, 'r')
            scriptText = f.read()
            self.w.script.setText(scriptText)

        self.w.nmrSpectrum.setCurrentIndex(10)
        # end open_script

    def ph_corr_plot_spc(self):
        xlim = self.w.MplWidget.canvas.axes.get_xlim()
        ylim = self.w.MplWidget.canvas.axes.get_ylim()
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        if (d.pos_col == "RGB"):
            pos_col = d.pos_col_rgb
        else:
            pos_col = d.pos_col

        if (d.neg_col == "RGB"):
            neg_col = d.neg_col_rgb
        else:
            neg_col = d.neg_col

        ref_col = d.ph_ref_col
        pos_col = matplotlib.colors.to_hex(pos_col)
        neg_col = matplotlib.colors.to_hex(neg_col)
        ref_col = matplotlib.colors.to_hex(ref_col)
        xlabel = d.x_label + " [" + d.axis_type1 + "]"
        ylabel = d.y_label + " [" + d.axis_type2 + "]"
        if (self.nd.nmrdat[self.nd.s][self.nd.e].dim == 1):
            self.w.MplWidget.canvas.axes.clear()
            if ((d.ph_ref_ds > 0) & (d.ph_ref_exp > 0) & (
                    ((d.ph_ref_ds - 1 == self.nd.s) & (d.ph_ref_exp - 1 == self.nd.e)) is False)):
                self.w.MplWidget.canvas.axes.plot(self.nd.nmrdat[d.ph_ref_ds - 1][d.ph_ref_exp - 1].ppm1,
                                                  self.nd.nmrdat[d.ph_ref_ds - 1][d.ph_ref_exp - 1].spc[0].real,
                                                  color=ref_col)

            self.w.MplWidget.canvas.axes.plot(self.nd.nmrdat[self.nd.s][self.nd.e].ppm1, self.ph_corr.spc[0].real,
                                              color=pos_col)
            self.w.MplWidget.canvas.axes.plot([self.ph_corr.pivot, self.ph_corr.pivot],
                                              [2.0 * self.ph_corr.spc_max, -2.0 * self.ph_corr.spc_max], color='r')
            self.w.MplWidget.canvas.axes.set_xlabel(xlabel)
            self.w.MplWidget.canvas.axes.invert_xaxis()
            self.w.MplWidget.canvas.axes.set_xlim(xlim)
            self.w.MplWidget.canvas.axes.set_ylim(ylim)

        self.set_proc_pars()
        self.w.MplWidget.canvas.draw()
        # This is a messy solution to force the matplotlib widget to update the plot by introducing an error (calling
        # a figure object and redirecting the error output
        code_err = io.StringIO()
        sys.stderr = code_err
        try:
            self.w.MplWidget.canvas.figure()
        except:
            pass

        sys.stderr = sys.__stderr__
        # end ph_corr_plot_spc

    def ph_corr_plot_spc_2d(self, keep_zoom=True):
        if keep_zoom:
            xlim = self.w.MplWidget.canvas.axes.get_xlim()
            ylim = self.w.MplWidget.canvas.axes.get_ylim()

        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        self.w.MplWidget.canvas.axes.set_prop_cycle(None)
        if self.ph_corr.dim == 0:
            xlabel = d.x_label + " [" + d.axis_type1 + "]"
        else:
            xlabel = d.y_label + " [" + d.axis_type2 + "]"

        self.w.MplWidget.canvas.axes.clear()
        self.ph_corr.spc_max = 0.0
        for k in range(len(self.ph_corr.spc)):
            self.w.MplWidget.canvas.axes.plot(self.ph_corr.ppm, self.ph_corr.spc[k].real)
            self.ph_corr.spc_max = max(self.ph_corr.spc_max, np.max(np.max(np.abs(self.ph_corr.spc[k].real))))

        self.w.MplWidget.canvas.axes.invert_xaxis()
        if not keep_zoom:
            xlim = self.w.MplWidget.canvas.axes.get_xlim()
            ylim = self.w.MplWidget.canvas.axes.get_ylim()

        self.w.MplWidget.canvas.axes.plot(
            [self.ph_corr.pivot2d[self.ph_corr.dim], self.ph_corr.pivot2d[self.ph_corr.dim]],
            [2.0 * self.ph_corr.spc_max, -2.0 * self.ph_corr.spc_max], color='r')
        self.w.MplWidget.canvas.axes.set_xlabel(xlabel)
        self.w.MplWidget.canvas.axes.invert_xaxis()
        self.w.MplWidget.canvas.axes.set_xlim(xlim)
        self.w.MplWidget.canvas.axes.set_ylim(ylim)
        self.w.MplWidget.canvas.draw()
        # This is a messy solution to force the matplotlib widget to update the plot by introducing an error (calling
        # a figure object and redirecting the error output
        code_err = io.StringIO()
        sys.stderr = code_err
        try:
            self.w.MplWidget.canvas.figure()
        except:
            pass

        sys.stderr = sys.__stderr__
        # end ph_corr_plot_spc_2d

    def phase1d(self, mat, ph0, ph1, piv):
        npts = len(mat[0])
        ph0 = -ph0 * math.pi / 180.0
        ph1 = -ph1 * math.pi / 180.0
        frac = np.linspace(0, 1, npts) - float(npts - piv) / float(npts)
        ph = ph0 + frac * ph1
        mat = np.cos(ph) * mat.real + np.sin(ph) * mat.imag + 1j * (-np.sin(ph) * mat.real + np.cos(ph) * mat.imag)
        return mat
        # end phase1d

    def phases_remove_pivot(self, phc0, phc1, piv, npts):
        phases = np.array([0.0, 0.0])
        frac = np.linspace(0, 1, npts) - float(npts - piv) / float(npts)
        ph = -phc0 - frac * phc1
        phases[0] = -ph[0]
        phases[1] = ph[0] - ph[len(ph) - 1]
        return phases
        # end phases_remove_pivot

    def pick_col_row(self):
        self.w.statusBar().clearMessage()
        self.w.statusBar().showMessage("Click to add row/col")
        self.show_acquisition_parameters()
        self.show_nmr_spectrum()
        self.ginput2d()
        # end pick_col_row

    def plot_2d_col_row(self):
        while len(self.w.MplWidget.canvas.axes.lines) > 0:
            self.w.MplWidget.canvas.axes.lines[0].remove()

        self.w.MplWidget.canvas.axes.set_prop_cycle(None)
        ppm1 = self.nd.nmrdat[self.nd.s][self.nd.e].ppm1
        ppm2 = self.nd.nmrdat[self.nd.s][self.nd.e].ppm2
        for k in range(len(self.ph_corr.spc_row)):
            pid = self.w.MplWidget.canvas.axes.plot([self.ph_corr.spc_col[k], self.ph_corr.spc_col[k]],
                                                    [np.min(ppm2), np.max(ppm2)])
            self.w.MplWidget.canvas.axes.plot([np.min(ppm1), np.max(ppm1)],
                                              [self.ph_corr.spc_row[k], self.ph_corr.spc_row[k]],
                                              color=pid[0].get_color())

        self.w.MplWidget.canvas.draw()

    def plot_spc(self, hide_pre_processing=False):
        self.keep_zoom = self.w.keepZoom.isChecked()
        xlim = self.w.MplWidget.canvas.axes.get_xlim()
        ylim = self.w.MplWidget.canvas.axes.get_ylim()
        self.w.nmrSpectrum.setCurrentIndex(0)
        if (len(self.nd.nmrdat[self.nd.s]) == 0):
            return

        if (len(self.nd.nmrdat[self.nd.s][self.nd.e].spc) == 0):
            return

        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        if d.pos_col == "RGB":
            pos_col = d.pos_col_rgb
        else:
            pos_col = d.pos_col

        if d.neg_col == "RGB":
            neg_col = d.neg_col_rgb
        else:
            neg_col = d.neg_col

        pos_col = matplotlib.colors.to_hex(pos_col)
        neg_col = matplotlib.colors.to_hex(neg_col)
        xlabel = d.x_label + " [" + d.axis_type1 + "]"
        ylabel = d.y_label + " [" + d.axis_type2 + "]"
        # print(self.nd.nmrdat[self.nd.s][self.nd.e].dim)
        if (self.nd.nmrdat[self.nd.s][self.nd.e].dim == 1):
            self.w.MplWidget.canvas.axes.clear()
            for k in range(len(self.nd.nmrdat[self.nd.s])):
                if ((k != self.nd.e) and (self.nd.nmrdat[self.nd.s][k].display.display_spc == True)):
                    d = self.nd.nmrdat[self.nd.s][k].display
                    if (d.pos_col == "RGB"):
                        pos_col = d.pos_col_rgb
                    else:
                        pos_col = d.pos_col

                    if (d.neg_col == "RGB"):
                        neg_col = d.neg_col_rgb
                    else:
                        neg_col = d.neg_col

                    pos_col = matplotlib.colors.to_hex(pos_col)
                    neg_col = matplotlib.colors.to_hex(neg_col)
                    self.w.MplWidget.canvas.axes.plot(self.nd.nmrdat[self.nd.s][k].ppm1,
                                                      self.nd.nmrdat[self.nd.s][k].spc[0].real, color=pos_col)

            d = self.nd.nmrdat[self.nd.s][self.nd.e].display
            if (d.pos_col == "RGB"):
                pos_col = d.pos_col_rgb
            else:
                pos_col = d.pos_col

            if (d.neg_col == "RGB"):
                neg_col = d.neg_col_rgb
            else:
                neg_col = d.neg_col

            pos_col = matplotlib.colors.to_hex(pos_col)
            neg_col = matplotlib.colors.to_hex(neg_col)
            xlabel = d.x_label + " [" + d.axis_type1 + "]"
            ylabel = d.y_label + " [" + d.axis_type2 + "]"
            if len(self.nd.nmrdat[self.nd.s][self.nd.e].start_peak) > 0:
                if self.w.peakPicking.isChecked() == True:
                    s = self.nd.s
                    e = self.nd.e
                    for k in range(len(self.nd.nmrdat[s][e].start_peak)):
                        self.w.MplWidget.canvas.axes.axvspan(self.nd.nmrdat[s][e].start_peak[k],
                                                             self.nd.nmrdat[s][e].end_peak[k],
                                                             alpha=self.nd.pp.alpha, color=self.nd.pp.colour)

            self.w.MplWidget.canvas.axes.plot(self.nd.nmrdat[self.nd.s][self.nd.e].ppm1,
                                              self.nd.nmrdat[self.nd.s][self.nd.e].spc[0].real, color=pos_col)
            self.w.MplWidget.canvas.axes.set_xlabel(xlabel)
            self.w.MplWidget.canvas.axes.autoscale()
            self.w.MplWidget.canvas.axes.invert_xaxis()
            if (self.keep_zoom == True):
                self.w.MplWidget.canvas.axes.set_xlim(xlim)
                self.w.MplWidget.canvas.axes.set_ylim(ylim)

            # self.w.MplWidget.canvas.toolbar.update()
            self.w.MplWidget.canvas.draw()
            if (self.keep_x_zoom == True):
                self.w.MplWidget.canvas.axes.set_xlim(xlim)
                self.vertical_auto_scale()
                self.keep_x_zoom = False


        else:
            mm = np.max(np.abs(self.nd.nmrdat[self.nd.s][self.nd.e].spc.real))
            pos_lev = np.linspace(d.min_level * mm, d.max_level * mm, d.n_levels)
            neg_lev = np.linspace(-d.max_level * mm, -d.min_level * mm, d.n_levels)
            self.w.MplWidget.canvas.axes.clear()
            self.w.MplWidget.canvas.axes.contour(self.nd.nmrdat[self.nd.s][self.nd.e].ppm1,
                                                 self.nd.nmrdat[self.nd.s][self.nd.e].ppm2,
                                                 self.nd.nmrdat[self.nd.s][self.nd.e].spc.real, pos_lev, colors=pos_col,
                                                 linestyles='solid', antialiased=True)
            self.w.MplWidget.canvas.axes.contour(self.nd.nmrdat[self.nd.s][self.nd.e].ppm1,
                                                 self.nd.nmrdat[self.nd.s][self.nd.e].ppm2,
                                                 self.nd.nmrdat[self.nd.s][self.nd.e].spc.real, neg_lev, colors=neg_col,
                                                 linestyles='solid', antialiased=True)
            self.w.MplWidget.canvas.axes.set_xlabel(xlabel)
            self.w.MplWidget.canvas.axes.set_ylabel(ylabel)
            self.w.MplWidget.canvas.axes.autoscale()
            self.w.MplWidget.canvas.axes.invert_xaxis()
            self.w.MplWidget.canvas.axes.invert_yaxis()
            if (self.keep_zoom == True):
                self.w.MplWidget.canvas.axes.set_xlim(xlim)
                self.w.MplWidget.canvas.axes.set_ylim(ylim)
            else:
                if (self.keep_x_zoom == True):
                    self.w.MplWidget.canvas.axes.set_xlim(xlim)
                    self.keep_x_zoom = False

            # self.w.MplWidget.canvas.toolbar.update()
            self.w.MplWidget.canvas.draw()

        self.keep_zoom = False
        if hide_pre_processing == False and self.w.peakPicking.isChecked() == False:
            if self.exited_peak_picking == False:
                #a = 3
                pyautogui.click(clicks=1)
            else:
                self.exited_peak_picking = True

        # end plot_spc

    def plot_spc_disp(self):
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.change_data_set_exp()
        if (self.ph_corr_active == False):
            self.plot_spc()
        else:
            self.ph_corr_plot_spc()

        # end plot_spc_disp

    def plot_spc_pre_proc(self):
        if (len(self.nd.pp.class_select) == 0):
            self.nd.pre_proc_init()

        # self.w.rDolphinExport.setChecked(self.nd.pp.rDolphinExport)
        self.fill_pre_processing_numbers()
        sel = self.w.selectClassTW.selectedIndexes()
        cls = np.array([])
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            cls = np.append(cls, self.w.selectClassTW.item(k, 1).text())

        self.nd.pp.class_select = cls
        cls2 = np.unique(cls)
        sel2 = np.array([], dtype='int')
        for k in range(len(sel)):
            if (sel[k].column() == 0):
                sel2 = np.append(sel2, int(sel[k].row()))

        self.nd.pp.plot_select = sel2
        self.keep_zoom = self.w.keepZoom.isChecked()
        xlim = self.w.MplWidget.canvas.axes.get_xlim()
        ylim = self.w.MplWidget.canvas.axes.get_ylim()
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.w.MplWidget.canvas.axes.clear()
        if (self.w.preProcessingWidget.currentIndex() == 1):
            for k in range(len(self.nd.pp.exclude_start)):
                self.w.MplWidget.canvas.axes.axvspan(self.nd.pp.exclude_start[k], self.nd.pp.exclude_end[k],
                                                     alpha=self.nd.pp.alpha, color=self.nd.pp.colour)

        if (self.w.preProcessingWidget.currentIndex() == 2):
            for k in range(len(self.nd.pp.seg_start)):
                self.w.MplWidget.canvas.axes.axvspan(self.nd.pp.seg_start[k], self.nd.pp.seg_end[k],
                                                     alpha=self.nd.pp.alpha, color=self.nd.pp.colour)

        for k in range(len(self.nd.pp.plot_select)):
            colIdx = np.where(cls2 == cls[self.nd.pp.plot_select[k]])[0][0]
            plotCol = matplotlib.colors.to_hex(self.nd.pp.plot_colours[colIdx])
            self.w.MplWidget.canvas.axes.plot(self.nd.nmrdat[self.nd.s][self.nd.pp.plot_select[k]].ppm1,
                                              self.nd.nmrdat[self.nd.s][self.nd.pp.plot_select[k]].spc[0].real,
                                              color=plotCol)

        if (self.w.preProcessingWidget.currentIndex() == 3):
            self.w.MplWidget.canvas.axes.axvspan(self.nd.pp.noise_start, self.nd.pp.noise_end, alpha=self.nd.pp.alpha,
                                                 color=self.nd.pp.colour)
            val = self.nd.pp.noise_threshold * self.nd.pp.std_val
            x = [self.nd.nmrdat[self.nd.s][0].ppm1[0], self.nd.nmrdat[self.nd.s][0].ppm1[-1]]
            y = [val, val]
            self.w.MplWidget.canvas.axes.plot(x, y, color=self.nd.pp.th_colour, linewidth=self.nd.pp.th_line_width)

        if (self.w.preProcessingWidget.currentIndex() == 5):
            for k in range(len(self.nd.pp.compress_start)):
                self.w.MplWidget.canvas.axes.axvspan(self.nd.pp.compress_start[k], self.nd.pp.compress_end[k],
                                                     alpha=self.nd.pp.alpha, color=self.nd.pp.colour)

        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        xlabel = d.x_label + " [" + d.axis_type1 + "]"
        self.w.MplWidget.canvas.axes.set_xlabel(xlabel)
        self.w.MplWidget.canvas.axes.autoscale()
        self.w.MplWidget.canvas.axes.invert_xaxis()
        if (self.keep_zoom == True):
            self.w.MplWidget.canvas.axes.set_xlim(xlim)
            self.w.MplWidget.canvas.axes.set_ylim(ylim)

        # self.w.MplWidget.canvas.toolbar.update()
        self.w.MplWidget.canvas.draw()

    def previous_command(self):
        if (self.w.cmdLine.hasFocus() == True):
            if (self.nd.cmd_idx > 0):
                self.nd.cmd_idx -= 1
                self.w.cmdLine.setText(self.nd.cmd_buffer[self.nd.cmd_idx])

        # end previous_command

    def quit_app(self):
        # some actions to perform before actually quitting:
        try:
            self.p.terminate()
            sleep(2)
        except:
            pass

        self.w.close()
        # end quit_app

    def read_nmr_spc(self):
        kz = self.w.keepZoom.isChecked()
        if (len(self.nd.nmrdat[0]) == 0):
            self.w.keepZoom.setChecked(False)

        selected_directory = QFileDialog.getExistingDirectory()
        if (len(selected_directory) > 0):
            # Use the selected directory...
            idx = selected_directory.rfind('/')
            ds_name = selected_directory[:idx]
            exp_name = selected_directory[idx + 1:]
            self.nd.read_spc(ds_name, exp_name)
            self.set_j_res()
            self.nd.ft()
            if self.nd.nmrdat[self.nd.s][self.nd.e].dim == 0:
                self.nd.auto_ref(True)
            else:
                self.nd.auto_ref(False)

            self.nd.e = len(self.nd.nmrdat[self.nd.s]) - 1
            self.plot_spc()
            self.w.keepZoom.setChecked(kz)
            self.set_proc_pars()
            self.set_acq_pars()
            self.set_title_file()
            self.set_pulse_program()
            self.w.expBox.setValue(self.nd.e + 1)
            self.set_disp_pars()
            self.update_gui()

        # end read_nmr_spc

    def read_nmrpipe_spc(self, sfile=False):
        if sfile == False:
            selected_file = QFileDialog.getOpenFileName()
            if len(selected_file[0]) == 0:
                return

        else:
            selected_file = (sfile, '')

        # print(selected_file)
        f_name = os.path.split(selected_file[0])[1]
        data_path = os.path.split(os.path.split(selected_file[0])[0])[0]
        exp_num = os.path.split(os.path.split(selected_file[0])[0])[1]
        if exp_num.find('.') > -1:
            exp_num = exp_num[:exp_num.find('.')]

        self.read_nmrpipe_spcs([data_path], [exp_num], f_name)
        self.set_standard_colours()
        self.update_gui()
        self.reset_plot()
        # end read_nmrpipe_spc

    def read_nmrpipe_spcs(self, data_path, data_sets, proc_data_name='test.dat'):
        z_fill = 25
        if (data_path[0] == 'interactive'):
            data_path = [QFileDialog.getExistingDirectory()]

        if (len(data_path) > 0):
            if (str(data_sets) == 'all'):
                folders = []
                for r, d, f in os.walk(data_path):
                    for folder in d:
                        if (os.path.isfile(os.path.join(r, folder, proc_data_name))):
                            folders.append(folder.z_fill(z_fill).rstrip('.proc'))

                folders.sort()
                data_sets = []
                for k in range(len(folders)):
                    data_sets.append(int(folders[k]))

            self.nd.read_nmrpipe_spcs(data_path, data_sets, proc_data_name)
        # end read_nmrpipe_spcs

    def read_spcs(self, data_path, data_sets):
        z_fill = 25
        if (data_path[0] == 'interactive'):
            data_path = [QFileDialog.getExistingDirectory()]

        if (len(data_path) > 0):
            if (str(data_sets[0]) == 'all'):
                folders = []
                for r, d, f in os.walk(data_path[0]):
                    for folder in d:
                        if (os.path.isfile(os.path.join(r, folder, 'fid'))):
                            if (folder != '99999'):
                                folders.append(folder.z_fill(z_fill))

                        if (os.path.isfile(os.path.join(r, folder, 'ser'))):
                            folders.append(folder.z_fill(z_fill))

                folders.sort()
                data_sets = []
                for k in range(len(folders)):
                    data_sets.append(int(folders[k]))

            if (str(data_sets[0]) == 'all1d'):
                folders = []
                for r, d, f in os.walk(data_path[0]):
                    for folder in d:
                        if (os.path.isfile(os.path.join(r, folder, 'fid'))):
                            if (folder != '99999'):
                                folders.append(folder.z_fill(z_fill))

                folders.sort()
                data_sets = []
                for k in range(len(folders)):
                    data_sets.append(int(folders[k]))

            if (str(data_sets[0]) == 'all2d'):
                folders = []
                for r, d, f in os.walk(data_path[0]):
                    for folder in d:
                        if (os.path.isfile(os.path.join(r, folder, 'ser'))):
                            folders.append(folder.z_fill(z_fill))

                folders.sort()
                data_sets = []
                for k in range(len(folders)):
                    data_sets.append(int(folders[k]))

            if len(data_path) > 1:
                dp = []
                for d in data_path:
                    if os.path.isfile(os.path.join(d, data_sets[0], 'fid')) or os.path.isfile(
                            os.path.join(d, data_sets[0], 'ser')):
                        dp.append(d)

                data_path = dp

            else:
                ds = []
                for d in data_sets:
                    if os.path.isfile(os.path.join(data_path[0], str(d), 'fid')) or os.path.isfile(
                            os.path.join(data_path[0], str(d), 'ser')):
                        ds.append(d)

                data_sets = ds

            if len(data_path) > 0 and len(data_sets) > 0:
                self.nd.read_spcs(data_path, data_sets)

        # end read_spcs

    def reference1d(self, ref_shift=0.0):
        self.temp_ref_shift = ref_shift
        self.w.MplWidget.canvas.setFocus()
        self.show_nmr_spectrum()
        self.ginput_ref_1d(1)
        # end reference1d

    def reference2d(self, ref_shift=[0.0, 0.0]):
        self.temp_ref_shift = ref_shift
        self.w.MplWidget.canvas.setFocus()
        self.show_nmr_spectrum()
        self.ginput_ref_2d(1)
        # end reference2d

    def remove_last_col_row(self):
        n_lines = len(self.w.MplWidget.canvas.axes.lines)
        if n_lines > 0:
            self.w.MplWidget.canvas.axes.lines[n_lines - 1].remove()
            self.ph_corr.spc_row = self.ph_corr.spc_row[:-1]
            self.ph_corr.spc_col = self.ph_corr.spc_col[:-1]
            self.ph_corr.spc_row_pts = self.ph_corr.spc_row_pts[:-1]
            self.ph_corr.spc_col_pts = self.ph_corr.spc_col_pts[:-1]
            self.plot_2d_col_row()
            self.show_acquisition_parameters()
            self.show_nmr_spectrum()

        # end remove_last_col_row

    def reset_config(self):
        self.cf = nmrConfig.NmrConfig()
        self.cf.save_config()
        self.load_config()
        # end reset_config

    def reset_data_pre_processing(self):
        self.nd.reset_data_pre_processing()
        self.plot_spc_pre_proc()
        self.vertical_auto_scale()
        self.w.MplWidget.canvas.flush_events()
        self.w.MplWidget.canvas.draw()
        # end data_pre_processing

    def reset_help(self):
        if self.cf.mode == 'dark':
            f_name = os.path.join(os.path.dirname(__file__), "nmr", "web", "introductionDark", "index.html")
        else:
            f_name = os.path.join(os.path.dirname(__file__), "nmr", "web", "introduction", "index.html")

        url = "file:///" + f_name.replace('\\', '/')
        self.w.helpView.setUrl(url)
        self.w.nmrSpectrum.setCurrentIndex(12)
        # end reset_help

    def reset_plot(self):
        zoom_checked = self.w.keepZoom.isChecked()
        self.w.keepZoom.setChecked(False)
        self.plot_spc()
        if (zoom_checked == True):
            self.w.keepZoom.setChecked(True)

        # end reset_plot

    def save_button(self):
        pf_name = QFileDialog.getSaveFileName(None, "Save MetaboLabPy DataSet", "", "*.mlpy", "*.mlpy")
        f_name = pf_name[0].rstrip('.mlpy').rstrip(' ').rstrip('/').rstrip('.mlpy') + '.mlpy'

        if (os.path.isfile(f_name)):
            os.remove(f_name)

        if (os.path.isdir(f_name)):
            shutil.rmtree(f_name)

        self.nd.script = self.w.script.toHtml()
        self.nd.console = self.w.console.toHtml()
        self.nd.save(f_name)
        # end save_button

    def save_config(self):
        self.cf.auto_plot = self.w.autoPlot.isChecked()
        self.cf.keep_zoom = self.w.keepZoom.isChecked()
        self.cf.font_size = self.w.fontSize.value()
        self.cf.phase_reference_colour = self.nd.nmrdat[0][0].display.ph_ref_col
        self.cf.pos_col10 = self.std_pos_col1[0]
        self.cf.pos_col11 = self.std_pos_col1[1]
        self.cf.pos_col12 = self.std_pos_col1[2]
        self.cf.neg_col10 = self.std_neg_col1[0]
        self.cf.neg_col11 = self.std_neg_col1[1]
        self.cf.neg_col12 = self.std_neg_col1[2]
        self.cf.pos_col20 = self.std_pos_col2[0]
        self.cf.pos_col21 = self.std_pos_col2[1]
        self.cf.pos_col22 = self.std_pos_col2[2]
        self.cf.neg_col20 = self.std_neg_col2[0]
        self.cf.neg_col21 = self.std_neg_col2[1]
        self.cf.neg_col22 = self.std_neg_col2[2]
        self.cf.save_config()
        # end save_config

    def save_mat(self):
        scipy.io.save_mat('/Users/ludwigc/metabolabpy.mat',
                          {'spc': self.nd.nmrdat[0][0].spc, 'fid': self.nd.nmrdat[0][0].fid})
        # end save_mat

    def save_script(self, f_name=""):
        if (f_name == False):
            f_name = ""

        if (len(f_name) == 0):
            f_name = QFileDialog.getSaveFileName(None, 'Save Script File', '', 'Python files (*.py)')
            f_name = f_name[0]

        if (len(f_name) > 0):
            scriptText = self.w.script.toPlainText()
            f = open(f_name, 'w')
            f.write(scriptText)

        # end open_script

    def scale_2d_spectrum_up(self):
        self.nd.nmrdat[self.nd.s][self.nd.e].display.min_level /= 1.1
        self.nd.nmrdat[self.nd.s][self.nd.e].display.max_level /= 1.1
        self.set_disp_pars()
        self.plot_spc()
        # end scale_2d_spectrum_up

    def scale_2d_spectrum_down(self):
        self.nd.nmrdat[self.nd.s][self.nd.e].display.min_level *= 1.1
        self.nd.nmrdat[self.nd.s][self.nd.e].display.max_level *= 1.1
        self.set_disp_pars()
        self.plot_spc()
        # end scale_2d_spectrum_down

    def scale_all_2d_spectra_up(self):
        self.nd.nmrdat[self.nd.s][self.nd.e].display.min_level /= 1.1
        self.nd.nmrdat[self.nd.s][self.nd.e].display.max_level /= 1.1
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            self.nd.nmrdat[self.nd.s][k].display.min_level = self.nd.nmrdat[self.nd.s][self.nd.e].display.min_level
            self.nd.nmrdat[self.nd.s][k].display.max_level = self.nd.nmrdat[self.nd.s][self.nd.e].display.max_level

        self.set_disp_pars()
        self.plot_spc()
        # end scaleAll_2d_spectra_up

    def scale_all_2d_spectra_down(self):
        self.nd.nmrdat[self.nd.s][self.nd.e].display.min_level *= 1.1
        self.nd.nmrdat[self.nd.s][self.nd.e].display.max_level *= 1.1
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            self.nd.nmrdat[self.nd.s][k].display.min_level = self.nd.nmrdat[self.nd.s][self.nd.e].display.min_level
            self.nd.nmrdat[self.nd.s][k].display.max_level = self.nd.nmrdat[self.nd.s][self.nd.e].display.max_level

        self.set_disp_pars()
        self.plot_spc()
        # end scaleAll_2d_spectra_down

    def script_editor(self):
        self.w.nmrSpectrum.setCurrentIndex(10)
        # end script_editor

    def set_datasets_exps(self):
        if self.w.intAllDS.isChecked() == True:
            self.nd.int_all_data_sets = True
        else:
            self.nd.int_all_data_sets = False

        if self.w.intAllExps.isChecked() == True:
            self.nd.int_all_exps = True
        else:
            self.nd.int_all_exps = False

        if self.w.exportFormatCB.currentIndex() == 0:
            self.nd.export_peak_excel = True
        else:
            self.nd.export_peak_excel = False

        start_peak = self.nd.nmrdat[self.nd.s][self.nd.e].start_peak
        end_peak = self.nd.nmrdat[self.nd.s][self.nd.e].end_peak
        peak_label = self.nd.nmrdat[self.nd.s][self.nd.e].peak_label
        self.nd.set_peak(start_peak, end_peak, peak_label)
        # end set_datasets_exps

    def select_add_compress_pre_proc(self):
        self.ginput_compress(2)
        # end select_add_exclude_pre_proc

    def select_add_exclude_pre_proc(self):
        self.ginput_exclude(2)
        # end select_add_exclude_pre_proc

    def select_add_seg_align_pre_proc(self):
        self.ginput_seg_align(2)
        # end select_add_exclude_pre_proc

    def select_all_pre_proc(self):
        n_spc = len(self.nd.pp.class_select)
        self.nd.pp.plot_select = np.arange(n_spc)
        self.fill_pre_processing_numbers()
        self.set_plot_pre_proc()
        self.plot_spc_pre_proc()
        self.w.selectClassTW.setFocus()
        # end select_all_pre_proc

    def select_class_pre_proc(self):
        cls = self.w.selectClassLE.text()
        cls2 = self.nd.pp.class_select
        sel = np.array([])
        for k in range(len(cls2)):
            if (cls2[k] == cls):
                sel = np.append(sel, k)

        if (len(sel) == 0):
            sel = np.arange(len(cls2))

        self.nd.pp.plot_select = sel
        self.fill_pre_processing_numbers()
        self.set_plot_pre_proc()
        self.plot_spc_pre_proc()
        self.w.selectClassTW.setFocus()
        # end select_class_pre_proc

    def select_clear_compress_pre_proc(self):
        self.nd.pp.pre_proc_fill = True
        for k in range(len(self.nd.pp.compress_start)):
            self.w.compressBucketsTW.item(k, 0).setText("")
            self.w.compressBucketsTW.setFocus()
            self.w.compressBucketsTW.item(k, 1).setText("")
            self.w.compressBucketsTW.setFocus()

        self.nd.pp.pre_proc_fill = False
        self.nd.pp.compress_start = np.array([])
        self.nd.pp.compress_end = np.array([])
        self.w.compressBucketsTW.setFocus()
        self.fill_pre_processing_numbers()
        self.w.compressBucketsTW.setFocus()
        self.set_plot_pre_proc()
        self.w.compressBucketsTW.setFocus()
        self.plot_spc_pre_proc()
        self.set_compress_pre_proc()
        self.w.MplWidget.canvas.flush_events()
        self.w.MplWidget.canvas.draw()
        # end select_clear_exclude_pre_proc

    def select_clear_exclude_pre_proc(self):
        self.nd.pp.pre_proc_fill = True
        for k in range(len(self.nd.pp.exclude_start)):
            self.w.excludeRegionTW.item(k, 0).setText("")
            self.w.excludeRegionTW.setFocus()
            self.w.excludeRegionTW.item(k, 1).setText("")
            self.w.excludeRegionTW.setFocus()

        self.nd.pp.pre_proc_fill = False
        self.nd.pp.exclude_start = np.array([])
        self.nd.pp.exclude_end = np.array([])
        self.w.excludeRegionTW.setFocus()
        self.fill_pre_processing_numbers()
        self.w.excludeRegionTW.setFocus()
        self.set_plot_pre_proc()
        self.w.excludeRegionTW.setFocus()
        self.plot_spc_pre_proc()
        self.set_exclude_pre_proc()
        self.w.MplWidget.canvas.flush_events()
        self.w.MplWidget.canvas.draw()
        # end select_clear_exclude_pre_proc

    def select_clear_seg_align_pre_proc(self):
        self.nd.pp.pre_proc_fill = True
        for k in range(len(self.nd.pp.seg_start)):
            self.w.segAlignTW.item(k, 0).setText("")
            self.w.segAlignTW.setFocus()
            self.w.segAlignTW.item(k, 1).setText("")
            self.w.segAlignTW.setFocus()

        self.nd.pp.pre_proc_fill = False
        self.nd.pp.seg_start = np.array([])
        self.nd.pp.seg_end = np.array([])
        self.w.segAlignTW.setFocus()
        self.fill_pre_processing_numbers()
        self.w.segAlignTW.setFocus()
        self.set_plot_pre_proc()
        self.w.segAlignTW.setFocus()
        self.plot_spc_pre_proc()
        self.set_seg_align_pre_proc()
        self.w.MplWidget.canvas.flush_events()
        self.w.MplWidget.canvas.draw()
        # end select_clear_exclude_pre_proc

    def select_even_pre_proc(self):
        n_spc = len(self.nd.pp.class_select)
        self.nd.pp.plot_select = np.arange(n_spc)
        self.nd.pp.plot_select = self.nd.pp.plot_select[1::2]
        self.fill_pre_processing_numbers()
        self.set_plot_pre_proc()
        self.plot_spc_pre_proc()
        self.w.selectClassTW.setFocus()
        # end select_even_pre_proc

    def select_odd_pre_proc(self):
        n_spc = len(self.nd.pp.class_select)
        self.nd.pp.plot_select = np.arange(n_spc)
        self.nd.pp.plot_select = self.nd.pp.plot_select[0::2]
        self.fill_pre_processing_numbers()
        self.set_plot_pre_proc()
        self.plot_spc_pre_proc()
        self.w.selectClassTW.setFocus()
        # end select_odd_pre_proc

    def select_plot_all(self):
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            self.nd.nmrdat[self.nd.s][k].display.display_spc = True

        self.w.nmrSpectrum.setCurrentIndex(0)
        self.change_data_set_exp()
        self.plot_spc()
        return "select_plot_all"
        # end select_plot_all

    def select_plot_clear(self):
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            self.nd.nmrdat[self.nd.s][k].display.display_spc = False

        # self.plot_spc()
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.change_data_set_exp()
        self.plot_spc()
        return "select_plot_clear"
        # end select_plot_clear

    def select_plot_list(self, plot_select, auto_plot_spc=True):
        plot_select = np.array(plot_select)
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            self.nd.nmrdat[self.nd.s][k].display.display_spc = False

        plot_select -= 1
        for k in range(len(plot_select)):
            if ((plot_select[k] > -1) and (plot_select[k] < len(self.nd.nmrdat[self.nd.s]))):
                self.nd.nmrdat[self.nd.s][plot_select[k]].display.display_spc = True

        # self.plot_spc()
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.change_data_set_exp()
        if auto_plot_spc:
            self.plot_spc()

        return "select_plot_list"
        # end select_plot_list

    def set_acq_pars(self):
        s = self.nd.s
        e = self.nd.e
        a = self.nd.nmrdat[s][e].acq
        acq_str = "originalDataset      " + self.nd.nmrdat[s][e].orig_data_set + "\n"
        acq_str += "___________________________________________________________________________________________________\n"
        acq_str += "\n"
        acq_str += "metaInfo             "
        for k in range(len(a.title)):
            acq_str += a.title[k] + " "

        acq_str += "\n                    "
        acq_str += " Origin\t" + a.origin + "\n                    "
        acq_str += " Owner\t" + a.owner + "\n"
        acq_str += "___________________________________________________________________________________________________\n"
        acq_str += "\n"
        acq_str += "probe                          " + a.probe + "\n"
        pp = a.pul_prog_name
        pp = pp[1:]
        pp = pp[:len(pp) - 1]
        acq_str += "pulseProgram                   " + pp + "\n\n"
        acq_str += "sw                   [ppm]    " + "% 9.2f" % a.sw[0] + "        |    % 9.2f" % a.sw[
            1] + "        |    % 9.2f\n" % a.sw[2]
        acq_str += "sw_h                 [Hz]     " + "% 9.2f" % a.sw_h[0] + "        |    % 9.2f" % a.sw_h[
            1] + "        |    % 9.2f\n" % a.sw_h[2]
        acq_str += "bf1/2/3              [MHz]    " + "% 9.2f" % a.bf1 + "        |    % 9.2f" % a.bf2 + "        |    % 9.2f\n" % a.bf3
        acq_str += "sfo1/2/3             [MHz]    " + "% 9.2f" % a.sfo1 + "        |    % 9.2f" % a.sfo2 + "        |    % 9.2f\n" % a.sfo3
        acq_str += "o1/2/3               [Hz]     " + "% 9.2f" % a.o1 + "        |    % 9.2f" % a.o2 + "        |    % 9.2f\n" % a.o3
        acq_str += "nPoints                       " + "% 6d" % a.n_data_points[0] + "           |    % 6d" % \
                   a.n_data_points[1] + "           |    % 6d\n" % a.n_data_points[2]
        acq_str += "transients                    " + "% 6d\n" % a.transients
        acq_str += "steadyStateScans              " + "% 6d\n\n" % a.steady_state_scans
        acq_str += "groupDelay           [us]     " + "% 9.2f\n" % a.group_delay
        acq_str += "decim                         " + "% 6d\n" % a.decim
        acq_str += "dspfvs                        " + "% 6d\n" % a.dspfvs
        acq_str += "temperature          [K]      " + "% 9.2f\n" % a.temperature
        self.w.acqPars.setText(acq_str)
        # end set_acq_pars

    def set_avoid_neg_values(self):
        if (self.nd.pp.pre_proc_fill == False):
            if (self.w.avoidNegValues.isChecked() == True):
                self.nd.pp.avoid_negative_values = True
            else:
                self.nd.pp.avoid_negative_values = False

        # end set_avoid_neg_values

    def set_bucket_ppm_pre_proc(self):
        try:
            bucket_ppm = float(self.w.bucketPpmLE.text())
        except:
            bucket_ppm = self.nd.pp.bucket_ppm

        ppm_per_point = abs(self.nd.nmrdat[self.nd.s][0].ppm1[0] - self.nd.nmrdat[self.nd.s][0].ppm1[1])
        bucket_points = round(bucket_ppm / ppm_per_point)
        bucket_ppm = np.round(1e4 * bucket_points * ppm_per_point) / 1e4
        self.w.bucketPpmLE.setText(str(bucket_ppm))
        self.w.bucketDataPointsLE.setText(str(int(bucket_points)))
        self.nd.pp.bucket_points = bucket_points
        self.nd.pp.bucket_ppm = bucket_ppm
        # end set_bucket_ppm_pre_proc

    def set_bucket_points_pre_proc(self):
        try:
            bucket_points = float(self.w.bucketDataPointsLE.text())
        except:
            bucket_points = self.nd.pp.bucket_points

        ppm_per_point = abs(self.nd.nmrdat[self.nd.s][0].ppm1[0] - self.nd.nmrdat[self.nd.s][0].ppm1[1])
        bucket_points = round(bucket_points)
        bucket_ppm = np.round(1e4 * bucket_points * ppm_per_point) / 1e4
        self.w.bucketPpmLE.setText(str(bucket_ppm))
        self.w.bucketDataPointsLE.setText(str(int(bucket_points)))
        self.nd.pp.bucket_points = bucket_points
        self.nd.pp.bucket_ppm = bucket_ppm
        # end set_bucket_points_pre_proc

    def set_bucket_spectra(self):
        if (self.nd.pp.pre_proc_fill == False):
            if (self.w.bucketSpectra.isChecked() == True):
                self.nd.pp.flag_bucket_spectra = True
                self.w.preProcessingSelect.setCurrentIndex(4)
            else:
                self.nd.pp.flag_bucket_spectra = False

        # end set_bucket_spectra

    def set_change_pre_proc(self):
        if (self.nd.pp.pre_proc_fill == False):
            cls = np.array([])
            for k in range(len(self.nd.pp.class_select)):
                cls = np.append(cls, self.w.selectClassTW.item(k, 1).text())

            self.nd.pp.class_select = cls

        # end set_change_pre_proc

    def set_compress_buckets(self):
        if (self.nd.pp.pre_proc_fill == False):
            if (self.w.compressBuckets.isChecked() == True):
                self.nd.pp.flag_compress_buckets = True
                self.w.preProcessingSelect.setCurrentIndex(5)
            else:
                self.nd.pp.flag_compress_buckets = False

        # end set_compress_buckets

    def set_compress_pre_proc(self):
        if (self.nd.pp.pre_proc_fill == False):
            n_rows = self.w.compressBucketsTW.rowCount()
            co_start = np.array([])
            co_end = np.array([])
            t_start = np.array([])
            t_end = np.array([])
            for k in range(n_rows):
                # t_start = np.array([])
                # t_end   = np.array([])
                try:
                    t_start = np.append(t_start, float(self.w.compressBucketsTW.item(k, 0).text()))
                    # self.w.compressBucketsTW.item(k,0).clearContents()
                except:
                    t_start = np.append(t_start, -10000.0)

                try:
                    t_end = np.append(t_end, float(self.w.compressBucketsTW.item(k, 1).text()))
                    # self.w.compressBucketsTW.item(k,1).clearContents()
                except:
                    t_end = np.append(t_end, -10000.0)

            # self.w.compressBucketsTW.clearContents()
            self.w.compressBucketsTW.setRowCount(0)
            self.w.compressBucketsTW.setRowCount(n_rows)
            self.nd.pp.pre_proc_fill = True
            for k in np.arange(len(t_start) - 1, -1, -1):  # range(len(t_start)):
                comp_number1 = QTableWidgetItem(2 * k)
                comp_number1.setTextAlignment(QtCore.Qt.AlignHCenter)
                self.w.compressBucketsTW.setItem(k, 0, comp_number1)
                comp_number2 = QTableWidgetItem(2 * k + 1)
                comp_number2.setTextAlignment(QtCore.Qt.AlignHCenter)
                self.w.compressBucketsTW.setItem(k, 1, comp_number2)
                if ((t_start[k] > -10000.0) & (t_end[k] > -10000.0)):
                    t_min = min(t_start[k], t_end[k])
                    t_end[k] = max(t_start[k], t_end[k])
                    t_start[k] = t_min
                    co_start = np.append(co_start, t_start[k])
                    co_end = np.append(co_end, t_end[k])
                    t_start = np.delete(t_start, k)
                    t_end = np.delete(t_end, k)

                if (t_start[k] > -10000.0):
                    self.w.compressBucketsTW.item(k, 0).setText(str(t_start[k]))
                    self.w.compressBucketsTW.setFocus()
                else:
                    self.w.compressBucketsTW.item(k, 0).setText("")
                    self.w.compressBucketsTW.setFocus()

                if (t_end[k] > -10000.0):
                    self.w.compressBucketsTW.item(k, 1).setText(str(t_end[k]))
                    self.w.compressBucketsTW.setFocus()
                else:
                    self.w.compressBucketsTW.item(k, 1).setText("")
                    self.w.compressBucketsTW.setFocus()

            self.nd.pp.pre_proc_fill = False
            sort_idx = np.argsort(co_start)
            self.nd.pp.compress_start = co_start[sort_idx]
            self.nd.pp.compress_end = co_end[sort_idx]
            self.plot_spc_pre_proc()

        # end set_compress_pre_proc

    def set_dark_mode(self):
        self.cf.read_config()
        self.cf.mode = 'dark'
        self.cf.save_config()
        # restart program
        # os.execv(sys.executable, ['python'] + sys.argv)
        # end save_config

    def set_disp_pars(self):
        # print("s: {}, e: {}".format(self.nd.s, self.nd.e))
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        self.w.posColR.setText(str(d.pos_col_rgb[0]))
        self.w.posColG.setText(str(d.pos_col_rgb[1]))
        self.w.posColB.setText(str(d.pos_col_rgb[2]))
        self.w.negColR.setText(str(d.neg_col_rgb[0]))
        self.w.negColG.setText(str(d.neg_col_rgb[1]))
        self.w.negColB.setText(str(d.neg_col_rgb[2]))
        self.w.nLevels.setText(str(d.n_levels))
        self.w.minLevel.setText(str(d.min_level))
        self.w.maxLevel.setText(str(d.max_level))
        self.w.spcOffset.setText(str(d.spc_offset))
        self.w.spcScale.setText(str(d.spc_scale))
        self.w.xLabel.setText(d.x_label)
        self.w.yLabel.setText(d.y_label)
        self.w.spcLabel.setText(d.spc_label)
        self.w.posCol.setCurrentIndex(d.colours.get(d.pos_col))
        self.w.negCol.setCurrentIndex(d.colours.get(d.neg_col))
        self.w.axisType1.setCurrentIndex(d.axes.get(d.axis_type1))
        self.w.axisType2.setCurrentIndex(d.axes.get(d.axis_type2))
        self.w.displaySpc.setCurrentIndex(d.false_true.get(d.display_spc))
        self.w.phRefColour.setCurrentIndex(d.colours2.get(d.ph_ref_col))
        self.w.phRefDS.setValue(d.ph_ref_ds)
        self.w.phRefExp.setValue(d.ph_ref_exp)
        # end set_disp_pars

    def set_add_peak(self):
        if (self.nd.peak_fill == False):
            n_rows = self.w.peakWidget.rowCount()
            start_peak = np.array([])
            end_peak = np.array([])
            peak_label = np.array([])
            s_peak = np.array([])
            e_peak = np.array([])
            p_label = np.array([])
            for k in range(n_rows):
                try:
                    s_peak = np.append(s_peak, float(self.w.peakWidget.item(k, 0).text()))
                except:
                    s_peak = np.append(s_peak, -10000.0)

                try:
                    e_peak = np.append(e_peak, float(self.w.peakWidget.item(k, 1).text()))
                except:
                    e_peak = np.append(e_peak, -10000.0)

                try:
                    p_label = np.append(p_label, self.w.peakWidget.item(k, 2).text())
                except:
                    p_label = np.append(p_label, '')

            print("s: {}, e: {}, l: {}".format(s_peak, e_peak, p_label))
            self.w.peakWidget.setRowCount(0)
            self.w.peakWidget.setRowCount(n_rows)
            self.nd.peak_fill = True
            for k in np.arange(len(s_peak) - 1, -1, -1):  # range(len(t_start)):
                p_number1 = QTableWidgetItem(3 * k)
                p_number1.setTextAlignment(QtCore.Qt.AlignHCenter)
                self.w.peakWidget.setItem(k, 0, p_number1)
                p_number2 = QTableWidgetItem(3 * k + 1)
                p_number2.setTextAlignment(QtCore.Qt.AlignHCenter)
                self.w.peakWidget.setItem(k, 1, p_number2)
                p_label1 = QTableWidgetItem(3 * k + 2)
                p_label1.setTextAlignment(QtCore.Qt.AlignHCenter)
                self.w.peakWidget.setItem(k, 3, p_label1)
                print("k: {}, s_peak[k]: {}, s_peak: {}".format(k, s_peak[k], s_peak))
                if ((s_peak[k] > -10000.0) & (e_peak[k] > -10000.0)):
                    p_min = min(s_peak[k], e_peak[k])
                    e_peak[k] = max(s_peak[k], e_peak[k])
                    s_peak[k] = p_min
                    start_peak = np.append(start_peak, s_peak[k])
                    end_peak = np.append(end_peak, e_peak[k])
                    peak_label = np.append(peak_label, p_label[k])
                    s_peak = np.delete(s_peak, k)
                    e_peak = np.delete(e_peak, k)

            sort_idx = np.argsort(start_peak)
            start_peak = start_peak[sort_idx]
            end_peak = end_peak[sort_idx]
            peak_label = peak_label[sort_idx]
            self.nd.set_peak(start_peak, end_peak, peak_label)
            self.set_peak_picking()
            self.nd.peak_fill = False
            self.plot_spc()

        # end set_exclude_pre_proc

    def set_exclude_pre_proc(self):
        if (self.nd.pp.pre_proc_fill == False):
            n_rows = self.w.excludeRegionTW.rowCount()
            ex_start = np.array([])
            ex_end = np.array([])
            t_start = np.array([])
            t_end = np.array([])
            for k in range(n_rows):
                # t_start = np.array([])
                # t_end   = np.array([])
                try:
                    t_start = np.append(t_start, float(self.w.excludeRegionTW.item(k, 0).text()))
                    # self.w.excludeRegionTW.item(k,0).clearContents()
                except:
                    t_start = np.append(t_start, -10000.0)

                try:
                    t_end = np.append(t_end, float(self.w.excludeRegionTW.item(k, 1).text()))
                    # self.w.excludeRegionTW.item(k,1).clearContents()
                except:
                    t_end = np.append(t_end, -10000.0)

            # self.w.excludeRegionTW.clearContents()
            self.w.excludeRegionTW.setRowCount(0)
            self.w.excludeRegionTW.setRowCount(n_rows)
            self.nd.pp.pre_proc_fill = True
            for k in np.arange(len(t_start) - 1, -1, -1):  # range(len(t_start)):
                excl_number1 = QTableWidgetItem(2 * k)
                excl_number1.setTextAlignment(QtCore.Qt.AlignHCenter)
                self.w.excludeRegionTW.setItem(k, 0, excl_number1)
                excl_number2 = QTableWidgetItem(2 * k + 1)
                excl_number2.setTextAlignment(QtCore.Qt.AlignHCenter)
                self.w.excludeRegionTW.setItem(k, 1, excl_number2)
                if ((t_start[k] > -10000.0) & (t_end[k] > -10000.0)):
                    t_min = min(t_start[k], t_end[k])
                    t_end[k] = max(t_start[k], t_end[k])
                    t_start[k] = t_min
                    ex_start = np.append(ex_start, t_start[k])
                    ex_end = np.append(ex_end, t_end[k])
                    t_start = np.delete(t_start, k)
                    t_end = np.delete(t_end, k)

                if (t_start[k] > -10000.0):
                    self.w.excludeRegionTW.item(k, 0).setText(str(t_start[k]))
                    self.w.excludeRegionTW.setFocus()
                else:
                    self.w.excludeRegionTW.item(k, 0).setText("")
                    self.w.excludeRegionTW.setFocus()

                if (t_end[k] > -10000.0):
                    self.w.excludeRegionTW.item(k, 1).setText(str(t_end[k]))
                    self.w.excludeRegionTW.setFocus()
                else:
                    self.w.excludeRegionTW.item(k, 1).setText("")
                    self.w.excludeRegionTW.setFocus()

            self.nd.pp.pre_proc_fill = False
            sort_idx = np.argsort(ex_start)
            self.nd.pp.exclude_start = ex_start[sort_idx]
            self.nd.pp.exclude_end = ex_end[sort_idx]
            self.plot_spc_pre_proc()

        # end set_exclude_pre_proc

    def set_exclude_region(self):
        if (self.nd.pp.pre_proc_fill == False):
            if (self.w.excludeRegion.isChecked() == True):
                self.nd.pp.flag_exclude_region = True
                self.w.preProcessingSelect.setCurrentIndex(1)
            else:
                self.nd.pp.flag_exclude_region = False

        # end set_exclude_region

    def set_export_character(self):
        tt = self.w.exportCharacter.text()
        if (len(tt) > 0):
            self.nd.pp.export_character = tt[0]
            self.w.exportCharacter.setText(tt[0])

        # end set_export_character

    def set_export_delimiter_tab(self):
        self.nd.pp.export_delimiter_tab = self.w.exportDelimiterTab.isChecked()
        # end set_export_delimiter_tab

    def set_export_file_name(self):
        if self.nd.pp.export_method == 0:
            self.nd.pp.export_excel = self.w.exportFileName.text()

        if self.nd.pp.export_method == 1:
            self.nd.pp.export_file_name = self.w.exportFileName.text()

        if self.nd.pp.export_method == 2:
            self.nd.pp.export_metabo_analyst = self.w.exportFileName.text()

        if self.nd.pp.export_method == 3:
            self.nd.pp.export_r_dolphin = self.w.exportFileName.text()

        if self.nd.pp.export_method == 4:
            self.nd.pp.export_batman = self.w.exportFileName.text()

        if self.nd.pp.export_method == 5:
            self.nd.pp.export_bruker = self.w.exportFileName.text()

        # end set_export_file_name

    def set_export_method(self):
        if self.nd.pp.export_method == 0:
            self.w.delimiterLabel.setHidden(True)
            self.w.exportDelimiterTab.setHidden(True)
            self.w.exportDelimiterCharacter.setHidden(True)
            self.w.exportCharacter.setHidden(True)
            self.w.samplesInRowsLabel.setHidden(False)
            self.w.samplesInComboBox.setHidden(False)
            self.w.exportPath.setText(self.nd.pp.export_excel_path)
            self.w.exportFileName.setText(self.nd.pp.export_excel)

        if self.nd.pp.export_method == 1:
            self.w.delimiterLabel.setHidden(False)
            self.w.exportDelimiterTab.setHidden(False)
            self.w.exportDelimiterCharacter.setHidden(False)
            self.w.exportCharacter.setHidden(False)
            self.w.samplesInRowsLabel.setHidden(False)
            self.w.samplesInComboBox.setHidden(False)
            self.w.exportPath.setText(self.nd.pp.export_path_name)
            self.w.exportFileName.setText(self.nd.pp.export_file_name)

        if self.nd.pp.export_method == 2:
            self.w.delimiterLabel.setHidden(True)
            self.w.exportDelimiterTab.setHidden(True)
            self.w.exportDelimiterCharacter.setHidden(True)
            self.w.exportCharacter.setHidden(True)
            self.w.samplesInRowsLabel.setHidden(True)
            self.w.samplesInComboBox.setHidden(True)
            self.w.exportPath.setText(self.nd.pp.export_metabo_analyst_path)
            self.w.exportFileName.setText(self.nd.pp.export_metabo_analyst)

        if self.nd.pp.export_method == 3:
            self.w.delimiterLabel.setHidden(True)
            self.w.exportDelimiterTab.setHidden(True)
            self.w.exportDelimiterCharacter.setHidden(True)
            self.w.exportCharacter.setHidden(True)
            self.w.samplesInRowsLabel.setHidden(True)
            self.w.samplesInComboBox.setHidden(True)
            self.w.exportPath.setText(self.nd.pp.export_r_dolphin_path)
            self.w.exportFileName.setText(self.nd.pp.export_r_dolphin)

        if self.nd.pp.export_method == 4:
            self.w.delimiterLabel.setHidden(True)
            self.w.exportDelimiterTab.setHidden(True)
            self.w.exportDelimiterCharacter.setHidden(True)
            self.w.exportCharacter.setHidden(True)
            self.w.samplesInRowsLabel.setHidden(True)
            self.w.samplesInComboBox.setHidden(True)
            self.w.exportPath.setText(self.nd.pp.export_batman_path)
            self.w.exportFileName.setText(self.nd.pp.export_batman)

        if self.nd.pp.export_method == 5:
            self.w.delimiterLabel.setHidden(True)
            self.w.exportDelimiterTab.setHidden(True)
            self.w.exportDelimiterCharacter.setHidden(True)
            self.w.exportCharacter.setHidden(True)
            self.w.samplesInRowsLabel.setHidden(True)
            self.w.samplesInComboBox.setHidden(True)
            self.w.exportPath.setText(self.nd.pp.export_bruker_path)
            self.w.exportFileName.setText(self.nd.pp.export_bruker)

        # end set_export_method

    def set_export_method_options(self):
        self.nd.pp.export_method = self.w.exportMethod.currentIndex()
        self.set_export_method()
        # end set_export_method_options

    def set_export_path(self):
        if self.nd.pp.export_method == 0:
            self.nd.pp.export_excel_path = self.w.exportPath.text()

        if self.nd.pp.export_method == 1:
            self.nd.pp.export_path_name = self.w.exportPath.text()

        if self.nd.pp.export_method == 2:
            self.nd.pp.export_metabo_analyst_path = self.w.exportPath.text()

        if self.nd.pp.export_method == 3:
            self.nd.pp.export_r_dolphin_path = self.w.exportPath.text()

        if self.nd.pp.export_method == 4:
            self.nd.pp.export_batman_path = self.w.exportPath.text()

        if self.nd.pp.export_method == 5:
            self.nd.pp.export_bruker_path = self.w.exportPath.text()

        # end set_export_path

    def set_export_data_set(self):
        if (self.nd.pp.pre_proc_fill == False):
            if (self.w.exportDataSet.isChecked() == True):
                self.nd.pp.flag_export_data_set = True
                self.w.preProcessingSelect.setCurrentIndex(8)
            else:
                self.nd.pp.flag_export_data_set = False

        # end set_export_data_set

    def set_export_table(self):
        p_name = QFileDialog.getExistingDirectory()
        # p_name = p_name[0]
        if (len(p_name) > 0):
            if self.nd.pp.export_method == 0:
                self.w.exportPath.setText(p_name)
                self.nd.pp.export_excel_path = p_name

            if self.nd.pp.export_method == 1:
                self.w.exportPath.setText(p_name)
                self.nd.pp.export_path_name = p_name

            if self.nd.pp.export_method == 2:
                self.w.exportPath.setText(p_name)
                self.nd.pp.export_metabo_analyst_path = p_name

            if self.nd.pp.export_method == 3:
                self.w.exportPath.setText(p_name)
                self.nd.pp.export_r_dolphin_path = p_name

            if self.nd.pp.export_method == 4:
                self.w.exportPath.setText(p_name)
                self.nd.pp.export_batman_path = p_name

            if self.nd.pp.export_method == 5:
                self.w.exportPath.setText(p_name)
                self.nd.pp.export_bruker_path = p_name

        # end set_export_table

    def set_font_size(self):
        font_size = self.w.fontSize.value()
        f = self.w.acqPars.font()
        f.setPointSize(font_size)
        self.w.acqPars.setFont(f)
        self.w.titleFile.setFont(f)
        cursor = self.w.script.textCursor()
        self.w.script.selectAll()
        self.w.script.setFontPointSize(font_size)
        self.w.script.setTextCursor(cursor)
        self.w.script.setCurrentFont(f)
        # self.w.script.setFont(f)
        cursor = self.w.console.textCursor()
        self.w.console.selectAll()
        self.w.console.setFontPointSize(font_size)
        self.w.console.setTextCursor(cursor)
        self.w.console.setCurrentFont(f)
        # self.w.console.setFont(f)
        self.w.pulseProgram.setFont(f)
        self.w.cmdLine.setFont(f)
        self.w.setStyleSheet("font-size: " + str(font_size) + "pt")
        # end set_font_size

    def set_help(self):
        url = []
        idx = self.w.helpComboBox.currentIndex()
        if self.cf.mode == 'dark':
            f_name = os.path.join(os.path.dirname(__file__), "nmr", "web", "introductionDark", "index.html")
        else:
            f_name = os.path.join(os.path.dirname(__file__), "nmr", "web", "introduction", "index.html")

        url.append("file:///" + f_name.replace('\\', '/'))
        url.append("https://www.hmdb.ca")
        url.append("https://www.smpdb.ca")
        url.append("https://bmrb.io/metabolomics/")
        url.append("https://www.genome.jp/kegg/pathway.html#metabolism")
        url.append("https://nmrshiftdb.nmr.uni-koeln.de")
        url.append("https://sdbs.db.aist.go.jp/sdbs/cgi-bin/cre_index.cgi")
        url.append("http://dmar.riken.jp/spincouple/")
        self.w.helpView.setUrl(url[idx])
        # end set_help

    def set_invert(self):
        s = self.nd.s
        e = self.nd.e
        self.nd.nmrdat[s][e].proc.invert_matrix[0] = self.w.invertMatrix_1.isChecked()
        self.nd.nmrdat[s][e].proc.invert_matrix[1] = self.w.invertMatrix_2.isChecked()
        # end set_invert

    def set_j_res(self):
        if (self.nd.nmrdat[self.nd.s][self.nd.e].acq.fn_mode == 1):
            self.nd.nmrdat[self.nd.s][self.nd.e].display.y_label = '1H'
            self.nd.nmrdat[self.nd.s][self.nd.e].display.axis_type2 = 'Hz'
            self.nd.nmrdat[self.nd.s][self.nd.e].proc.window_type = np.array([5, 3, 0])
            self.nd.nmrdat[self.nd.s][self.nd.e].proc.lb[0] = 0.5

        # end set_j_res

    def set_light_mode(self):
        self.cf.read_config()
        self.cf.mode = 'light'
        self.cf.save_config()
        # end save_config

    def set_noise_filtering(self):
        if (self.nd.pp.pre_proc_fill == False):
            if (self.w.noiseFiltering.isChecked() == True):
                self.nd.pp.flag_noise_filtering = True
                self.w.preProcessingSelect.setCurrentIndex(3)
            else:
                self.nd.pp.flag_noise_filtering = False

        # end set_noise_filtering

    def set_noise_reg_pre_proc(self):
        try:
            th = float(self.w.noiseThresholdLE.text())
        except:
            th = self.nd.pp.noise_threshold

        try:
            ns = float(self.w.noiseRegionStartLE.text())
        except:
            ns = self.nd.pp.noise_start

        try:
            ne = float(self.w.noiseRegionEndLE.text())
        except:
            ne = self.nd.pp.noise_end

        try:
            lw = float(self.w.thLineWidthLE.text())
        except:
            lw = self.nd.pp.th_line_width

        tm = min(ns, ne)
        ne = max(ns, ne)
        ns = tm
        self.nd.pp.noise_threshold = th
        self.nd.pp.noise_start = ns
        self.nd.pp.noise_end = ne
        self.nd.pp.th_line_width = lw
        self.w.noiseThresholdLE.setText(str(th))
        self.w.noiseRegionStartLE.setText(str(ns))
        self.w.noiseRegionEndLE.setText(str(ne))
        self.w.thLineWidthLE.setText(str(lw))
        self.plot_spc_pre_proc()
        # end set_noise_reg_pre_proc

    def set_ph_ref_exp(self, phRefExp, phRefDS=1):
        self.w.phRefDS.setValue(phRefDS)
        self.w.phRefExp.setValue(phRefExp)
        # end set_ph_ref_exp

    def set_plot_pre_proc(self):
        if (self.nd.pp.pre_proc_fill == False):
            sel = np.array([])
            sel = self.w.selectClassTW.selectedIndexes()
            sel2 = np.array([])
            for k in range(len(sel)):
                if (sel[k].column() == 0):
                    sel2 = np.append(sel2, sel[k].row())

            self.nd.pp.plot_select = sel2
            self.plot_spc_pre_proc()

        # end set_plot_pre_proc

    def set_pqn_tsa_scaling(self):
        if self.w.pqnButton.isChecked() is True:
            self.nd.pp.scale_pqn = True
        else:
            self.nd.pp.scale_pqn = False

        self.w.preserveOverallScale.setDisabled(self.nd.pp.scale_pqn)

        # end set_pqn_tsa_scaling

    def set_pre_processing(self):
        if (self.w.preprocessing.isChecked() == True):
            self.w.peakPicking.setChecked(False)
            # self.w.peakPickingTab.setHidden(True)
            self.w.preProcPeak.setCurrentIndex(0)
            if len(self.nd.nmrdat[self.nd.s]) != len(self.nd.pp.class_select):
                self.nd.pre_proc_init()

            self.show_pre_processing()
            self.fill_pre_processing_numbers()
            self.nd.noise_filtering_init()
        else:
            self.hide_pre_processing()
            self.w.preProcPeak.setHidden(True)

        # end set_pre_processing

    def set_peak_picking(self):
        if (self.w.peakPicking.isChecked() == True):
            self.w.preprocessing.setChecked(False)
            # self.w.preProcessingTab.setHidden(True)
            self.w.preProcPeak.setCurrentIndex(1)
            if self.nd.int_all_data_sets == True:
                self.w.intAllDS.setChecked(True)
            else:
                self.w.intAllDS.setChecked(False)

            if self.nd.int_all_exps == True:
                self.w.intAllExps.setChecked(True)
            else:
                self.w.intAllExps.setChecked(False)

            if self.nd.export_peak_excel == True:
                self.w.exportFormatCB.setCurrentIndex(0)
            else:
                self.w.exportFormatCB.setCurrentIndex(1)

            self.show_peak_picking()
            self.fill_peak_numbers()
        else:
            self.exited_peak_picking = True
            self.hide_peak_picking()
            self.w.preProcPeak.setHidden(True)
            self.w.peakPicking.setChecked(False)

        self.plot_spc()
        # end set_pre_processing

    def set_standard_colours(self):
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            if self.cf.mode == 'dark':
                self.nd.nmrdat[self.nd.s][k].display.pos_col_rgb = self.std_pos_col2
                self.nd.nmrdat[self.nd.s][k].display.neg_col_rgb = self.std_neg_col2
            else:
                self.nd.nmrdat[self.nd.s][k].display.pos_col_rgb = self.std_pos_col1
                self.nd.nmrdat[self.nd.s][k].display.neg_col_rgb = self.std_neg_col1

        self.plot_spc()

    def set_hsqc_analysis(self):
        if (self.w.hsqcAnalysis.isChecked() == True):
            self.w.multipletAnalysis.setVisible(True)
            self.w.isotopomerAnalysis.setVisible(True)
            self.w.nmrSpectrum.setTabEnabled(1, True)
            self.w.nmrSpectrum.setTabEnabled(2, True)
            self.w.nmrSpectrum.setStyleSheet(
                "QTabBar::tab::disabled {width: 0; height: 0; margin: 0; padding: 0; border: none;} ")
            self.w.nmrSpectrum.setCurrentIndex(1)
            self.activate_command_line()
            self.activate_command_line()
        else:
            self.w.multipletAnalysis.setChecked(False)
            self.w.isotopomerAnalysis.setChecked(False)
            self.w.multipletAnalysis.setVisible(False)
            self.w.isotopomerAnalysis.setVisible(False)
            self.w.nmrSpectrum.setTabEnabled(1, False)
            self.w.nmrSpectrum.setTabEnabled(2, False)
            self.w.nmrSpectrum.setStyleSheet(
                "QTabBar::tab::disabled {width: 0; height: 0; margin: 0; padding: 0; border: none;} ")
            self.w.nmrSpectrum.setCurrentIndex(0)

        # end set_hsqc_analysis

    def set_multiplet_analysis(self):
        if (self.w.multipletAnalysis.isChecked() == True):
            self.w.nmrSpectrum.setTabEnabled(3, True)
            self.w.nmrSpectrum.setStyleSheet(
                "QTabBar::tab::disabled {width: 0; height: 0; margin: 0; padding: 0; border: none;} ")
            self.w.nmrSpectrum.setCurrentIndex(3)
        else:
            self.w.nmrSpectrum.setTabEnabled(3, False)
            self.w.nmrSpectrum.setStyleSheet(
                "QTabBar::tab::disabled {width: 0; height: 0; margin: 0; padding: 0; border: none;} ")
            self.w.nmrSpectrum.setCurrentIndex(1)

        # end set_multiplet_analysis

    def set_isotopomer_analysis(self):
        if (self.w.isotopomerAnalysis.isChecked() == True):
            self.w.nmrSpectrum.setTabEnabled(4, True)
            self.w.nmrSpectrum.setStyleSheet(
                "QTabBar::tab::disabled {width: 0; height: 0; margin: 0; padding: 0; border: none;} ")
            self.w.nmrSpectrum.setCurrentIndex(4)
        else:
            self.w.nmrSpectrum.setTabEnabled(4, False)
            self.w.nmrSpectrum.setStyleSheet(
                "QTabBar::tab::disabled {width: 0; height: 0; margin: 0; padding: 0; border: none;} ")
            self.w.nmrSpectrum.setCurrentIndex(1)

        # end set_isotopomer_analysis

    def set_pre_processing_options(self):
        cur_idx = self.w.preProcessingSelect.currentIndex()
        self.w.preProcessingWidget.setCurrentIndex(cur_idx)
        if self.nd.nmrdat[self.nd.s][0].acq.manufacturer == 'Bruker':
            if self.w.exportMethod.count() == 5:
                self.w.exportMethod.addItem('Bruker Dataset')

        else:
            if self.w.exportMethod.count() == 6:
                self.w.exportMethod.removeItem(5)

        self.plot_spc_pre_proc()
        # end set_pre_processingOption

    def set_preserve_overall_scale(self):
        self.nd.pp.preserve_overall_scale = self.w.preserveOverallScale.isChecked()
        # end set_preserve_overall_scale

    def set_proc_pars(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        a = self.nd.nmrdat[self.nd.s][self.nd.e].apc
        self.w.zeroFilling.setText(str(p.n_points[0]))
        self.w.zeroFilling_2.setText(str(p.n_points[1]))
        self.w.lb.setText(str(p.lb[0]))
        self.w.gb.setText(str(p.gb[0]))
        self.w.ssb.setText(str(p.ssb[0]))
        self.w.lb_2.setText(str(p.lb[1]))
        self.w.gb_2.setText(str(p.gb[1]))
        self.w.ssb_2.setText(str(p.ssb[1]))
        self.w.ph0.setText(str(p.ph0[0]))
        self.w.ph1.setText(str(p.ph1[0]))
        self.w.ph0_2.setText(str(p.ph0[1]))
        self.w.ph1_2.setText(str(p.ph1[1]))
        self.w.polyOrder.setText(str(p.poly_order))
        self.w.extrapolationSize.setText(str(p.conv_extrapolation_size[0]))
        self.w.windowSize.setText(str(p.conv_window_size[0]))
        self.w.fidOffsetCorrection.setText(str(p.fid_offset_correction))
        self.w.windowFunction.setCurrentIndex(p.window_type[0])
        self.w.windowFunction_2.setCurrentIndex(p.window_type[1])
        self.w.phaseCorrection.setCurrentIndex(p.ph_corr[0])
        self.w.phaseCorrection_2.setCurrentIndex(p.ph_corr[1])
        self.w.waterSuppression.setCurrentIndex(p.water_suppression)
        self.w.stripTransformStart.setText(str(p.strip_start))
        self.w.stripTransformEnd.setText(str(p.strip_end))
        self.w.winType.setCurrentIndex(p.conv_window_type[0])
        self.w.gibbs.setCurrentIndex(p.gibbs_p.get(p.gibbs[0]))
        self.w.gibbs_2.setCurrentIndex(p.gibbs_p.get(p.gibbs[1]))
        self.w.rSpc_p0.setText(str(a.r_spc[0]))
        self.w.rSpc_p1.setText(str(a.r_spc[1]))
        self.w.rSpc_p2.setText(str(a.r_spc[2]))
        self.w.rSpc_p3.setText(str(a.r_spc[3]))
        self.w.rSpc_p4.setText(str(a.r_spc[4]))
        self.w.rSpc_p5.setText(str(a.r_spc[5]))
        self.w.rSpc_p6.setText(str(a.r_spc[6]))
        self.w.iSpc_p0.setText(str(a.i_spc[0]))
        self.w.iSpc_p1.setText(str(a.i_spc[1]))
        self.w.iSpc_p2.setText(str(a.i_spc[2]))
        self.w.iSpc_p3.setText(str(a.i_spc[3]))
        self.w.iSpc_p4.setText(str(a.i_spc[4]))
        self.w.iSpc_p5.setText(str(a.i_spc[5]))
        self.w.iSpc_p6.setText(str(a.i_spc[6]))
        self.w.baselineOrder.setCurrentIndex(a.n_order)
        self.w.baselineCorrection.setCurrentIndex(a.correct_baseline)
        # end set_proc_pars

    def set_pulse_program(self):
        self.w.pulseProgram.setText(self.nd.nmrdat[self.nd.s][self.nd.e].pulse_program)
        # end set_pulse_program

    # def setrDolphinExport(self):
    #    self.nd.pp.rDolphinExport = self.w.rDolphinExport.isChecked()
    #
    def set_samples_in_combo_box(self):
        self.nd.pp.export_samples_in_rows_cols = self.w.samplesInComboBox.currentIndex()
        # end set_samples_in_combo_box

    def set_scale_spectra(self):
        if (self.nd.pp.pre_proc_fill == False):
            if (self.w.scaleSpectra.isChecked() == True):
                self.nd.pp.flag_scale_spectra = True
                self.w.preProcessingSelect.setCurrentIndex(6)
            else:
                self.nd.pp.flag_scale_spectra = False

        # end set_scale_spectra

    def set_seg_align_pre_proc(self):
        if (self.nd.pp.pre_proc_fill == False):
            n_rows = self.w.segAlignTW.rowCount()
            segStart = np.array([])
            segEnd = np.array([])
            t_start = np.array([])
            t_end = np.array([])
            for k in range(n_rows):
                # t_start = np.array([])
                # t_end   = np.array([])
                try:
                    t_start = np.append(t_start, float(self.w.segAlignTW.item(k, 0).text()))
                    # self.w.segAlignTW.item(k,0).clearContents()
                except:
                    t_start = np.append(t_start, -10000.0)

                try:
                    t_end = np.append(t_end, float(self.w.segAlignTW.item(k, 1).text()))
                    # self.w.segAlignTW.item(k,1).clearContents()
                except:
                    t_end = np.append(t_end, -10000.0)

            # self.w.segAlignTW.clearContents()
            self.w.segAlignTW.setRowCount(0)
            self.w.segAlignTW.setRowCount(n_rows)
            self.nd.pp.pre_proc_fill = True
            for k in np.arange(len(t_start) - 1, -1, -1):  # range(len(t_start)):
                seg_number1 = QTableWidgetItem(2 * k)
                seg_number1.setTextAlignment(QtCore.Qt.AlignHCenter)
                self.w.segAlignTW.setItem(k, 0, seg_number1)
                seg_number2 = QTableWidgetItem(2 * k + 1)
                seg_number2.setTextAlignment(QtCore.Qt.AlignHCenter)
                self.w.segAlignTW.setItem(k, 1, seg_number2)
                if ((t_start[k] > -10000.0) & (t_end[k] > -10000.0)):
                    t_min = min(t_start[k], t_end[k])
                    t_end[k] = max(t_start[k], t_end[k])
                    t_start[k] = t_min
                    segStart = np.append(segStart, t_start[k])
                    segEnd = np.append(segEnd, t_end[k])
                    t_start = np.delete(t_start, k)
                    t_end = np.delete(t_end, k)

                if (t_start[k] > -10000.0):
                    self.w.segAlignTW.item(k, 0).setText(str(t_start[k]))
                    self.w.segAlignTW.setFocus()
                else:
                    self.w.segAlignTW.item(k, 0).setText("")
                    self.w.segAlignTW.setFocus()

                if (t_end[k] > -10000.0):
                    self.w.segAlignTW.item(k, 1).setText(str(t_end[k]))
                    self.w.segAlignTW.setFocus()
                else:
                    self.w.segAlignTW.item(k, 1).setText("")
                    self.w.segAlignTW.setFocus()

            self.nd.pp.pre_proc_fill = False
            sort_idx = np.argsort(segStart)
            self.nd.pp.seg_start = segStart[sort_idx]
            self.nd.pp.seg_end = segEnd[sort_idx]
            self.plot_spc_pre_proc()

        # end set_seg_align_pre_proc

    def set_segmental_alignment(self):
        if (self.nd.pp.pre_proc_fill == False):
            if (self.w.segmentalAlignment.isChecked() == True):
                self.nd.pp.flag_segmental_alignment = True
                self.w.preProcessingSelect.setCurrentIndex(2)
            else:
                self.nd.pp.flag_segmental_alignment = False

        # end set_segmental_alignment

    def set_select_class(self):
        for k in range(len(self.nd.pp.class_select)):
            self.w.selectClassTW.item(k, 1).setText(self.nd.pp.class_select[k])

        # end set_select_class

    def set_sym_j(self):
        cur_idx = self.w.symJ.currentIndex()
        if (cur_idx == 0):
            self.nd.nmrdat[self.nd.s][self.nd.e].proc.symj = True
            self.nd.nmrdat[self.nd.s][self.nd.e].proc.tilt = True
            self.w.tilt.setCurrentIndex(0)
        else:
            self.nd.nmrdat[self.nd.s][self.nd.e].proc.symj = False

        # end set_tilt

    def set_tilt(self):
        cur_idx = self.w.tilt.currentIndex()
        if (cur_idx == 0):
            self.nd.nmrdat[self.nd.s][self.nd.e].proc.tilt = True
        else:
            self.nd.nmrdat[self.nd.s][self.nd.e].proc.tilt = False
            self.nd.nmrdat[self.nd.s][self.nd.e].proc.symj = False
            self.w.symJ.setCurrentIndex(1)

        # end set_tilt

    def set_title_file(self):
        self.w.titleFile.setText(self.nd.nmrdat[self.nd.s][self.nd.e].title)
        # end set_title_file

    def setup_processing_parameters(self):
        self.w.nmrSpectrum.setCurrentIndex(5)
        # end setup_processing_parameters

    def set_variance_stabilisation(self):
        if (self.nd.pp.pre_proc_fill == False):
            if (self.w.varianceStabilisation.isChecked() == True):
                self.nd.pp.flag_variance_stabilisation = True
                self.w.preProcessingSelect.setCurrentIndex(7)
            else:
                self.nd.pp.flag_variance_stabilisation = False

        # end set_variance_stabilisation

    def set_variance_stabilisation_options(self):
        if self.w.autoScaling.isChecked():
            self.nd.pp.auto_scaling = True
            self.nd.pp.pareto_scaling = False
            self.nd.pp.g_log_transform = False
        elif self.w.paretoScaling.isChecked():
            self.nd.pp.auto_scaling = False
            self.nd.pp.pareto_scaling = True
            self.nd.pp.g_log_transform = False
        else:
            self.nd.pp.auto_scaling = False
            self.nd.pp.pareto_scaling = False
            self.nd.pp.g_log_transform = True

        self.w.lambdaText.setEnabled(self.nd.pp.g_log_transform)
        self.w.y0Text.setEnabled(self.nd.pp.g_log_transform)
        self.w.lambdaLE.setEnabled(self.nd.pp.g_log_transform)
        self.w.y0LE.setEnabled(self.nd.pp.g_log_transform)
        # end set_variance_stabilisation_options

    def set_var_lambda(self):
        self.nd.pp.var_lambda = float(self.w.lambdaLE.text())
        # end set_var_lambda

    def set_var_y0(self):
        self.nd.pp.var_y0 = float(self.w.y0LE.text())
        # end set_var_lambda

    def set_pan(self):
        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.set_zoom_release)
        cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
        try:
            self.w.MplWidget.canvas.figure.canvas.toolbar.pan()
        except:
            pass

        self.zoom_was_on = False
        self.pan_was_on = True

    def set_zoom(self):
        try:
            self.w.MplWidget.canvas.figure.canvas.toolbar.zoom()
        except:
            pass

        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.set_zoom_release)
        self.zoom_was_on = True
        self.pan_was_on = False

    def set_zoom_off(self):
        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.set_zoom_release)
        cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)

    def set_zoom_release(self, event):
        if (event.button > 1):
            # Right MB click will unzoom the plot
            try:
                self.w.MplWidget.canvas.figure.canvas.toolbar.home()
            except:
                pass

            pyautogui.click(clicks=1)
            #self.w.MplWidget.setFocus()


    def show(self):
        self.w.show()
        # end show

    def show_acquisition_parameters(self):
        self.w.nmrSpectrum.setCurrentIndex(7)
        # end show_acquisition_parameters

    def show_auto_baseline(self):
        self.w.statusBar().clearMessage()
        self.w.statusBar().showMessage("Automatic baseline correction in progress...")
        self.show_acquisition_parameters()
        self.show_nmr_spectrum()
        # end show_auto_baseline

    def show_auto_phase(self):
        self.w.statusBar().clearMessage()
        self.w.statusBar().showMessage("Automatic phase correction in progress...")
        self.show_acquisition_parameters()
        self.show_nmr_spectrum()
        # end show_auto_phase

    def show_console(self):
        self.w.nmrSpectrum.setCurrentIndex(11)
        # end show_console

    def show_display_parameters(self):
        self.w.nmrSpectrum.setCurrentIndex(6)
        # end show_display_parameters

    def show_help(self):
        self.w.nmrSpectrum.setCurrentIndex(12)
        # end show_help

    def show_main_window(self):
        if (self.w.isFullScreen() == True):
            self.w.showNormal()
        else:
            self.w.showFullScreen()

        # end show_main_window

    def show_nmr_spectrum(self):
        self.w.nmrSpectrum.setCurrentIndex(0)
        # if (self.w.preprocessing.isChecked() == False):
        #    self.plot_spc()
        # end show_nmr_spectrum

    def show_ph_corr(self):
        self.w.statusBar().clearMessage()
        self.w.statusBar().showMessage(
            "Left Mouse Button (MB) for ph0, Right MB or Left MB + shift for ph1, Middle MB or Left MB + Cmd to set pivot")
        #    #"Left Mouse Button (MB) for ph0, Right MB or Left MB + shift for ph1, Middle MB or Left MB + Cmd to set pivot        |        Press Alt+p to exit    |   Press Alt+z to zoom")
        self.show_acquisition_parameters()
        self.show_nmr_spectrum()
        # end show_ph_corr

    def show_ph_corr2d(self):
        self.w.statusBar().clearMessage()
        self.w.statusBar().showMessage(
            "2D Interactive Phase Correction")
        #    "Press: Alt+k to pick row/col | Alt+e to empty selection | Alt+r to remove last row/col | Alt+1 for horizontal phase correction | Alt+2 for vertical phase correction | Alt+x to eXit")
        self.show_acquisition_parameters()
        self.show_nmr_spectrum()
        # end show_ph_corr2d

    def show_ph_corr2d_1d(self, dim=0):
        self.w.statusBar().clearMessage()
        self.w.statusBar().showMessage(
            "Left Mouse Button (MB) for ph0, Right MB or Left MB + shift for ph1, Middle MB or Left MB + Cmd to set pivot")
        #    "Left Mouse Button (MB) for ph0, Right MB or Left MB + shift for ph1, Middle MB or Left MB + Cmd to set pivot | Press: Alt+Shift+p to apply phCorr | Alt+Shift+x to cancel | Alt+z to zoom")
        self.show_acquisition_parameters()
        self.show_nmr_spectrum()
        # end show_ph_corr2d

    def show_ph_zoom(self):
        self.w.statusBar().clearMessage()
        self.w.statusBar().showMessage(
            "Left Mouse Button (MB) for rectangular zoom, Right MB to unzoom")
        #    "Left Mouse Button (MB) for rectangular zoom, Right MB to unzoom        |        Press Alt+z to exit to phase correction")
        self.show_acquisition_parameters()
        self.show_nmr_spectrum()
        # end show_ph_zoom

    def show_pre_processing(self):
        self.w.preProcPeak.setHidden(False)
        self.w.preProcessingTab.setHidden(False)
        self.w.peakPickingTab.setHidden(True)
        self.w.preProcPeak.setTabEnabled(0, True)
        self.w.preProcPeak.setTabEnabled(1, False)
        self.w.preProcessingGroupBox.setHidden(False)
        self.w.preProcessingSelect.setHidden(False)
        self.w.preProcessingWidget.setHidden(False)
        self.w.runPreProcessingButton.setHidden(False)
        self.w.resetPreProcessingButton.setHidden(False)
        self.w.writeScriptButton.setHidden(False)
        self.set_export_method()
        # self.set_select_class()
        self.plot_spc_pre_proc()
        # end show_pre_processing

    def show_pulse_program(self):
        self.w.nmrSpectrum.setCurrentIndex(9)
        # end show_pulse_program

    def show_title_file_information(self):
        self.w.nmrSpectrum.setCurrentIndex(8)
        # end show_title_file_information

    def show_version(self):
        self.w.statusBar().clearMessage()
        self.w.statusBar().showMessage("MetaboLabPy " + self.__version__)
        self.show_acquisition_parameters()
        self.show_nmr_spectrum()
        # end show_version

    def splash(self):
        p_name = os.path.join(os.path.dirname(__file__), "png")
        cf = nmrConfig.NmrConfig()
        cf.read_config()
        if cf.mode == 'dark':
            splash_pix = QPixmap(os.path.join(p_name, "metabolabpy_dark.png"))
        else:
            splash_pix = QPixmap(os.path.join(p_name, "metabolabpy.png"))

        splash = QSplashScreen(splash_pix)
        splash.setMask(splash_pix.mask())
        # adding progress bar
        splash.show()
        QCoreApplication.processEvents()
        max_time = 5
        max_range = 100
        time_inc = max_range
        for i in range(max_range):
            # Simulate something that takes time
            time.sleep(max_time / float(max_range))

        splash.close()
        # end splash

    # def startNotebook(self):
    #    try:
    #        self.p.terminate()
    #        sleep(2)
    #    except:
    #        pass
    #
    #    if self.cf.mode == 'dark':
    #        jupyterthemes.install_theme('chesterish')
    #        #subprocess.run(["jt", "-tchesterish"])
    #    else:
    #        jupyterthemes.install_theme('grade3')
    #        #subprocess.run(["jt", "-r"])
    #
    #    jobs = []
    #    self.p = multiprocess.Process(target=notebookapp.main, args=(['/Users/ludwigc/jupyter', '--no-browser', '--ip=127.0.0.1', '--port=9997'],))
    #    jobs.append(self.p)
    #    self.p.start()
    #    sleep(2)
    #    if self.cf.mode == 'dark':
    #        self.w.helpView.setUrl("http://127.0.0.1:9997/notebooks/test2d_dark.ipynb")
    #    else:
    #        self.w.helpView.setUrl("http://127.0.0.1:9997/notebooks/test2d_light.ipynb")
    #
    #   self.w.nmrSpectrum.setCurrentIndex(12)
    #    # end startNotebook

    def start_stop_ph_corr(self):
        s = self.nd.s
        e = self.nd.e
        # self.set_zoom_off()

        if (self.nd.nmrdat[s][e].dim == 1):
            try:
                self.zoom_was_on = True
                self.w.MplWidget.canvas.figure.canvas.toolbar.zoom()
            except:
                pass

            self.set_zoom_off()
            if (self.ph_corr_active == False):
                self.ph_corr.spc = self.nd.nmrdat[s][e].spc
                self.ph_corr.spc_max = max(max(abs(self.ph_corr.spc)))
                self.ph_corr.piv_points = self.nd.nmrdat[s][e].ppm2points(self.ph_corr.pivot, 0)
                cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.on_ph_corr_click)
                cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.on_ph_corr_release)
                self.ph_corr_active = True
                self.show_ph_corr()
                # self.w.MplWidget.canvas.figure.canvas.toolbar.setEnabled(False)
                self.w.exitPhCorr1d.setVisible(True)
                self.w.zoomPhCorr1d.setVisible(True)
                self.w.exitZoomPhCorr1d.setVisible(False)
                self.update_gui()
                self.ph_corr_plot_spc()
            else:
                cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.on_ph_corr_click)
                cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.on_ph_corr_release)
                cid = self.w.MplWidget.canvas.mpl_disconnect(cid)
                cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
                self.ph_corr_active = False
                # self.w.MplWidget.canvas.figure.canvas.toolbar.setEnabled(True)
                self.show_version()
                self.w.exitPhCorr1d.setVisible(False)
                self.w.zoomPhCorr1d.setVisible(False)
                self.w.exitZoomPhCorr1d.setVisible(False)
                self.update_gui()
                self.plot_spc()
                self.set_zoom_off()
                # self.set_zoom()
                # print("setted Zoom!")
                cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.set_zoom_release)

        else:  # dim == 2
            # try:
            #    self.zoom_was_on = True
            #    self.w.MplWidget.canvas.figure.canvas.toolbar.zoom()
            # except:
            #    pass
            #
            # self.set_zoom_off()
            if (self.ph_corr_active == False):
                self.w.pickRowColPhCorr2d.setVisible(True)
                self.w.emptyRowColPhCorr2d.setVisible(True)
                self.w.removeRowColPhCorr2d.setVisible(True)
                self.w.horzPhCorr2d.setVisible(True)
                self.w.vertPhCorr2d.setVisible(True)
                self.w.exitPhCorr2d.setVisible(True)
                self.ph_corr_active = True
                self.show_ph_corr2d()
            else:
                self.empty_col_row()
                self.w.pickRowColPhCorr2d.setVisible(False)
                self.w.emptyRowColPhCorr2d.setVisible(False)
                self.w.removeRowColPhCorr2d.setVisible(False)
                self.w.horzPhCorr2d.setVisible(False)
                self.w.vertPhCorr2d.setVisible(False)
                self.w.exitPhCorr2d.setVisible(False)
                self.ph_corr_active = False
                self.show_version()
                # self.set_zoom_off()
                # self.set_zoom()
                cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.set_zoom_release)

        self.show_acquisition_parameters()
        self.show_nmr_spectrum()
        # end start_stop_ph_corr

    @contextlib.contextmanager
    def stdoutIO(self, stdout=None):
        old = sys.stdout
        if stdout is None:
            stdout = StringIO()
        sys.stdout = stdout
        yield stdout
        sys.stdout = old
        # end stdoutIO

    def tutorials(self):
        # url = "http://beregond.bham.ac.uk/~ludwigc/tutorials"
        f_name = os.path.join(os.path.dirname(__file__), "nmr", "web", "tutorials", "index.html")
        url = "file:///" + f_name.replace('\\', '/')
        self.w.helpView.setUrl(url)
        self.w.nmrSpectrum.setCurrentIndex(12)
        # end tutorials

    def update_gui(self):
        s = self.nd.s
        e = self.nd.e
        self.w.setBox.valueChanged.disconnect()
        self.w.expBox.valueChanged.disconnect()
        self.w.expBox.setValue(e + 1)
        self.w.setBox.setValue(s + 1)
        self.w.setBox.valueChanged.connect(lambda: self.change_data_set_exp())
        self.w.expBox.valueChanged.connect(lambda: self.change_data_set_exp())
        self.set_disp_pars()
        self.set_proc_pars()
        self.set_acq_pars()
        self.set_title_file()
        self.set_pulse_program()
        self.w.expBox.setValue(e + 1)
        self.w.invertMatrix_1.setChecked(self.nd.nmrdat[s][e].proc.invert_matrix[0])
        self.w.invertMatrix_2.setChecked(self.nd.nmrdat[s][e].proc.invert_matrix[1])
        if (self.nd.nmrdat[s][e].dim == 1):
            self.w.preprocessing.setVisible(True)
            self.w.peakPicking.setVisible(True)
        else:
            self.w.preprocessing.setVisible(False)
            self.w.peakPicking.setVisible(False)

        if self.nd.nmrdat[s][e].dim > 1:
            if self.nd.nmrdat[self.nd.s][self.nd.e].acq.pul_prog_name.find("hsqc") > 0 or self.nd.nmrdat[self.nd.s][
                self.nd.e].acq.pul_prog_name.find("hmqc") > 0:
                self.w.hsqcAnalysis.setVisible(True)  # develop set true
            else:
                self.w.hsqcAnalysis.setVisible(False)

        else:
            self.w.hsqcAnalysis.setVisible(False)

        self.w.multipletAnalysis.setVisible(False)
        self.w.isotopomerAnalysis.setVisible(False)
        return "updated GUI"
        # end update_gui

    def vertical_auto_scale(self):
        if (self.nd.nmrdat[self.nd.s][self.nd.e].dim == 1):
            lines = self.w.MplWidget.canvas.axes.get_lines()
            bottom, top = np.inf, -np.inf
            for line in lines:
                newBottom, newTop = self.get_bottom_top(line)
                if (newBottom < bottom): bottom = newBottom
                if (newTop > top): top = newTop

            if bottom != np.inf and top != -np.inf:
                self.w.MplWidget.canvas.axes.set_ylim(bottom, top)

            self.w.MplWidget.canvas.draw()

        # end vertical_auto_scale

    def vert_ph_corr_2d(self):
        s = self.nd.s
        e = self.nd.e
        self.ph_corr.n_dims = 2
        self.ph_corr.dim = 1
        n_lines = len(self.ph_corr.spc_col_pts)
        if n_lines > 0:
            npts0 = len(self.nd.nmrdat[s][e].spc)
            npts = len(self.nd.nmrdat[s][e].spc[0])
            self.ph_corr.spc = np.zeros((n_lines, npts0), dtype='complex')
            spc1 = np.copy(self.nd.nmrdat[s][e].spc)
            spc1 = np.ndarray.transpose(spc1)
            for k in range(n_lines):
                spc = np.array([spc1[npts - self.ph_corr.spc_col_pts[k]]])
                spc = self.hilbert(spc)
                self.ph_corr.spc[k] = spc[0]

            self.ph_corr.ppm = self.nd.nmrdat[s][e].ppm2
            if self.ph_corr.pivot_points2d[1] < 0:
                self.ph_corr.pivot_points2d[1] = int(len(self.ph_corr.ppm) / 2)
                self.ph_corr.pivot2d[1] = self.nd.nmrdat[s][e].points2ppm(self.ph_corr.pivot_points2d[1], 1)

        self.show_ph_corr2d_1d(self.ph_corr.dim)
        self.ph_corr.spc_max = np.max(np.max(np.abs(self.ph_corr.spc)))
        try:
            zwo = True
            self.w.MplWidget.canvas.figure.canvas.toolbar.zoom()
        except:
            pass

        self.set_zoom_off()
        self.ph_corr.max_ph0 = 90.0
        self.ph_corr.max_ph1 = 90.0
        cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.on_ph_corr_click_2d)
        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.on_ph_corr_release_2d)
        # self.w.actionApplyPhCorr.triggered.connect(self.apply_2d_ph_corr)
        # self.w.actionCancelPhCorr.triggered.connect(self.cancel2dPhCorr)
        self.ph_corr_plot_spc_2d(False)
        self.w.pickRowColPhCorr2d.setVisible(False)
        self.w.emptyRowColPhCorr2d.setVisible(False)
        self.w.removeRowColPhCorr2d.setVisible(False)
        self.w.horzPhCorr2d.setVisible(False)
        self.w.vertPhCorr2d.setVisible(False)
        self.w.zoomPhCorr2d.setVisible(True)
        self.w.applyPhCorr2d.setVisible(True)
        self.w.cancelPhCorr2d.setVisible(True)
        self.w.exitPhCorr2d.setVisible(False)
        self.w.exitZoomPhCorr2d.setVisible(False)
        self.ph_corr_plot_spc_2d(False)
        self.show_acquisition_parameters()
        self.show_nmr_spectrum()
        # end vert_ph_corr_2d

    def zero_acq_pars(self):
        self.w.acqPars.setText("")
        # end zero_acq_pars

    def zero_console(self):
        self.w.console.setText("")
        # end zero_console

    def zero_disp_pars(self):
        self.w.posColR.setText("")
        self.w.posColG.setText("")
        self.w.posColB.setText("")
        self.w.negColR.setText("")
        self.w.negColG.setText("")
        self.w.negColB.setText("")
        self.w.nLevels.setText("")
        self.w.minLevel.setText("")
        self.w.maxLevel.setText("")
        self.w.spcOffset.setText("")
        self.w.spcScale.setText("")
        self.w.xLabel.setText("")
        self.w.yLabel.setText("")
        self.w.spcLabel.setText("")
        self.w.posCol.setCurrentIndex(0)
        self.w.negCol.setCurrentIndex(0)
        self.w.axisType1.setCurrentIndex(0)
        self.w.axisType2.setCurrentIndex(0)
        self.w.displaySpc.setCurrentIndex(0)
        self.w.phRefDS.setValue(0)
        self.w.phRefExp.setValue(0)
        # end zero_disp_pars

    def zero_proc_pars(self):
        self.w.zeroFilling.setText("")
        self.w.zeroFilling_2.setText("")
        self.w.lb.setText("")
        self.w.gb.setText("")
        self.w.ssb.setText("")
        self.w.lb_2.setText("")
        self.w.gb_2.setText("")
        self.w.ssb_2.setText("")
        self.w.ph0.setText("")
        self.w.ph1.setText("")
        self.w.ph0_2.setText("")
        self.w.ph1_2.setText("")
        self.w.polyOrder.setText("")
        self.w.extrapolationSize.setText("")
        self.w.windowSize.setText("")
        self.w.fidOffsetCorrection.setText("")
        self.w.windowFunction.setCurrentIndex(0)
        self.w.windowFunction_2.setCurrentIndex(0)
        self.w.phaseCorrection.setCurrentIndex(0)
        self.w.phaseCorrection_2.setCurrentIndex(0)
        self.w.waterSuppression.setCurrentIndex(0)
        self.w.winType.setCurrentIndex(0)
        self.w.gibbs.setCurrentIndex(0)
        self.w.gibbs_2.setCurrentIndex(0)
        # end zero_proc_pars

    def zero_pulse_program(self):
        self.w.pulseProgram.setText("")
        # end zero_pulse_program

    def zero_script(self):
        self.w.script.setText("")
        # end zero_console

    def zero_title_file(self):
        self.w.titleFile.setText("")
        # end zero_title_file

    def zoom_ph_corr(self):
        if (self.ph_corr_active == True):
            try:
                self.set_zoom()
            except:
                pass

            if (self.zoom == False):
                # Enable zoom
                self.zoom = True
                self.show_ph_zoom()
                if self.ph_corr.n_dims == 1:
                    self.w.exitPhCorr1d.setVisible(False)
                    self.w.zoomPhCorr1d.setVisible(False)
                    self.w.exitZoomPhCorr1d.setVisible(True)
                else:
                    self.w.zoomPhCorr2d.setVisible(False)
                    self.w.applyPhCorr2d.setVisible(False)
                    self.w.cancelPhCorr2d.setVisible(False)
                    self.w.exitZoomPhCorr2d.setVisible(True)

            else:
                # Disable zoom
                self.zoom = False
                if self.ph_corr.n_dims == 1:
                    self.show_ph_corr()
                    self.w.exitPhCorr1d.setVisible(True)
                    self.w.zoomPhCorr1d.setVisible(True)
                    self.w.exitZoomPhCorr1d.setVisible(False)
                    self.set_zoom_off()
                else:
                    self.show_ph_corr2d_1d()
                    self.w.zoomPhCorr2d.setVisible(True)
                    self.w.applyPhCorr2d.setVisible(True)
                    self.w.cancelPhCorr2d.setVisible(True)
                    self.w.exitZoomPhCorr2d.setVisible(False)
                    self.set_zoom_off()

        self.show_acquisition_parameters()
        self.show_nmr_spectrum()
        # end zoom_ph_corr

    def _download_requested(self, download_item) -> None:
        code_out = io.StringIO()
        code_err = io.StringIO()
        sys.stdout = code_out
        sys.stderr = code_err
        dialog = QFileDialog(self.w)
        path = dialog.getSaveFileName(dialog, "Save File", download_item.path())
        print(path)
        if path[0]:
            download_item.setPath(path[0])
            print(f"downloading file to:( {download_item.path()} )")
            download_item.accept()
            self.download_item = download_item
            self.download_item.finished.connect(self._download_finished)
        else:
            print("Download canceled")

        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        if self.cf.mode == 'dark':
            txt_col = QColor.fromRgbF(1.0, 1.0, 1.0, 1.0)
            err_col = QColor.fromRgbF(1.0, 0.5, 0.5, 1.0)
        else:
            txt_col = QColor.fromRgbF(0.0, 0.0, 0.0, 1.0)
            err_col = QColor.fromRgbF(1.0, 0.0, 0.0, 1.0)

        self.w.console.setTextColor(txt_col)
        self.w.console.append(code_out.getvalue())
        self.w.console.setTextColor(err_col)
        self.w.console.append(code_err.getvalue())
        code_out.close()
        code_err.close()
        self.w.nmrSpectrum.setCurrentIndex(11)

    def _download_finished(self) -> None:
        code_out = io.StringIO()
        code_err = io.StringIO()
        sys.stdout = code_out
        sys.stderr = code_err
        print("Download complete")
        if self.w.autoUnzip.isChecked() == True:
            f_name, fExt = os.path.splitext(self.download_item.downloadFileName())
            if fExt == '.zip':
                print('Extracting .zip-file')
                with zipfile.ZipFile(
                        os.path.join(self.download_item.downloadDirectory(), self.download_item.downloadFileName()),
                        'r') as zip_ref:
                    zip_ref.extractall(os.path.join(self.download_item.downloadDirectory(), f_name))

                print('.zip-file extraction finished')

        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        if self.cf.mode == 'dark':
            txt_col = QColor.fromRgbF(1.0, 1.0, 1.0, 1.0)
            err_col = QColor.fromRgbF(1.0, 0.5, 0.5, 1.0)
        else:
            txt_col = QColor.fromRgbF(0.0, 0.0, 0.0, 1.0)
            err_col = QColor.fromRgbF(1.0, 0.0, 0.0, 1.0)

        self.w.console.setTextColor(txt_col)
        self.w.console.append(code_out.getvalue())
        self.w.console.setTextColor(err_col)
        self.w.console.append(code_err.getvalue())
        code_out.close()
        code_err.close()
        self.w.nmrSpectrum.setCurrentIndex(11)


def main():  # pragma: no cover
    sys.argv.append('None')
    ap = argparse.ArgumentParser()
    ap.add_argument("-s", "--script", required=False, help="optional script argument")
    ap.add_argument("-ns", "--noSplash", required=False, help="turn splash screen off", action="store_true")
    ap.add_argument("-fs", "--FullScreen", required=False, help="open applicatin in full screen mode",
                    action="store_true")
    ap.add_argument("fileName", metavar="fileName", type=str, help="load MetaboLabPy DataSet File")
    dd = ap.parse_known_args()
    # dd = ap.parse_known_intermixed_args()
    if len(dd[1]) > 0:
        sys.argv.pop()

    args = vars(ap.parse_args())
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_ShareOpenGLContexts)
    app = QApplication(['MetaboLabPy'])
    icon = QIcon()
    p_name = os.path.join(os.path.dirname(__file__), "icon")
    icon.addFile(os.path.join(p_name, "icon-16.png"), QtCore.QSize(16, 16))
    icon.addFile(os.path.join(p_name, "icon-24.png"), QtCore.QSize(24, 24))
    icon.addFile(os.path.join(p_name, "icon-32.png"), QtCore.QSize(32, 32))
    icon.addFile(os.path.join(p_name, "icon-48.png"), QtCore.QSize(48, 48))
    icon.addFile(os.path.join(p_name, "icon-256.png"), QtCore.QSize(256, 256))
    app.setWindowIcon(icon)
    app.setApplicationDisplayName("MetaboLabPy")
    w = main_w()
    if args["FullScreen"] == True:
        w.w.showFullScreen()

    if args["noSplash"] == False:
        ##
        # Create and display the splash screen
        p_name = os.path.join(os.path.dirname(__file__), "png")
        cf = nmrConfig.NmrConfig()
        cf.read_config()
        if cf.mode == 'dark':
            splash_pix = QPixmap(os.path.join(p_name, "metabolabpy_dark.png"))
        else:
            splash_pix = QPixmap(os.path.join(p_name, "metabolabpy.png"))

        splash = QSplashScreen(splash_pix)
        splash.setMask(splash_pix.mask())
        # adding progress bar
        splash.show()
        app.processEvents()
        max_time = 2
        max_range = 30
        time_inc = max_range
        for i in range(max_range):
            # Simulate something that takes time
            time.sleep(max_time / float(max_range))

        splash.close()
        ## End of splash screen

    if args["fileName"] != "None":
        try:
            w.load_file(args["fileName"])
        except:
            if (args["script"] != None):
                w.open_script(args["script"])
                w.script_editor()
                w.exec_script()

    else:
        if (args["script"] != None):
            w.open_script(args["script"])
            w.script_editor()
            w.exec_script()

    if cf.mode == 'light':
        qtmodern.styles.light(app)
    else:
        qtmodern.styles.dark(app)

    w.show()
    sys.exit(app.exec_())


if __name__ == "__main__":  # pragma: no cover
    main()
