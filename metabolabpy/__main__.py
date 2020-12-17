#!/usr/bin/env python
import sys   # pragma: no cover
import matplotlib   # pragma: no cover
if "linux" in sys.platform:  # pragma: no cover
    gui_env = ['TkAgg', 'GTKAgg', 'Qt5Agg', 'WXAgg']  # pragma: no cover
elif sys.platform == "darwin":  # pragma: no cover
    try:  # pragma: no cover
        import PySide2  # pragma: no cover
        gui_env = ['Qt5Agg']  # pragma: no cover
    except ImportError:  # pragma: no cover
        gui_env = ['TkAgg', 'GTKAgg', 'Qt5Agg', 'WXAgg']  # pragma: no cover
else:  # pragma: no cover
    pass  # pragma: no cover

if sys.platform != "win32":  # pragma: no cover
    for gui in gui_env:  # pragma: no cover
        try:  # pragma: no cover
            matplotlib.use(gui, warn=False, force=True)  # pragma: no cover
            break  # pragma: no cover
        except:  # pragma: no cover
            continue  # pragma: no cover


#matplotlib.use('Qt5Agg')  # pragma: no cover
try:
    from matplotlib.backends.backend_qt5agg import (FigureCanvasQTAgg as FigureCanvas,
                                                    NavigationToolbar2QT as NavigationToolbar)  # pragma: no cover
except:
    pass

from matplotlib.figure import Figure  # pragma: no cover
import matplotlib.pyplot as pl  # pragma: no cover
import argparse  # pragma: no cover
from PySide2.QtUiTools import QUiLoader  # pragma: no cover
from PySide2.QtCore import QFile  # pragma: no cover
from PySide2.QtWidgets import *  # pragma: no cover
from PySide2 import QtWidgets  # pragma: no cover
from PySide2.QtGui import *  # pragma: no cover
from PySide2 import QtGui  # pragma: no cover
from PySide2 import QtCore  # pragma: no cover
from PySide2.QtWidgets import QFileDialog # pragma: no cover
from PySide2.QtCore import SIGNAL  # pragma: no cover
from time import sleep
from PySide2.QtCore import QUrl, Qt
from PySide2.QtWebEngineCore import QWebEngineUrlSchemeHandler
from PySide2.QtWebEngineWidgets import QWebEngineView, QWebEngineProfile, QWebEnginePage, QWebEngineSettings
try:  # pragma: no cover
    import pyautogui  # pragma: no cover
except:  # pragma: no cover
    pass  # pragma: no cover

import qtmodern.styles
import numpy as np  # pragma: no cover
import io  # pragma: no cover
from metabolabpy.nmr import nmrDataSet  # pragma: no cover
from metabolabpy.GUI import phCorr  # pragma: no cover
import time  # pragma: no cover
import platform  # pragma: no cover
import math  # pragma: no cover
from metabolabpy.nmr import nmrConfig  # pragma: no cover
import os  # pragma: no cover
import traceback  # pragma: no cover
import shutil  # pragma: no cover
import scipy.io  # pragma: no cover
import inspect
from io import StringIO
import contextlib
import zipfile

# import pandas as pd                       # pragma: no cover

# ------------------ MplWidget ------------------
class MplWidget(QWidget):  # pragma: no cover

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

        fig = Figure()
        self.canvas = FigureCanvas(fig)

        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)
        self.toolbar = NavigationToolbar(self.canvas, self)
        vertical_layout.addWidget(self.toolbar)

        self.canvas.axes = self.canvas.figure.add_subplot(111)
        self.setLayout(vertical_layout)
        home = NavigationToolbar.home

        def new_home(self, *args, **kwargs):
            self.canvas.axes.autoscale()
            self.canvas.draw()
            self.canvas.toolbar.update()
            home(self, *args, **kwargs)

        NavigationToolbar.home = new_home
        self.phCorr = phCorr.PhCorr()
        # end __init__


# ------------------ MplWidget ------------------

# ------------------ MplWidget ------------------
class MplWidget2(QWidget):  # pragma: no cover

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

        self.canvas = FigureCanvas(Figure())

        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)
        self.toolbar = NavigationToolbar(self.canvas, self)
        vertical_layout.addWidget(self.toolbar)

        self.canvas.axes = self.canvas.figure.add_subplot(111)
        self.setLayout(vertical_layout)
        home = NavigationToolbar.home

        def new_home(self, *args, **kwargs):
            self.canvas.axes.autoscale()
            self.canvas.draw()
            self.canvas.toolbar.update()
            home(self, *args, **kwargs)

        NavigationToolbar.home = new_home
        self.phCorr = phCorr.PhCorr()
        # end __init__

class QWebEngineView2(QWebEngineView):
    def printCmd(self):
        print("QWebEngineView2")

# ------------------ MplWidget2 ------------------

class main_w(object):  # pragma: no cover
    def __init__(self):
        self.__version__ = '0.6.18'
        self.zoomWasOn = True
        self.panWasOn = False
        self.stdPosCol1 = (0.0, 0.0, 1.0)
        self.stdNegCol1 = (1.0, 0.0, 0.0)
        self.stdPosCol2 = (0.8, 0.8, 1.0)
        self.stdNegCol2 = (1.0, 0.8, 0.8)
        self.nClicks = 1
        self.curClicks = 0
        self.xy = [[]]
        self.xdata = []
        self.ydata = []
        self.nd = nmrDataSet.NmrDataSet()
        self.phCorr = phCorr.PhCorr()
        # load ui; create w
        fName = os.path.join(os.path.dirname(__file__), "ui", "metabolabpy_mainwindow.ui")
        self.file = QFile(fName)
        self.file.open(QFile.ReadOnly)
        self.loader = QUiLoader()
        self.loader.registerCustomWidget(QWebEngineView2)
        self.loader.registerCustomWidget(MplWidget)
        #self.loader.registerCustomWidget(MplWidget2)
        self.w = self.loader.load(self.file)
        self.zoom = False

        self.hidePreProcessing()
        self.w.preprocessing.setVisible(False)
        self.w.hsqcAnalysis.setVisible(False)
        self.w.multipletAnalysis.setVisible(False)
        self.w.isotopomerAnalysis.setVisible(False)
        self.w.nmrSpectrum.setTabEnabled(1, False)
        self.w.nmrSpectrum.setTabEnabled(2, False)
        self.w.nmrSpectrum.setTabEnabled(3, False)
        self.w.nmrSpectrum.setTabEnabled(4, False)
        self.w.nmrSpectrum.setStyleSheet("QTabBar::tab::disabled {width: 0; height: 0; margin: 0; padding: 0; border: none;} ")
        # connections
        # self.w.rDolphinExport.clicked.connect(self.setrDolphinExport)
        self.w.exportPath.textChanged.connect(self.setExportPath)
        self.w.exportFileName.textChanged.connect(self.setExportFileName)
        self.w.exportDelimiterTab.toggled.connect(self.setExportDelimiterTab)
        self.w.exportCharacter.returnPressed.connect(self.setExportCharacter)
        self.w.samplesInComboBox.currentIndexChanged.connect(self.setSamplesInComboBox)
        self.w.runPreProcessingButton.clicked.connect(self.dataPreProcessing)
        self.w.resetPreProcessingButton.clicked.connect(self.resetDataPreProcessing)
        self.w.avoidNegValues.stateChanged.connect(self.setAvoidNegValues)
        self.w.excludeRegion.stateChanged.connect(self.setExcludeRegion)
        self.w.segmentalAlignment.stateChanged.connect(self.setSegmentalAlignment)
        self.w.compressBuckets.stateChanged.connect(self.setCompressBuckets)
        self.w.noiseFiltering.stateChanged.connect(self.setNoiseFiltering)
        self.w.bucketSpectra.stateChanged.connect(self.setBucketSpectra)
        self.w.scaleSpectraRefSpc.valueChanged.connect(self.changeScaleSpectraRefSpc)
        self.w.segAlignRefSpc.valueChanged.connect(self.changeSegAlignRefSpc)
        self.w.scaleSpectra.stateChanged.connect(self.setScaleSpectra)
        self.w.pqnButton.clicked.connect(self.setPqnTsaScaling)
        self.w.tsaButton.clicked.connect(self.setPqnTsaScaling)
        self.w.autoScaling.clicked.connect(self.setVarianceStabilisationOptions)
        self.w.paretoScaling.clicked.connect(self.setVarianceStabilisationOptions)
        self.w.gLogTransform.clicked.connect(self.setVarianceStabilisationOptions)
        self.w.varianceStabilisation.stateChanged.connect(self.setVarianceStabilisation)
        self.w.exportDataSet.stateChanged.connect(self.setExportDataSet)
        self.w.excludeRegionTW.cellChanged.connect(self.setExcludePreProc)
        self.w.segAlignTW.cellChanged.connect(self.setSegAlignPreProc)
        self.w.selectClassTW.itemSelectionChanged.connect(self.setPlotPreProc)
        self.w.selectClassTW.cellChanged.connect(self.setChangePreProc)
        self.w.excludeClearButton.clicked.connect(self.selectClearExcludePreProc)
        self.w.segAlignClearButton.clicked.connect(self.selectClearSegAlignPreProc)
        self.w.compressClearButton.clicked.connect(self.selectClearCompressPreProc)
        self.w.excludeAddButton.clicked.connect(self.selectAddExcludePreProc)
        self.w.segAlignAddButton.clicked.connect(self.selectAddSegAlignPreProc)
        self.w.compressAddButton.clicked.connect(self.selectAddCompressPreProc)
        self.w.selectAllButton.clicked.connect(self.selectAllPreProc)
        self.w.selectEvenButton.clicked.connect(self.selectEvenPreProc)
        self.w.selectOddButton.clicked.connect(self.selectOddPreProc)
        self.w.selectClassButton.clicked.connect(self.selectClassPreProc)
        self.w.selectClassLE.returnPressed.connect(self.selectClassPreProc)
        self.w.cmdLine.returnPressed.connect(self.execCmd)
        self.w.noiseThresholdLE.returnPressed.connect(self.setnoiseRegPreProc)
        self.w.noiseRegionStartLE.returnPressed.connect(self.setnoiseRegPreProc)
        self.w.noiseRegionEndLE.returnPressed.connect(self.setnoiseRegPreProc)
        self.w.thLineWidthLE.returnPressed.connect(self.setnoiseRegPreProc)
        self.w.bucketPpmLE.returnPressed.connect(self.setBucketPPMPreProc)
        self.w.bucketDataPointsLE.returnPressed.connect(self.setBucketPointsPreProc)
        self.w.actionVertical_AutoScale.triggered.connect(self.verticalAutoScale)
        self.w.actionZoom.triggered.connect(self.setZoom)
        self.w.actionPan.triggered.connect(self.setPan)
        self.w.actionShow_Next_Tab.triggered.connect(self.nextTab)
        self.w.actionShow_Previous_Tab.triggered.connect(self.previousTab)
        self.w.actionPlot_spc.triggered.connect(self.plotSpc)
        self.w.actionSave.triggered.connect(self.saveButton)
        self.w.actionLoad.triggered.connect(self.loadButton)
        self.w.actionOpen_NMRPipe.triggered.connect(self.readNMRPipeSpc)
        self.w.actionActivate_Command_Line.triggered.connect(self.activateCommandLine)
        self.w.actionPrevious_command.triggered.connect(self.previousCommand)
        self.w.actionNext_command.triggered.connect(self.nextCommand)
        self.w.actionCorrect_Phase.triggered.connect(self.startStopPhCorr)
        #self.w.actionZoomCorrect_Phase.triggered.connect(self.zoomPhCorr)
        self.w.zoomPhCorr1d.clicked.connect(self.zoomPhCorr)
        self.w.exitZoomPhCorr1d.clicked.connect(self.zoomPhCorr)
        self.w.exitPhCorr1d.clicked.connect(self.startStopPhCorr)
        self.w.actionClear.triggered.connect(self.clear)
        self.w.lambdaLE.textChanged.connect(self.setVarLambda)
        self.w.y0LE.textChanged.connect(self.setVary0)
        self.w.actionRead_NMR_Spectrum.triggered.connect(self.readNMRSpc)
        self.w.preprocessing.stateChanged.connect(self.setPreProcessing)
        self.w.hsqcAnalysis.stateChanged.connect(self.setHsqcAnalysis)
        self.w.multipletAnalysis.stateChanged.connect(self.setMultipletAnalysis)
        self.w.isotopomerAnalysis.stateChanged.connect(self.setIsotopomerAnalysis)
        self.w.preserveOverallScale.stateChanged.connect(self.setPreserveOverallScale)
        self.w.actionReset.triggered.connect(self.resetPlot)
        self.w.actionShow_NMR_Spectrum.triggered.connect(self.showNMRSpectrum)
        self.w.actionSetup_Processing_Parameters.triggered.connect(self.setupProcessingParameters)
        self.w.actionShow_Display_Parameters.triggered.connect(self.showDisplayParameters)
        self.w.actionShow_Acquisition_Parameters.triggered.connect(self.showAcquisitionParameters)
        self.w.actionShow_Title_File_Information.triggered.connect(self.showTitleFileInformation)
        self.w.actionShow_pulseProgram.triggered.connect(self.showPulseProgram)
        self.w.actionFourier_Transform.triggered.connect(self.ft)
        self.w.actionScript_Editor.triggered.connect(self.scriptEditor)
        self.w.actionChange_to_next_Exp.triggered.connect(self.changeToNextExp)
        self.w.actionChange_to_previous_Exp.triggered.connect(self.changeToPreviousExp)
        self.w.actionChange_to_next_DS.triggered.connect(self.changeToNextDS)
        self.w.actionChange_to_previous_DS.triggered.connect(self.changeToPreviousDS)
        self.w.exampleScripts.view().pressed.connect(self.loadExampleScript)
        self.w.actionAutomatic_Phase_Correction.triggered.connect(self.autophase1d)
        self.w.actionAutomatic_Baseline_Correction.triggered.connect(self.autobaseline1d)
        self.w.actionScale_2D_Spectrum_Up.triggered.connect(self.scale2DSpectrumUp)
        self.w.actionScale_2D_Spectrum_Down.triggered.connect(self.scale2DSpectrumDown)
        self.w.actionScale_all_2D_Spectra_Up.triggered.connect(self.scaleAll2DSpectraUp)
        self.w.actionScale_all_2D_Spectra_Down.triggered.connect(self.scaleAll2DSpectraDown)
        self.w.actionSelect_All.triggered.connect(self.selectPlotAll)
        self.w.actionClear_All.triggered.connect(self.selectPlotClear)
        self.w.actionConsole.triggered.connect(self.showConsole)
        self.w.actionHelp.triggered.connect(self.showHelp)
        self.w.actionToggle_FullScreen.triggered.connect(self.showMainWindow)
        self.w.setBox.valueChanged.connect(self.changeDataSetExp)
        self.w.expBox.valueChanged.connect(self.changeDataSetExp)
        self.w.posCol.currentIndexChanged.connect(self.getDispPars1)
        self.w.negCol.currentIndexChanged.connect(self.getDispPars2)
        self.w.posColR.returnPressed.connect(self.getDispPars3)
        self.w.posColG.returnPressed.connect(self.getDispPars3)
        self.w.posColB.returnPressed.connect(self.getDispPars3)
        self.w.negColR.returnPressed.connect(self.getDispPars3)
        self.w.negColG.returnPressed.connect(self.getDispPars3)
        self.w.negColB.returnPressed.connect(self.getDispPars3)
        self.w.nLevels.returnPressed.connect(self.getDispPars4)
        self.w.minLevel.returnPressed.connect(self.getDispPars5)
        self.w.maxLevel.returnPressed.connect(self.getDispPars6)
        self.w.axisType1.currentIndexChanged.connect(self.getDispPars7)
        self.w.axisType2.currentIndexChanged.connect(self.getDispPars8)
        self.w.displaySpc.currentIndexChanged.connect(self.getDispPars9)
        self.w.baselineCorrection.currentIndexChanged.connect(self.checkBaselineCorrection)
        self.w.baselineOrder.currentIndexChanged.connect(self.checkBaselineOrder)
        self.w.spcOffset.returnPressed.connect(self.getDispPars10)
        self.w.spcScale.returnPressed.connect(self.getDispPars11)
        self.w.fontSize.valueChanged.connect(self.setFontSize)
        self.w.xLabel.returnPressed.connect(self.getDispPars12)
        self.w.yLabel.returnPressed.connect(self.getDispPars13)
        self.w.spcLabel.returnPressed.connect(self.getDispPars14)
        self.w.preProcessingSelect.currentIndexChanged.connect(self.setPreProcessingOptions)
        self.w.exportMethod.currentIndexChanged.connect(self.setExportMethodOptions)
        self.w.tilt.currentIndexChanged.connect(self.setTilt)
        self.w.symJ.currentIndexChanged.connect(self.setSymJ)
        self.w.windowFunction.currentIndexChanged.connect(self.getProcPars1)
        self.w.windowFunction_2.currentIndexChanged.connect(self.getProcPars2)
        self.w.phaseCorrection.currentIndexChanged.connect(self.getProcPars3)
        self.w.phaseCorrection_2.currentIndexChanged.connect(self.getProcPars4)
        self.w.waterSuppression.currentIndexChanged.connect(self.getProcPars5)
        self.w.winType.currentIndexChanged.connect(self.getProcPars6)
        self.w.gibbs.currentIndexChanged.connect(self.getProcPars7)
        self.w.gibbs_2.currentIndexChanged.connect(self.getProcPars8)
        self.w.zeroFilling.returnPressed.connect(self.getProcPars9)
        self.w.zeroFilling_2.returnPressed.connect(self.getProcPars10)
        self.w.lb.returnPressed.connect(self.getProcPars11)
        self.w.gb.returnPressed.connect(self.getProcPars12)
        self.w.ssb.returnPressed.connect(self.getProcPars13)
        self.w.lb_2.returnPressed.connect(self.getProcPars14)
        self.w.gb_2.returnPressed.connect(self.getProcPars15)
        self.w.ssb_2.returnPressed.connect(self.getProcPars16)
        self.w.ph0.returnPressed.connect(self.getProcPars17)
        self.w.ph1.returnPressed.connect(self.getProcPars18)
        self.w.ph0_2.returnPressed.connect(self.getProcPars19)
        self.w.ph1_2.returnPressed.connect(self.getProcPars20)
        self.w.polyOrder.returnPressed.connect(self.getProcPars21)
        self.w.extrapolationSize.returnPressed.connect(self.getProcPars22)
        self.w.windowSize.returnPressed.connect(self.getProcPars23)
        self.w.fidOffsetCorrection.returnPressed.connect(self.getProcPars24)
        self.w.stripTransformStart.returnPressed.connect(self.getProcPars25)
        self.w.stripTransformEnd.returnPressed.connect(self.getProcPars26)
        self.w.phRefDS.valueChanged.connect(self.changeDataSetExpPhRef)
        self.w.phRefExp.valueChanged.connect(self.changeDataSetExpPhRef)
        self.w.phRefColour.currentIndexChanged.connect(self.getDispPars15)
        self.w.fourierTransformButton.clicked.connect(self.ft)
        self.w.executeScript.clicked.connect(self.execScript)
        self.w.openScript.clicked.connect(self.openScript)
        self.w.saveScript.clicked.connect(self.saveScript)
        self.w.actionOpen_Script.triggered.connect(self.openScript)
        self.w.actionSave_Script.triggered.connect(self.saveScript)
        self.w.actionExecute_Script.triggered.connect(self.execScript)
        #self.w.helpComboBox.currentIndexChanged.connect(self.setHelp)
        self.w.helpComboBox.activated.connect(self.setHelp)
        # Quit Button
        self.w.quitButton.clicked.connect(self.quit_app)
        self.w.saveButton.clicked.connect(self.saveButton)
        self.w.loadButton.clicked.connect(self.loadButton)
        self.w.exportPathSelectButton.clicked.connect(self.setExportTable)
        self.w.actionQuit.triggered.connect(self.quit_app)
        self.w.dispPlotButton.clicked.connect(self.plotSpcDisp)
        self.showVersion()
        self.keepZoom = False
        self.keepXZoom = False
        self.phCorrActive = False
        self.setFontSize()
        self.cf = nmrConfig.NmrConfig()
        self.cf.readConfig()
        self.w.autoPlot.setChecked(self.cf.autoPlot)
        self.w.keepZoom.setChecked(self.cf.keepZoom)
        self.w.fontSize.setValue(self.cf.fontSize)
        self.stdPosCol1 = (self.cf.posCol10,self.cf.posCol11,self.cf.posCol12)
        self.stdNegCol1 = (self.cf.negCol10,self.cf.negCol11,self.cf.negCol12)
        self.stdPosCol2 = (self.cf.posCol20,self.cf.posCol21,self.cf.posCol22)
        self.stdNegCol2 = (self.cf.negCol20,self.cf.negCol21,self.cf.negCol22)
        self.w.actionSave_as_Default.triggered.connect(self.saveConfig)
        self.w.actionLoad_Default.triggered.connect(self.loadConfig)
        self.w.actionReset_Config.triggered.connect(self.resetConfig)
        self.w.rSpc_p0.returnPressed.connect(self.get_rSpc_p0)
        self.w.rSpc_p1.returnPressed.connect(self.get_rSpc_p1)
        self.w.rSpc_p2.returnPressed.connect(self.get_rSpc_p2)
        self.w.rSpc_p3.returnPressed.connect(self.get_rSpc_p3)
        self.w.rSpc_p4.returnPressed.connect(self.get_rSpc_p4)
        self.w.rSpc_p5.returnPressed.connect(self.get_rSpc_p5)
        self.w.rSpc_p6.returnPressed.connect(self.get_rSpc_p6)
        self.w.iSpc_p0.returnPressed.connect(self.get_iSpc_p0)
        self.w.iSpc_p1.returnPressed.connect(self.get_iSpc_p1)
        self.w.iSpc_p2.returnPressed.connect(self.get_iSpc_p2)
        self.w.iSpc_p3.returnPressed.connect(self.get_iSpc_p3)
        self.w.iSpc_p4.returnPressed.connect(self.get_iSpc_p4)
        self.w.iSpc_p5.returnPressed.connect(self.get_iSpc_p5)
        self.w.iSpc_p6.returnPressed.connect(self.get_iSpc_p6)
        self.setFontSize()
        self.w.MplWidget.toolbar.setVisible(False)
        self.w.MplWidget2.toolbar.setVisible(False)
        self.w.MplWidget.setFocus()
        self.setZoom()
        self.w.pickRowColPhCorr2d.clicked.connect(self.pickColRow)
        self.w.emptyRowColPhCorr2d.clicked.connect(self.emptyColRow)
        self.w.removeRowColPhCorr2d.clicked.connect(self.removeLastColRow)
        self.w.horzPhCorr2d.clicked.connect(self.horzPhCorr2d)
        self.w.vertPhCorr2d.clicked.connect(self.vertPhCorr2d)
        self.w.exitPhCorr2d.clicked.connect(self.startStopPhCorr)
        self.w.applyPhCorr2d.clicked.connect(self.apply2dPhCorr)
        self.w.cancelPhCorr2d.clicked.connect(self.cancel2dPhCorr)
        self.w.zoomPhCorr2d.clicked.connect(self.zoomPhCorr)
        self.w.exitZoomPhCorr2d.clicked.connect(self.zoomPhCorr)
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
        self.w.actionSet_light_mode_requires_restart.triggered.connect(self.setLightMode)
        self.w.actionSet_dark_mode_requires_restart.triggered.connect(self.setDarkMode)
        self.w.MplWidget.canvas.draw()
        self.w.setStyleSheet("font-size: " + str(self.cf.fontSize) + "pt")
        if self.cf.mode == 'dark':
            self.loadDarkMode()
        else:
            self.loadLightMode()

        self.w.helpView.page().profile().downloadRequested.connect(self._download_requested)
        # end __init__

    def activateCommandLine(self):
        if (self.w.cmdLine.hasFocus() == True):
            self.w.cmdLine.clearFocus()
        else:
            self.w.cmdLine.setFocus()

        # end activateCommandLine

    def apply2dPhCorr(self):
        s = self.nd.s
        e = self.nd.e
        if self.nd.nmrdat[s][e].proc.phaseInversion:
            self.phCorr.ph0_2d[self.phCorr.dim] *= -1
            self.phCorr.ph1_2d[self.phCorr.dim] *= -1

        cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.onPhCorrClick2d)
        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.onPhCorrRelease2d)
        cid = self.w.MplWidget.canvas.mpl_disconnect(cid)
        cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
        #self.w.actionApplyPhCorr.triggered.disconnect()
        #self.w.actionCancelPhCorr.triggered.disconnect()
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
        ph0 = ((self.phCorr.ph0_2d[self.phCorr.dim] + 180.0) % 360.0) - 180.0
        ph1 = self.phCorr.ph1_2d[self.phCorr.dim]
        if self.nd.nmrdat[s][e].proc.phaseInversion is False:
            self.nd.nmrdat[s][e].phase2a(ph0, ph1, self.phCorr.dim)
        else:
            self.nd.nmrdat[s][e].phase2a(-ph0, -ph1, self.phCorr.dim)

        ph0 = ((ph0 + self.nd.nmrdat[s][e].proc.ph0[self.phCorr.dim] + 180.0) % 360.0) - 180.0
        ph1 = ph1 + self.nd.nmrdat[s][e].proc.ph1[self.phCorr.dim]

        self.nd.nmrdat[s][e].proc.ph0[self.phCorr.dim] = ph0
        self.nd.nmrdat[s][e].proc.ph1[self.phCorr.dim] = ph1
        self.phCorr.ph0_2d[self.phCorr.dim] = 0
        self.phCorr.ph1_2d[self.phCorr.dim] = 0
        self.phCorr.spc = np.array([[]], dtype='complex')
        self.phCorr.spc2 = np.array([[]], dtype='complex')
        zoomStatus = self.w.keepZoom.isChecked()
        self.w.keepZoom.setChecked(False)
        self.plotSpc()
        self.w.keepZoom.setChecked(zoomStatus)
        self.plot2dColRow()
        if (self.zoomWasOn == True):
            self.setZoomOff()
            self.setZoom()

        if (self.panWasOn == True):
            self.setPan()

        self.showPhCorr2d()
        self.setProcPars()
        self.showAcquisitionParameters()
        self.showNMRSpectrum()
        # end apply2dPhCorr

    def autobaseline1d(self):
        codeOut = io.StringIO()
        codeErr = io.StringIO()
        sys.stdout = codeOut
        sys.stderr = codeErr
        self.showAutoBaseline()
        self.nd.ft()
        self.nd.autoref()
        self.nd.autobaseline1d()
        self.w.baselineCorrection.setCurrentIndex(1)
        self.nd.ft()
        self.nd.baseline1d()
        # self.w.baselineCorrection.setCurrentIndex(1)
        self.setProcPars()
        self.showVersion()
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.changeDataSetExp()
        self.plotSpc()
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        self.w.console.setTextColor('Black')
        self.w.console.append(codeOut.getvalue())
        self.w.console.setTextColor('Red')
        self.w.console.append(codeErr.getvalue())
        codeOut.close()
        codeErr.close()
        # end autobaseline1d

    def autobaseline1dAll(self):
        codeOut = io.StringIO()
        codeErr = io.StringIO()
        sys.stdout = codeOut
        sys.stderr = codeErr
        self.showAutoBaseline()
        self.nd.ft()
        self.nd.autoref()
        self.nd.autobaseline1dAll()
        self.w.baselineCorrection.setCurrentIndex(1)
        self.nd.ft()
        self.nd.baseline1d()
        # self.w.baselineCorrection.setCurrentIndex(1)
        self.setProcPars()
        self.showVersion()
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.changeDataSetExp()
        self.plotSpc()
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        self.w.console.setTextColor('Black')
        self.w.console.append(codeOut.getvalue())
        self.w.console.setTextColor('Red')
        self.w.console.append(codeErr.getvalue())
        codeOut.close()
        codeErr.close()
        # end autobaseline1dAll

    def autophase1d(self):
        codeOut = io.StringIO()
        codeErr = io.StringIO()
        sys.stdout = codeOut
        sys.stderr = codeErr
        self.showAutoPhase()
        self.nd.ft()
        self.nd.autoref()
        self.nd.autophase1d()
        self.w.baselineCorrection.setCurrentIndex(1)
        self.nd.ft()
        self.nd.baseline1d()
        self.plotSpc()
        self.setProcPars()
        self.showVersion()
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.changeDataSetExp()
        self.plotSpc()
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        self.w.console.setTextColor('Black')
        self.w.console.append(codeOut.getvalue())
        self.w.console.setTextColor('Red')
        self.w.console.append(codeErr.getvalue())
        codeOut.close()
        codeErr.close()
        # end autophase1d

    def autophase1dAll(self):
        codeOut = io.StringIO()
        codeErr = io.StringIO()
        sys.stdout = codeOut
        sys.stderr = codeErr
        self.showAutoPhase()
        self.nd.ft()
        self.nd.autoref()
        self.nd.autophase1dAll()
        self.w.baselineCorrection.setCurrentIndex(1)
        self.nd.ft()
        self.nd.baseline1d()
        self.plotSpc()
        self.setProcPars()
        self.showVersion()
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.changeDataSetExp()
        self.plotSpc()
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        self.w.console.setTextColor('Black')
        self.w.console.append(codeOut.getvalue())
        self.w.console.setTextColor('Red')
        self.w.console.append(codeErr.getvalue())
        codeOut.close()
        codeErr.close()
        # end autophase1dAll

    def autoref(self, tmsp=True):
        self.nd.autoref(tmsp)
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.changeDataSetExp()
        self.plotSpc()
        return "Autoref"
        # end autoref

    def baseline1d(self):
        self.nd.baseline1d()
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.changeDataSetExp()
        self.plotSpc()
        # end baseline1d

    def cancel2dPhCorr(self):
        cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.onPhCorrClick2d)
        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.onPhCorrRelease2d)
        cid = self.w.MplWidget.canvas.mpl_disconnect(cid)
        cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
        #self.w.actionApplyPhCorr.triggered.disconnect()
        #self.w.actionCancelPhCorr.triggered.disconnect()
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
        self.phCorr.ph0_2d[self.phCorr.dim] = 0
        self.phCorr.ph1_2d[self.phCorr.dim] = 0
        zoomStatus = self.w.keepZoom.isChecked()
        self.w.keepZoom.setChecked(False)
        self.plotSpc()
        self.w.keepZoom.setChecked(zoomStatus)
        self.plot2dColRow()
        if (self.zoomWasOn == True):
            self.setZoomOff()
            self.setZoom()

        if (self.panWasOn == True):
            self.setPan()

        self.showPhCorr2d()
        self.showAcquisitionParameters()
        self.showNMRSpectrum()
        # end cancel2dPhCorr

    def changeDataSetExp(self):
        cidx = self.w.nmrSpectrum.currentIndex()
        if (len(self.nd.nmrdat) > 0):
            if (len(self.nd.nmrdat[self.nd.s]) > 0):
                self.keepZoom = self.w.keepZoom.isChecked()
                oldSet = self.nd.s
                oldExp = self.nd.e
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
                keepZoom = self.w.keepZoom.isChecked()
                if not ((oldSet == self.nd.s) and (oldExp == self.nd.e)):
                    if (self.nd.nmrdat[oldSet][oldExp].dim != self.nd.nmrdat[self.nd.s][self.nd.e].dim):
                        self.keepXZoom = True
                        self.keepZoom = False
                        self.w.keepZoom.setChecked(False)

                    self.setDispPars()
                    self.setProcPars()
                    self.setAcqPars()
                    self.setTitleFile()
                    self.setPulseProgram()
                    if (self.phCorrActive == False):
                        if (self.w.autoPlot.isChecked()):
                            self.plotSpc()
                        elif (self.w.nmrSpectrum.currentIndex() == 0):
                            self.plotSpc()

                    else:
                        if self.nd.nmrdat[self.nd.s][self.nd.e].dim == 1:
                            self.phCorr.spc = self.nd.nmrdat[self.nd.s][self.nd.e].spc
                            self.phCorrPlotSpc()
                        else:
                            self.plotSpc()
                            self.plot2dColRow()

                    self.keepZoom = keepZoom
                    self.w.keepZoom.setChecked(keepZoom)
                # else:
                #    if (self.phCorrActive == False):
                #        if (self.w.autoPlot.isChecked()):
                #            self.plotSpc()
                #        elif (self.w.nmrSpectrum.currentIndex() == 0):
                #            self.plotSpc()
                #
                #    else:
                #        self.phCorr.spc = self.nd.nmrdat[self.nd.s][self.nd.e].spc
                #        self.phCorrPlotSpc()

                self.keepZoom = False

            else:
                self.w.setBox.valueChanged.disconnect()
                self.w.expBox.valueChanged.disconnect()
                self.w.expBox.setValue(0)
                self.w.setBox.setValue(0)
                self.w.setBox.valueChanged.connect(lambda: self.changeDataSetExp())
                self.w.expBox.valueChanged.connect(lambda: self.changeDataSetExp())

            self.updateGUI()
        else:
            self.w.setBox.valueChanged.disconnect()
            self.w.expBox.valueChanged.disconnect()
            self.w.expBox.setValue(0)
            self.w.setBox.setValue(0)
            self.w.setBox.valueChanged.connect(lambda: self.changeDataSetExp())
            self.w.expBox.valueChanged.connect(lambda: self.changeDataSetExp())

        if (self.w.autoPlot.isChecked() is False):
            self.w.nmrSpectrum.setCurrentIndex(cidx)

        # end changeDataSetExp

    def changeDataSetExpPhRef(self):
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
                        self.nd.nmrdat[k][l].display.phRefDS = self.w.phRefDS.value()
                        self.nd.nmrdat[k][l].display.phRefExp = self.w.phRefExp.value()

        # end changeDataSetExpPhRef

    def changeScaleSpectraRefSpc(self):
        self.nd.pp.scaleSpectraRefSpc = self.w.scaleSpectraRefSpc.value()
        # end changeScaleSpectraRefSpc

    def changeSegAlignRefSpc(self):
        self.nd.pp.segAlignRefSpc = self.w.segAlignRefSpc.value()
        # end changeSegAlignRefSpc

    def changeStandardColours(self, posCol1 = (), negCol1 = (), posCol2 = (), negCol2 = ()):
        if len(posCol1) != 3:
            posCol1 = self.stdPosCol1

        if len(negCol1) != 3:
            negCol1 = self.stdNegCol1

        if len(posCol2) != 3:
            posCol2 = self.stdPosCol2

        if len(negCol2) != 3:
            negCol2 = self.stdNegCol2

        self.stdPosCol1 = posCol1
        self.stdNegCol1 = negCol1
        self.stdPosCol2 = posCol2
        self.stdNegCol2 = negCol2
        self.setStandardColours()
    # end changeStandardColours

    def changeToNextDS(self):
        self.w.setBox.setValue(self.w.setBox.value() + 1)
        # end changeToNextDS

    def changeToNextExp(self):
        self.w.expBox.setValue(self.w.expBox.value() + 1)
        # end changeToNextExp

    def changeToPreviousDS(self):
        self.w.setBox.setValue(self.w.setBox.value() - 1)
        # end changeToPreviousDS

    def changeToPreviousExp(self):
        self.w.expBox.setValue(self.w.expBox.value() - 1)
        # end changeToPreviousExp

    def checkBaselineCorrection(self):
        cbl = self.w.baselineCorrection.currentIndex()
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.correctBaseline = cbl
        if (cbl == 1):
            self.w.baselineOrder.setEnabled(True)
        else:
            self.w.baselineOrder.setEnabled(False)

        self.checkBaselineOrder()

        # end checkBaselineCorrection

    def checkBaselineOrder(self):
        blo = self.w.baselineOrder.currentIndex()
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.nOrder = blo
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

        # end checkBaselineOrder

    def clear(self):
        self.w.MplWidget.canvas.axes.clear()
        self.w.MplWidget.canvas.draw()
        self.zeroDispPars()
        self.zeroProcPars()
        self.zeroAcqPars()
        self.zeroConsole()
        self.zeroTitleFile()
        self.zeroPulseProgram()
        self.nd.nmrdat = [[]]
        self.nd.s = 0
        self.nd.e = -1
        self.w.setBox.valueChanged.disconnect()
        self.w.expBox.valueChanged.disconnect()
        self.w.expBox.setValue(0)
        self.w.setBox.setValue(0)
        self.w.setBox.valueChanged.connect(lambda: self.changeDataSetExp())
        self.w.expBox.valueChanged.connect(lambda: self.changeDataSetExp())
        return "Workspace cleared"
        # end clear

    def dataPreProcessing(self):
        self.resetDataPreProcessing()
        self.nd.dataPreProcessing()
        self.plotSpcPreProc()
        self.verticalAutoScale()
        self.w.MplWidget.canvas.flush_events()
        self.w.MplWidget.canvas.draw()
        # end dataPreProcessing

    def emptyColRow(self):
        while len(self.w.MplWidget.canvas.axes.lines) > 0:
            self.w.MplWidget.canvas.axes.lines[0].remove()

        self.w.MplWidget.canvas.draw()
        self.phCorr.spcRow = []
        self.phCorr.spcCol = []
        self.phCorr.spcRowPts = []
        self.phCorr.spcColPts = []
        self.showAcquisitionParameters()
        self.showNMRSpectrum()
        # end emptyColRow

    def enableBaseline(self):
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            self.nd.nmrdat[self.nd.s][k].apc.correctBaseline = 1

        self.w.baselineOrder.setCurrentIndex(self.nd.nmrdat[self.nd.s][0].apc.nOrder)
        self.w.baselineCorrection.setCurrentIndex(1)
        return "baselineCorrection enabled"
        # end enableBaseline

    def execCmd(self):
        if self.cf.mode == 'dark':
            txtCol = QColor.fromRgbF(1.0, 1.0, 1.0, 1.0)
            errCol = QColor.fromRgbF(1.0, 0.5, 0.5, 1.0)
        else:
            txtCol = QColor.fromRgbF(0.0, 0.0, 0.0, 1.0)
            errCol = QColor.fromRgbF(1.0, 0.0, 0.0, 1.0)

        cmdText = self.w.cmdLine.text()
        if (len(cmdText) > 0):
            self.w.nmrSpectrum.setCurrentIndex(11)
            self.w.cmdLine.setText("")
            self.nd.cmdBuffer = np.append(self.nd.cmdBuffer, cmdText)
            self.nd.cmdIdx = len(self.nd.cmdBuffer)
            codeOut = io.StringIO()
            codeErr = io.StringIO()
            sys.stdout = codeOut
            sys.stderr = codeErr
            print(">>> " + cmdText)
            try:
                output = eval(cmdText)
                print(output)
                self.w.console.setTextColor(txtCol)
                self.w.console.append(codeOut.getvalue())
            except:  # (SyntaxError, NameError, TypeError, ZeroDivisionError, AttributeError, ArithmeticError, BufferError, LookupError):
                cmdText2 = "self." + cmdText
                try:
                    output = eval(cmdText2)
                    print(output)
                    self.w.console.setTextColor(txtCol)
                    self.w.console.append(codeOut.getvalue())
                except:
                    traceback.print_exc()
                    self.w.console.setTextColor(errCol)
                    self.w.console.append(codeOut.getvalue())
                    self.w.console.append(codeErr.getvalue())

            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__
            codeOut.close()
            codeErr.close()
            self.w.console.verticalScrollBar().setValue(self.w.console.verticalScrollBar().maximum())

        # end execCmd

    def execScript(self):
        if self.cf.mode == 'dark':
            txtCol = QColor.fromRgbF(1.0, 1.0, 1.0, 1.0)
            errCol = QColor.fromRgbF(1.0, 0.5, 0.5, 1.0)
            scrCol = QColor.fromRgbF(0.5, 0.5, 1.0, 1.0)
            scrCol2 = QColor.fromRgbF(0.4, 0.4, 1.0, 1.0)
            admCol = QColor.fromRgbF(1.0, 1.0, 0.5, 1.0)
        else:
            txtCol = QColor.fromRgbF(0.0, 0.0, 0.0, 1.0)
            errCol = QColor.fromRgbF(1.0, 0.0, 0.0, 1.0)
            scrCol = QColor.fromRgbF(0.0, 0.0, 1.0, 1.0)
            scrCol2 = QColor.fromRgbF(0.0, 0.0, 0.6, 1.0)
            admCol = QColor.fromRgbF(0.4, 0.4, 0.4, 1.0)

        zoomChecked = self.w.keepZoom.isChecked()
        self.w.keepZoom.setChecked(False)
        codeOut = io.StringIO()
        codeErr = io.StringIO()
        sys.stdout = codeOut
        sys.stderr = codeErr
        code = self.w.script.toPlainText()
        code = code.replace('\\', '\\' * 2)
        try:
            exec(code)

        except:  # (SyntaxError, NameError, TypeError, ZeroDivisionError, AttributeError):
            self.w.nmrSpectrum.setCurrentIndex(11)
            traceback.print_exc()

        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        self.w.console.setTextColor(admCol)
        self.w.console.append('--- ScriptStart -------------------------\n')
        self.w.console.setTextColor(scrCol2)
        self.w.console.append('Executing script...\n')
        self.w.console.setTextColor(scrCol)
        codeSplit = code.split('\n')
        for k in range(len(codeSplit)):
            self.w.console.append(str(k+1) + ': ' + str(codeSplit[k]))

        self.w.console.setTextColor(admCol)
        self.w.console.append('\n--- ScriptOutput ------------------------\n')
        self.w.console.setTextColor(txtCol)
        self.w.console.append(codeOut.getvalue())
        self.w.console.setTextColor(errCol)
        self.w.console.append(codeErr.getvalue())
        self.w.console.setTextColor(admCol)
        self.w.console.append('--- ScriptEnd ---------------------------\n')
        self.w.console.setTextColor(txtCol)
        codeOut.close()
        codeErr.close()
        if (len(self.nd.nmrdat[0]) > 0):
            self.updateGUI()

        self.w.setBox.valueChanged.disconnect()
        self.w.expBox.valueChanged.disconnect()
        self.w.expBox.setValue(self.nd.e + 1)
        self.w.setBox.setValue(self.nd.s + 1)
        self.w.setBox.valueChanged.connect(lambda: self.changeDataSetExp())
        self.w.expBox.valueChanged.connect(lambda: self.changeDataSetExp())
        if (zoomChecked == True):
            self.w.keepZoom.setChecked(True)
        # end execScript

    def fillPreProcessingNumbers(self):
        self.nd.pp.preProcFill = True
        nSpc = len(self.nd.pp.classSelect)
        self.w.selectClassTW.setRowCount(nSpc)
        for k in range(nSpc):
            spcNumber = QTableWidgetItem(str(k))
            spcNumber.setTextAlignment(QtCore.Qt.AlignHCenter)
            self.w.selectClassTW.setItem(k, 0, spcNumber)
            # self.w.selectClassTW.setItemSelected(spcNumber, False)
            classNumber = QTableWidgetItem(self.nd.pp.classSelect[k])
            classNumber.setTextAlignment(QtCore.Qt.AlignHCenter)
            self.w.selectClassTW.setItem(k, 1, classNumber)

        self.w.selectClassTW.selectAll()
        selIt = self.w.selectClassTW.selectedItems()
        for k in np.arange(len(selIt) - 1, -1, -1):
            if (self.w.selectClassTW.selectedItems()[k].column() == 1):
                self.w.selectClassTW.selectedItems()[k].setSelected(False)

        selIt = self.w.selectClassTW.selectedItems()
        for k in np.arange(len(selIt) - 1, -1, -1):
            if (np.isin(self.w.selectClassTW.selectedItems()[k].row(), self.nd.pp.plotSelect)):
                self.w.selectClassTW.selectedItems()[k].setSelected(True)
            else:
                self.w.selectClassTW.selectedItems()[k].setSelected(False)

        for k in range(len(self.nd.pp.excludeStart)):
            exclNumber1 = QTableWidgetItem(str(2 * k))
            exclNumber1.setTextAlignment(QtCore.Qt.AlignHCenter)
            exclNumber2 = QTableWidgetItem(str(2 * k + 1))
            exclNumber2.setTextAlignment(QtCore.Qt.AlignHCenter)
            self.w.excludeRegionTW.setItem(k, 0, exclNumber1)
            self.w.excludeRegionTW.setItem(k, 1, exclNumber2)
            self.w.excludeRegionTW.item(k, 0).setText(str(self.nd.pp.excludeStart[k]))
            self.w.excludeRegionTW.item(k, 1).setText(str(self.nd.pp.excludeEnd[k]))

        for k in range(len(self.nd.pp.segStart)):
            segNumber1 = QTableWidgetItem(str(2 * k))
            segNumber1.setTextAlignment(QtCore.Qt.AlignHCenter)
            segNumber2 = QTableWidgetItem(str(2 * k + 1))
            segNumber2.setTextAlignment(QtCore.Qt.AlignHCenter)
            self.w.segAlignTW.setItem(k, 0, segNumber1)
            self.w.segAlignTW.setItem(k, 1, segNumber2)
            self.w.segAlignTW.item(k, 0).setText(str(self.nd.pp.segStart[k]))
            self.w.segAlignTW.item(k, 1).setText(str(self.nd.pp.segEnd[k]))

        for k in range(len(self.nd.pp.compressStart)):
            compNumber1 = QTableWidgetItem(str(2 * k))
            compNumber1.setTextAlignment(QtCore.Qt.AlignHCenter)
            compNumber2 = QTableWidgetItem(str(2 * k + 1))
            compNumber2.setTextAlignment(QtCore.Qt.AlignHCenter)
            self.w.compressBucketsTW.setItem(k, 0, compNumber1)
            self.w.compressBucketsTW.setItem(k, 1, compNumber2)
            self.w.compressBucketsTW.item(k, 0).setText(str(self.nd.pp.compressStart[k]))
            self.w.compressBucketsTW.item(k, 1).setText(str(self.nd.pp.compressEnd[k]))

        self.w.noiseThresholdLE.setText(str(self.nd.pp.noiseThreshold))
        self.w.noiseRegionStartLE.setText(str(self.nd.pp.noiseStart))
        self.w.noiseRegionEndLE.setText(str(self.nd.pp.noiseEnd))
        self.w.thLineWidthLE.setText(str(self.nd.pp.thLineWidth))
        self.w.bucketPpmLE.setText(str(self.nd.pp.bucketPPM))
        self.setBucketPPMPreProc()
        self.w.excludeRegion.setChecked(self.nd.pp.flagExcludeRegion)
        self.w.segmentalAlignment.setChecked(self.nd.pp.flagSegmentalAlignment)
        self.w.noiseFiltering.setChecked(self.nd.pp.flagNoiseFiltering)
        self.w.bucketSpectra.setChecked(self.nd.pp.flagBucketSpectra)
        self.w.compressBuckets.setChecked(self.nd.pp.flagCompressBuckets)
        self.w.scaleSpectra.setChecked(self.nd.pp.flagScaleSpectra)
        if self.nd.pp.scalePQN is True:
            self.w.pqnButton.setChecked(True)
        else:
            self.w.tsaButton.setChecked(True)

        self.w.varianceStabilisation.setChecked(self.nd.pp.flagVarianceStabilisation)
        self.w.exportDataSet.setChecked(self.nd.pp.flagExportDataSet)
        self.w.exportDelimiterTab.setChecked(self.nd.pp.exportDelimiterTab)
        self.w.exportDelimiterCharacter.setChecked(not self.nd.pp.exportDelimiterTab)
        self.w.exportCharacter.setText(self.nd.pp.exportCharacter)
        self.w.exportMethod.setCurrentIndex(self.nd.pp.exportMethod)
        self.w.samplesInComboBox.setCurrentIndex(self.nd.pp.exportSamplesInRowsCols)
        self.w.segAlignRefSpc.setMaximum(len(self.nd.nmrdat[self.nd.s]))
        self.w.scaleSpectraRefSpc.setMaximum(len(self.nd.nmrdat[self.nd.s]))
        self.w.segAlignRefSpc.setValue(self.nd.pp.segAlignRefSpc)
        self.w.scaleSpectraRefSpc.setValue(self.nd.pp.scaleSpectraRefSpc)
        self.w.preserveOverallScale.setChecked(self.nd.pp.preserveOverallScale)
        self.w.preserveOverallScale.setDisabled(self.nd.pp.scalePQN)
        self.w.autoScaling.setChecked(self.nd.pp.autoScaling)
        self.w.paretoScaling.setChecked(self.nd.pp.paretoScaling)
        self.w.gLogTransform.setChecked(self.nd.pp.gLogTransform)
        self.w.lambdaText.setEnabled(self.nd.pp.gLogTransform)
        self.w.y0Text.setEnabled(self.nd.pp.gLogTransform)
        self.w.lambdaLE.setEnabled(self.nd.pp.gLogTransform)
        self.w.y0LE.setEnabled(self.nd.pp.gLogTransform)
        self.w.lambdaLE.setText(str(self.nd.pp.varLambda))
        self.w.y0LE.setText(str(self.nd.pp.varY0))
        self.nd.pp.preProcFill = False
        # end fillPreProcessingNumbers

    def ft(self):
        self.nd.ft()
        if (self.w.baselineCorrection.currentIndex() > 0):
            self.baseline1d()

        self.w.nmrSpectrum.setCurrentIndex(0)
        self.changeDataSetExp()
        self.plotSpc()
        # end ft

    def ftAll(self):
        self.nd.ftAll()
        if (self.w.baselineCorrection.currentIndex() > 0):
            self.baseline1dAll()

        self.w.nmrSpectrum.setCurrentIndex(0)
        self.changeDataSetExp()
        self.plotSpc()
        # end ftAll

    def getDispPars1(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.posCol = d.colours.get(self.w.posCol.currentIndex())
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        # end getDispPars1

    def getDispPars2(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.negCol = d.colours.get(self.w.negCol.currentIndex())
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        # end getDispPars2

    def getDispPars3(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        posR = float(self.w.posColR.text())
        posG = float(self.w.posColG.text())
        posB = float(self.w.posColB.text())
        negR = float(self.w.negColR.text())
        negG = float(self.w.negColG.text())
        negB = float(self.w.negColB.text())
        d.posColRGB = (posR, posG, posB)
        d.negColRGB = (negR, negG, negB)
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        # end getDispPars3

    def getDispPars4(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.nLevels = round(float(self.w.nLevels.text()))
        # end getDispPars4

        self.nd.nmrdat[self.nd.s][self.nd.e].display = d

    def getDispPars5(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.minLevel = float(self.w.minLevel.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        # end getDispPars5

    def getDispPars6(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.maxLevel = float(self.w.maxLevel.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        # end getDispPars6

    def getDispPars7(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.axisType1 = d.axes.get(self.w.axisType1.currentIndex())
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        self.nd.nmrdat[self.nd.s][self.nd.e].calcPPM()
        # end getDispPars7

    def getDispPars8(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.axisType2 = d.axes.get(self.w.axisType2.currentIndex())
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        self.nd.nmrdat[self.nd.s][self.nd.e].calcPPM()
        # end getDispPars8

    def getDispPars9(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.displaySpc = d.falseTrue.get(self.w.displaySpc.currentIndex())
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        # end getDispPars9

    def getDispPars10(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.spcOffset = float(self.w.spcOffset.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        # end getDispPars10

    def getDispPars11(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.spcScale = float(self.w.spcScale.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        # end getDispPars11

    def getDispPars12(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.xLabel = self.w.xLabel.text()
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        # end getDispPars12

    def getDispPars13(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.yLabel = self.w.yLabel.text()
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        # end getDispPars13

    def getDispPars14(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.spcLabel = self.w.spcLabel.text()
        self.nd.nmrdat[self.nd.s][self.nd.e].display = d
        # end getDispPars14

    def getDispPars15(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        d.phRefCol = d.colours2.get(self.w.phRefColour.currentIndex())
        for k in range(len(self.nd.nmrdat)):
            for l in range(len(self.nd.nmrdat[k])):
                self.nd.nmrdat[k][l].display.phRefCol = d.phRefCol

        # end getDispPars15

    def getProcPars1(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.windowType[0] = self.w.windowFunction.currentIndex()
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars1

    def getProcPars2(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.windowType[1] = self.w.windowFunction_2.currentIndex()
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars2

    def getProcPars3(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.phCorr[0] = self.w.phaseCorrection.currentIndex()
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars3

    def getProcPars4(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.phCorr[1] = self.w.phaseCorrection_2.currentIndex()
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars4

    def getProcPars5(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.waterSuppression = self.w.waterSuppression.currentIndex()
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars5

    def getProcPars6(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.convWindowType[0] = self.w.winType.currentIndex()
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars6

    def getProcPars7(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.gibbs[0] = p.gibbsP.get(self.w.gibbs.currentIndex())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars7

    def getProcPars8(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.gibbs[1] = p.gibbsP.get(self.w.gibbs_2.currentIndex())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars8

    def getProcPars9(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.nPoints[0] = int(self.w.zeroFilling.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars9

    def getProcPars10(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.nPoints[1] = int(self.w.zeroFilling_2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars10

    def getProcPars11(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.lb[0] = float(self.w.lb.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars11

    def getProcPars12(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.gb[0] = float(self.w.gb.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars12

    def getProcPars13(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.ssb[0] = float(self.w.ssb.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars13

    def getProcPars14(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.lb[1] = float(self.w.lb_2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars14

    def getProcPars15(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.gb[1] = float(self.w.gb_2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars15

    def getProcPars16(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.ssb[1] = float(self.w.ssb_2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars16

    def getProcPars17(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.ph0[0] = float(self.w.ph0.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars17

    def getProcPars18(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.ph1[0] = float(self.w.ph1.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars18

    def getProcPars19(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.ph0[1] = float(self.w.ph0_2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars19

    def getProcPars20(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.ph1[1] = float(self.w.ph1_2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars20

    def getProcPars21(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.polyOrder = int(self.w.polyOrder.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars21

    def getProcPars22(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.convExtrapolationSize[0] = int(self.w.extrapolationSize.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars22

    def getProcPars23(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.convWindowSize[0] = int(self.w.windowSize.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars23

    def getProcPars24(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.fidOffsetCorrection = int(self.w.fidOffsetCorrection.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars24

    def getProcPars25(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.stripStart = int(self.w.stripTransformStart.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars25

    def getProcPars26(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.stripEnd = int(self.w.stripTransformEnd.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars25

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

    def get_rSpc_p0(self):
        r = self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc
        r[0] = float(self.w.rSpc_p0.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = r
        # end get_rSpc_p0

    def get_rSpc_p1(self):
        r = self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc
        r[1] = float(self.w.rSpc_p1.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = r
        # end get_rSpc_p1

    def get_rSpc_p2(self):
        r = self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc
        r[2] = float(self.w.rSpc_p2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = r
        # end get_rSpc_p2

    def get_rSpc_p3(self):
        r = self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc
        r[3] = float(self.w.rSpc_p3.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = r
        # end get_rSpc_p3

    def get_rSpc_p4(self):
        r = self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc
        r[4] = float(self.w.rSpc_p4.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = r
        # end get_rSpc_p4

    def get_rSpc_p5(self):
        r = self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc
        r[5] = float(self.w.rSpc_p5.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = r
        # end get_rSpc_p5

    def get_rSpc_p6(self):
        r = self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc
        r[6] = float(self.w.rSpc_p6.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = r
        # end get_rSpc_p6

    def get_iSpc_p0(self):
        i = self.nd.nmrdat[self.nd.s][self.nd.e].apc.iSpc
        i[0] = float(self.w.iSpc_p0.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = i
        # end get_iSpc_p0

    def get_iSpc_p1(self):
        i = self.nd.nmrdat[self.nd.s][self.nd.e].apc.iSpc
        i[1] = float(self.w.iSpc_p1.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = i
        # end get_iSpc_p1

    def get_iSpc_p2(self):
        i = self.nd.nmrdat[self.nd.s][self.nd.e].apc.iSpc
        i[2] = float(self.w.iSpc_p2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = i
        # end get_iSpc_p2

    def get_iSpc_p3(self):
        i = self.nd.nmrdat[self.nd.s][self.nd.e].apc.iSpc
        i[3] = float(self.w.iSpc_p3.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = i
        # end get_iSpc_p3

    def get_iSpc_p4(self):
        i = self.nd.nmrdat[self.nd.s][self.nd.e].apc.iSpc
        i[4] = float(self.w.iSpc_p4.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = i
        # end get_iSpc_p4

    def get_iSpc_p5(self):
        i = self.nd.nmrdat[self.nd.s][self.nd.e].apc.iSpc
        i[5] = float(self.w.iSpc_p5.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = i
        # end get_iSpc_p5

    def get_iSpc_p6(self):
        i = self.nd.nmrdat[self.nd.s][self.nd.e].apc.iSpc
        i[6] = float(self.w.iSpc_p6.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = i
        # end get_iSpc_p6

    def onGinputClick(self, event):
        self.curClicks += 1
        if self.curClicks < self.nClicks:
            self.xdata.append(event.xdata)
            self.ydata.append(event.ydata)
        else:
            self.xdata.append(event.xdata)
            self.ydata.append(event.ydata)
            self.nClicks = 1
            self.curClicks = 0
            cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.onGinputClick)
            cid = self.w.MplWidget.canvas.mpl_disconnect(cid)
            cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.onGinputClick)
            cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
            codeOut = io.StringIO()
            codeErr = io.StringIO()
            sys.stdout = codeOut
            sys.stderr = codeErr
            print("x-values: {} / xDiff [ppm]: {} / xDiff [Hz]: {}".format(self.xdata, np.abs(np.diff(self.xdata)),
                                                                           np.abs(np.diff(self.xdata)) *
                                                                           self.nd.nmrdat[self.nd.s][self.nd.e].acq.sfo1))
            if self.nd.nmrdat[self.nd.s][self.nd.e].dim == 1:
                print("y-values: {} / yDiff: {}".format(self.ydata, -np.diff(self.ydata)))
            else:
                print("y-values: {} / yDiff: {} / yDiff [Hz]: {}".format(self.ydata, np.abs(np.diff(self.ydata)),
                                                                         np.abs(np.diff(self.ydata)) *
                                                                         self.nd.nmrdat[self.nd.s][self.nd.e].acq.sfo2))

            self.w.console.append(codeOut.getvalue())
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__
            codeOut.close()
            codeErr.close()
            self.xdata = []
            self.ydata = []
            self.showConsole()

    def onGinput2dClick(self, event):
        self.curClicks += 1
        if self.curClicks < self.nClicks:
            self.xdata.append(event.xdata)
            self.ydata.append(event.ydata)
        else:
            self.xdata.append(event.xdata)
            self.ydata.append(event.ydata)
            self.nClicks = 1
            self.curClicks = 0
            cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.onGinput2dClick)
            cid = self.w.MplWidget.canvas.mpl_disconnect(cid)
            cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.onGinput2dClick)
            cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
            self.xy = []
            self.xy = np.resize(self.xy,(1,2))
            self.xy[0][0] = self.xdata[0]
            self.xy[0][1] = self.ydata[0]
            self.xdata = []
            self.ydata = []
            xy = self.xy
            self.showPhCorr2d()
            xyPts = []
            xy2 = []
            xyPts.append(self.nd.nmrdat[self.nd.s][self.nd.e].ppm2points(xy[0][0], 0))
            xyPts.append(self.nd.nmrdat[self.nd.s][self.nd.e].ppm2points(xy[0][1], 1))
            self.phCorr.spcRowPts.append(xyPts[1])
            self.phCorr.spcColPts.append(xyPts[0])
            xy2.append(self.nd.nmrdat[self.nd.s][self.nd.e].points2ppm(xyPts[0], 0))
            xy2.append(self.nd.nmrdat[self.nd.s][self.nd.e].points2ppm(xyPts[1], 1))
            self.phCorr.spcRow.append(xy2[1])
            self.phCorr.spcCol.append(xy2[0])
            self.plot2dColRow()

    def ginput(self, nClicks=1):
        self.w.MplWidget.canvas.setFocus()
        self.showNMRSpectrum()
        self.nClicks = nClicks
        cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.onGinputClick)
        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.onGinputClick)
        cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
        # end ginput

    def ginput2d(self):
        self.w.MplWidget.canvas.setFocus()
        self.showNMRSpectrum()
        self.nClicks = 1
        cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.onGinput2dClick)
        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.onGinput2dClick)
        cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
        # end ginput2d

    def h(self):
        print("Command history: ")
        print(">>><<<")
        for k in range(len(self.nd.cmdBuffer)):
            print(self.nd.cmdBuffer[k])

        return (">>><<<")
        # end h

    def hidePreProcessing(self):
        self.w.preProcessingGroupBox.setHidden(True)
        self.w.preProcessingSelect.setHidden(True)
        self.w.preProcessingWidget.setHidden(True)
        self.w.runPreProcessingButton.setHidden(True)
        self.w.resetPreProcessingButton.setHidden(True)
        self.w.writeScriptButton.setHidden(True)
        self.plotSpc(True)
        # end hidePreProcessing

    def hilbert(self, mat):
        npts = len(mat[0])
        npts1 = len(mat)
        v1 = np.ones(npts1)
        mat1 = np.array([[]], dtype='complex')
        mat1 = np.resize(mat1, (npts1, npts))
        bMat = np.zeros(int(2 * npts), dtype='complex')
        bMat[:(npts + 1)] = np.ones(npts + 1)
        bMat[1:npts] += bMat[1:npts]
        zMat = np.zeros(int(2 * npts), dtype='complex')
        bMat = np.outer(v1, bMat)
        zMat = np.outer(v1, zMat)
        zMat[:, :npts] = mat
        zMat = np.fft.ifft(bMat * np.fft.fft(zMat))
        mat = zMat[:, :npts]
        return mat
        # end hilbert

    def horzPhCorr2d(self):
        s = self.nd.s
        e = self.nd.e
        self.phCorr.nDims = 2
        self.phCorr.dim = 0
        nLines = len(self.phCorr.spcRowPts)
        if nLines > 0:
            npts0 = len(self.nd.nmrdat[s][e].spc)
            npts = len(self.nd.nmrdat[s][e].spc[0])
            self.phCorr.spc = np.zeros((nLines, npts), dtype='complex')
            spc1 = np.copy(self.nd.nmrdat[s][e].spc)
            for k in range(nLines):
                spc = np.array([spc1[npts0 - self.phCorr.spcRowPts[k]]])
                spc = self.hilbert(spc)
                self.phCorr.spc[k] = spc[0]

            self.phCorr.ppm = self.nd.nmrdat[s][e].ppm1
            if self.phCorr.pivotPoints2d[0] < 0:
                self.phCorr.pivotPoints2d[0] = int(len(self.phCorr.ppm) / 2)
                self.phCorr.pivot2d[0] = self.nd.nmrdat[s][e].points2ppm(self.phCorr.pivotPoints2d[0], 0)

        self.showPhCorr2d_1d(self.phCorr.dim)
        self.phCorr.spcMax = np.max(np.max(np.abs(self.phCorr.spc)))
        zwo = False
        pwo = False
        if self.zoomWasOn == True:
            try:
                zwo = True
                self.w.MplWidget.canvas.figure.canvas.toolbar.zoom()
            except:
                pass

            self.setZoomOff()

        if self.w.panWasOn == True:
            try:
                pwo = True
                self.w.MplWidget.canvas.figure.canvas.toolbar.pan()
            except:
                pass

        self.phCorr.maxPh0 = 90.0
        self.phCorr.maxPh1 = 90.0
        cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.onPhCorrClick2d)
        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.onPhCorrRelease2d)
        #self.w.actionApplyPhCorr.triggered.connect(self.apply2dPhCorr)
        #self.w.actionCancelPhCorr.triggered.connect(self.cancel2dPhCorr)
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
        self.phCorrPlotSpc2d(False)
        self.showAcquisitionParameters()
        self.showNMRSpectrum()
        # end horzPhCorr2d

    def loadButton(self):
        selectedDirectory = QFileDialog.getExistingDirectory()
        if (len(selectedDirectory) > 0):
            self.clear()
            self.zeroScript()
        else:
            return

        self.loadFile(selectedDirectory)
        # end saveButton

    def loadConfig(self):
        self.cf.readConfig()
        self.w.phRefColour.setCurrentIndex(self.nd.nmrdat[0][0].display.colours2.get(self.cf.phaseReferenceColour))
        self.w.autoPlot.setChecked(self.cf.autoPlot)
        self.w.keepZoom.setChecked(self.cf.keepZoom)
        self.w.fontSize.setValue(self.cf.fontSize)
        self.stdPosCol1 = (self.cf.posCol10,self.cf.posCol11,self.cf.posCol12)
        self.stdNegCol1 = (self.cf.negCol10,self.cf.negCol11,self.cf.negCol12)
        self.stdPosCol2 = (self.cf.posCol20,self.cf.posCol21,self.cf.posCol22)
        self.stdNegCol2 = (self.cf.negCol20,self.cf.negCol21,self.cf.negCol22)
        self.setStandardColours()
        # end loadConfig

    def loadExampleScript(self):
        idx = self.w.exampleScripts.view().selectedIndexes()[0].row()
        self.w.exampleScripts.setCurrentIndex(idx)
        if (idx == 0):
            fName = os.path.join(os.path.dirname(__file__), "exampleScripts", "example1DScript.py")

        if (idx == 1):
            fName = os.path.join(os.path.dirname(__file__), "exampleScripts", "exampleAutoPhaseScript.py")

        if (idx == 2):
            fName = os.path.join(os.path.dirname(__file__), "exampleScripts", "example2DJresScript.py")

        if (idx == 3):
            fName = os.path.join(os.path.dirname(__file__), "exampleScripts", "examplePreprocessingScript.py")

        if (idx == 4):
            fName = os.path.join(os.path.dirname(__file__), "exampleScripts", "example2DNMRPipeScript.py")

        f = open(fName, 'r')
        scriptText = f.read()
        self.w.script.setText(scriptText)
        # end loadExampleScript

    def loadDarkMode(self):
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            self.nd.nmrdat[self.nd.s][k].display.posColRGB = self.stdPosCol2
            self.nd.nmrdat[self.nd.s][k].display.negColRGB = self.stdNegCol2

        idx = self.w.helpComboBox.currentIndex()
        url = []
        fName = os.path.join(os.path.dirname(__file__), "nmr", "web", "introductionDark", "index.html")
        url.append("file:///" + fName.replace('\\', '/'))
        url.append("https://www.hmdb.ca")
        url.append("https://www.smpdb.ca")
        url.append("https://bmrb.io/metabolomics/")
        url.append("https://www.genome.jp/kegg/pathway.html#metabolism")
        url.append("https://nmrshiftdb.nmr.uni-koeln.de")
        url.append("https://sdbs.db.aist.go.jp/sdbs/cgi-bin/cre_index.cgi")
        url.append("http://dmar.riken.jp/spincouple/")
        self.w.helpView.setUrl(url[idx])
        bg     = ( 42/255,  42/255,  42/255)
        fg     = (255/255, 255/255, 255/255)
        self.w.MplWidget.canvas.figure.set_facecolor(bg)
        self.w.MplWidget.canvas.axes.set_facecolor(bg)
        self.w.MplWidget.canvas.axes.xaxis.label.set_color(fg)
        self.w.MplWidget.canvas.axes.yaxis.label.set_color(fg)
        self.w.MplWidget.canvas.axes.tick_params(axis = 'x', colors = fg)
        self.w.MplWidget.canvas.axes.tick_params(axis = 'y', colors = fg)
        self.w.MplWidget.canvas.axes.spines['bottom'].set_color(fg)
        self.w.MplWidget.canvas.axes.spines['top'].set_color(fg)
        self.w.MplWidget.canvas.axes.spines['left'].set_color(fg)
        self.w.MplWidget.canvas.axes.spines['right'].set_color(fg)
        # end loadDarkMode

    def loadFile(self, fileName):
        self.nd.load(fileName)
        self.w.script.insertHtml(self.nd.script)
        self.w.console.insertHtml(self.nd.console)
        self.resetPlot()
        self.updateGUI()
        self.w.console.verticalScrollBar().setValue(self.w.console.verticalScrollBar().maximum())
        self.showTitleFileInformation()
        self.showAcquisitionParameters()
        self.showNMRSpectrum()
        # end loadFile

    def loadLightMode(self):
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            self.nd.nmrdat[self.nd.s][k].display.posColRGB = self.stdPosCol1
            self.nd.nmrdat[self.nd.s][k].display.negColRGB = self.stdNegCol1

        idx = self.w.helpComboBox.currentIndex()
        url = []
        fName = os.path.join(os.path.dirname(__file__), "nmr", "web", "introduction", "index.html")
        url.append("file:///" + fName.replace('\\', '/'))
        url.append("https://www.hmdb.ca")
        url.append("https://www.smpdb.ca")
        url.append("https://bmrb.io/metabolomics/")
        url.append("https://www.genome.jp/kegg/pathway.html#metabolism")
        url.append("https://nmrshiftdb.nmr.uni-koeln.de")
        url.append("https://sdbs.db.aist.go.jp/sdbs/cgi-bin/cre_index.cgi")
        url.append("http://dmar.riken.jp/spincouple/")
        self.w.helpView.setUrl(url[idx])
        bg     = (255/255, 255/255, 255/255)
        fg     = (  0/255,   0/255,   0/255)
        self.w.MplWidget.canvas.figure.set_facecolor(bg)
        self.w.MplWidget.canvas.axes.set_facecolor(bg)
        self.w.MplWidget.canvas.axes.xaxis.label.set_color(fg)
        self.w.MplWidget.canvas.axes.yaxis.label.set_color(fg)
        self.w.MplWidget.canvas.axes.tick_params(axis = 'x', colors = fg)
        self.w.MplWidget.canvas.axes.tick_params(axis = 'y', colors = fg)
        self.w.MplWidget.canvas.axes.spines['bottom'].set_color(fg)
        self.w.MplWidget.canvas.axes.spines['top'].set_color(fg)
        self.w.MplWidget.canvas.axes.spines['left'].set_color(fg)
        self.w.MplWidget.canvas.axes.spines['right'].set_color(fg)
        # end loadLightMode

    def nextCommand(self):
        if (self.w.cmdLine.hasFocus() == True):
            if (self.nd.cmdIdx < len(self.nd.cmdBuffer)):
                self.nd.cmdIdx += 1
                if (self.nd.cmdIdx == len(self.nd.cmdBuffer)):
                    self.w.cmdLine.setText("")
                else:
                    self.w.cmdLine.setText(self.nd.cmdBuffer[self.nd.cmdIdx])

        # end nextCommand

    def nextTab(self):
        cidx = self.w.nmrSpectrum.currentIndex()
        while self.w.nmrSpectrum.isTabEnabled(cidx + 1) is False and cidx < 11:
            cidx += 1

        if cidx < 12:
            self.w.nmrSpectrum.setCurrentIndex(cidx + 1)
        else:
            self.w.nmrSpectrum.setCurrentIndex(0)


        #end nextTab

    def previousTab(self):
        cidx = self.w.nmrSpectrum.currentIndex()
        while self.w.nmrSpectrum.isTabEnabled(cidx - 1) is False and cidx > 1:
            cidx -= 1

        if cidx > 0:
            self.w.nmrSpectrum.setCurrentIndex(cidx - 1)
            self.w.nmrSpectrum.setFocus()
        else:
            self.w.nmrSpectrum.setCurrentIndex(12)

        # end previousTab

    def onPhCorrClick(self, event):
        s = self.nd.s
        e = self.nd.e
        if (self.zoom == False):
            self.phCorr.spc = self.nd.nmrdat[s][e].spc
            self.phCorr.spcMax = max(max(abs(self.phCorr.spc)))
            #self.w.MplWidget.canvas.toolbar._zoom_mode.__init__()
            if (event.button == 1):
                mods = QApplication.queryKeyboardModifiers()
                if (mods == QtCore.Qt.ControlModifier):
                    # set pivot for phase correction
                    self.phCorr.start = event.xdata
                    self.phCorr.pivot = event.xdata
                    self.phCorr.pivPoints = self.nd.nmrdat[s][e].ppm2points(self.phCorr.pivot, 0)

                if (mods == QtCore.Qt.ShiftModifier):
                    # first order phase correction
                    self.phCorr.start = event.ydata

                if (mods == QtCore.Qt.NoModifier):
                    # zero order phase correction
                    self.phCorr.start = event.ydata

                if (mods == QtCore.Qt.AltModifier):
                    self.w.MplWidget.canvas.manager.toolbar.zoom()

            else:
                if (event.button == 2):
                    # set pivot for phase correction
                    self.phCorr.start = event.xdata
                    self.phCorr.pivot = event.xdata
                    self.phCorr.pivPoints = self.nd.nmrdat[s][e].ppm2points(self.phCorr.pivot, 0)
                else:
                    # first order phase correction
                    self.phCorr.start = event.ydata

            cid3 = self.w.MplWidget.canvas.mpl_connect('motion_notify_event', self.onPhCorrDraw)

        # end onPhCorrClick

    def onPhCorrClick2d(self, event):
        s = self.nd.s
        e = self.nd.e
        if (self.zoom == False):
            self.phCorr.spc2 = np.copy(self.phCorr.spc)
            # self.phCorr.spc = self.nd.nmrdat[s][e].spc
            # self.phCorr.spcMax = max(max(abs(self.phCorr.spc)))
            if (event.button == 1):
                mods = QApplication.queryKeyboardModifiers()
                if (mods == QtCore.Qt.ControlModifier):
                    # set pivot for phase correction
                    self.phCorr.start = event.xdata
                    self.phCorr.pivot = event.xdata
                    self.phCorr.pivotPoints2d[self.phCorr.dim] = self.nd.nmrdat[s][e].ppm2points(
                        self.phCorr.pivot2d[self.phCorr.dim], self.phCorr.dim)

                if (mods == QtCore.Qt.ShiftModifier):
                    # first order phase correction
                    self.phCorr.start = event.ydata

                if (mods == QtCore.Qt.NoModifier):
                    # zero order phase correction
                    self.phCorr.start = event.ydata

                if (mods == QtCore.Qt.AltModifier):
                    self.w.MplWidget.canvas.manager.toolbar.zoom()

            else:
                if (event.button == 2):
                    # set pivot for phase correction
                    self.phCorr.start = event.xdata
                    self.phCorr.pivot2d[self.phCorr.dim] = event.xdata
                    self.phCorr.pivotPoints2d[self.phCorr.dim] = self.nd.nmrdat[s][e].ppm2points(
                        self.phCorr.pivot2d[self.phCorr.dim], self.phCorr.dim)
                else:
                    # first order phase correction
                    self.phCorr.start = event.ydata

            cid3 = self.w.MplWidget.canvas.mpl_connect('motion_notify_event', self.onPhCorrDraw2d)

        # end onPhCorrClick2d

    def onPhCorrDraw(self, event):
        if (self.zoom == False):
            s = self.nd.s
            e = self.nd.e
            if ((event.xdata != None) & (event.ydata != None)):
                self.phCorr.xData = event.xdata
                self.phCorr.yData = event.ydata
                if (event.button == 1):
                    mods = QApplication.queryKeyboardModifiers()
                    if (mods == QtCore.Qt.ControlModifier):
                        # set pivot for phase correction
                        self.phCorr.pivot = event.xdata
                        self.phCorr.pivPoints = self.nd.nmrdat[s][e].ppm2points(self.phCorr.pivot, 0)

                    if (mods == QtCore.Qt.ShiftModifier):
                        # first order phase correction
                        ph0 = 0
                        ph1 = self.phCorr.maxPh1 * (event.ydata - self.phCorr.start) / self.phCorr.spcMax
                        self.phCorr.spc = self.phase1d(self.nd.nmrdat[s][e].spc, ph0, ph1, self.phCorr.pivPoints)

                    if (mods == QtCore.Qt.NoModifier):
                        # zero order phase correction
                        ph0 = self.phCorr.maxPh0 * (event.ydata - self.phCorr.start) / self.phCorr.spcMax
                        ph1 = 0
                        self.phCorr.spc = self.phase1d(self.nd.nmrdat[s][e].spc, ph0, ph1, self.phCorr.pivPoints)

                else:
                    if (event.button == 2):
                        # set pivot for phase correction
                        self.phCorr.xData = event.xdata
                        self.phCorr.yData = event.ydata
                        self.phCorr.pivot = event.xdata
                        self.phCorr.pivPoints = self.nd.nmrdat[s][e].ppm2points(self.phCorr.pivot, 0)
                    else:
                        # first order phase correction
                        self.phCorr.xData = event.xdata
                        self.phCorr.yData = event.ydata
                        ph0 = 0
                        ph1 = self.phCorr.maxPh1 * (event.ydata - self.phCorr.start) / self.phCorr.spcMax
                        self.phCorr.spc = self.phase1d(self.nd.nmrdat[s][e].spc, ph0, ph1, self.phCorr.pivPoints)

            self.phCorrPlotSpc()

        # end onPhCorrDraw

    def onPhCorrDraw2d(self, event):
        if (self.zoom == False):
            s = self.nd.s
            e = self.nd.e
            if ((event.xdata != None) & (event.ydata != None)):
                self.phCorr.xData = event.xdata
                self.phCorr.yData = event.ydata
                if (event.button == 1):
                    mods = QApplication.queryKeyboardModifiers()
                    if (mods == QtCore.Qt.ControlModifier):
                        # set pivot for phase correction
                        self.phCorr.pivot2d[self.phCorr.dim] = event.xdata
                        self.phCorr.pivotPoints2d[self.phCorr.dim] = self.nd.nmrdat[s][e].ppm2points(self.phCorr.pivot2d[self.phCorr.dim], self.phCorr.dim)

                    if (mods == QtCore.Qt.ShiftModifier):
                        # first order phase correction
                        ph0 = 0
                        ph1 = self.phCorr.maxPh1 * (event.ydata - self.phCorr.start) / self.phCorr.spcMax
                        self.phCorr.spc = self.phase1d(self.phCorr.spc2, ph0, ph1, self.phCorr.pivotPoints2d[self.phCorr.dim])

                    if (mods == QtCore.Qt.NoModifier):
                        # zero order phase correction
                        ph0 = self.phCorr.maxPh0 * (event.ydata - self.phCorr.start) / self.phCorr.spcMax
                        ph1 = 0
                        self.phCorr.spc = self.phase1d(self.phCorr.spc2, ph0, ph1, self.phCorr.pivotPoints2d[self.phCorr.dim])

                else:
                    if (event.button == 2):
                        # set pivot for phase correction
                        self.phCorr.xData = event.xdata
                        self.phCorr.yData = event.ydata
                        self.phCorr.pivot2d[self.phCorr.dim] = event.xdata
                        self.phCorr.pivotPoints2d[self.phCorr.dim] = self.nd.nmrdat[s][e].ppm2points(
                            self.phCorr.pivot2d[self.phCorr.dim], self.phCorr.dim)
                    else:
                        # first order phase correction
                        self.phCorr.xData = event.xdata
                        self.phCorr.yData = event.ydata
                        ph0 = 0
                        ph1 = self.phCorr.maxPh1 * (event.ydata - self.phCorr.start) / self.phCorr.spcMax
                        self.phCorr.spc = self.phase1d(self.phCorr.spc2, ph0, ph1,
                                                       self.phCorr.pivotPoints2d[self.phCorr.dim])

            self.phCorrPlotSpc2d()

        # end onPhCorrDrawd

    def onPhCorrRelease(self, event):
        s = self.nd.s
        e = self.nd.e
        if ((event.xdata != None) & (event.ydata != None)):
            xdata = event.xdata
            ydata = event.ydata
        else:
            xdata = self.phCorr.xData
            ydata = self.phCorr.yData

        if (self.zoom == False):
            if (event.button == 1):
                mods = QApplication.queryKeyboardModifiers()
                if (mods == QtCore.Qt.ControlModifier):
                    # set pivot for phase correction
                    self.phCorr.pivot = xdata
                    self.phCorr.pivPoints = self.nd.nmrdat[s][e].ppm2points(self.phCorr.pivot, 0)

                if (mods == QtCore.Qt.ShiftModifier):
                    # first order phase correction
                    ph1 = (self.phCorr.maxPh1 * (ydata - self.phCorr.start) / self.phCorr.spcMax)
                    ph = self.phasesRemovePivot(0.0, ph1, self.phCorr.pivPoints, len(self.phCorr.spc[0]))
                    ph0 = ((self.nd.nmrdat[s][e].proc.ph0[0] + ph[0] + 180.0) % 360.0) - 180.0
                    ph1 = self.nd.nmrdat[s][e].proc.ph1[0] + ph[1]
                    self.nd.nmrdat[s][e].proc.ph0[0] = ph0
                    self.nd.nmrdat[s][e].proc.ph1[0] = ph1

                if (mods == QtCore.Qt.NoModifier):
                    # zero order phase correction
                    ph0a = (self.phCorr.maxPh0 * (ydata - self.phCorr.start) / self.phCorr.spcMax) % 360.0
                    ph1a = 0.0
                    ph = self.phasesRemovePivot(ph0a, ph1a, self.phCorr.pivPoints, len(self.phCorr.spc[0]))
                    ph0 = ((self.nd.nmrdat[s][e].proc.ph0[0] + ph[0] + 180.0) % 360.0) - 180.0
                    ph1 = self.nd.nmrdat[s][e].proc.ph1[0] + ph[1]
                    self.nd.nmrdat[s][e].proc.ph0[0] = ph0
                    self.nd.nmrdat[s][e].proc.ph1[0] = ph1

            else:
                if (event.button == 2):
                    # set pivot for phase correction
                    self.phCorr.pivot = xdata
                    self.phCorr.pivPoints = self.nd.nmrdat[s][e].ppm2points(self.phCorr.pivot, 0)
                else:
                    # first order phase correction
                    ph1 = (self.phCorr.maxPh1 * (ydata - self.phCorr.start) / self.phCorr.spcMax)
                    ph = self.phasesRemovePivot(0.0, ph1, self.phCorr.pivPoints, len(self.phCorr.spc[0]))
                    ph0 = ((self.nd.nmrdat[s][e].proc.ph0[0] + ph[0] + 180.0) % 360.0) - 180.0
                    ph1 = self.nd.nmrdat[s][e].proc.ph1[0] + ph[1]
                    self.nd.nmrdat[s][e].proc.ph0[0] = ph0
                    self.nd.nmrdat[s][e].proc.ph1[0] = ph1

            cid3 = self.w.MplWidget.canvas.mpl_connect('motion_notify_event', self.onPhCorrDraw)
            cid3 = self.w.MplWidget.canvas.mpl_disconnect(cid3)
            self.nd.nmrdat[s][e].spc = self.phCorr.spc
            self.setProcPars()
            self.nd.ft()
            self.phCorrPlotSpc()
        else:
            # zoom mode activated
            if (event.button > 1):
                # Right MB click will unzoom the plot
                try:
                    self.w.MplWidget.canvas.figure.canvas.toolbar.home()
                except:
                    pass

        # end onPhCorrRelease

    def onPhCorrRelease2d(self, event):
        s = self.nd.s
        e = self.nd.e
        if ((event.xdata != None) & (event.ydata != None)):
            xdata = event.xdata
            ydata = event.ydata
        else:
            xdata = self.phCorr.xData
            ydata = self.phCorr.yData

        if (self.zoom == False):
            if event.button == 1:
                mods = QApplication.queryKeyboardModifiers()
                if mods == QtCore.Qt.ControlModifier:
                    # set pivot for phase correction
                    self.phCorr.pivot2d[self.phCorr.dim] = xdata
                    self.phCorr.pivotPoints2d[self.phCorr.dim] = self.nd.nmrdat[s][e].ppm2points(self.phCorr.pivot2d[self.phCorr.dim], self.phCorr.dim)

                if mods == QtCore.Qt.ShiftModifier:
                    # first order phase correction
                    ph1 = (self.phCorr.maxPh1 * (ydata - self.phCorr.start) / self.phCorr.spcMax)
                    ph = self.phasesRemovePivot(0.0, ph1, self.phCorr.pivotPoints2d[self.phCorr.dim], len(self.phCorr.spc[0]))
                    ph0 = ((self.phCorr.ph0_2d[self.phCorr.dim] + ph[0] + 180.0) % 360.0) - 180.0
                    ph1 = self.phCorr.ph1_2d[self.phCorr.dim] + ph[1]
                    self.phCorr.ph0_2d[self.phCorr.dim] = ph0
                    self.phCorr.ph1_2d[self.phCorr.dim] = ph1

                if mods == QtCore.Qt.NoModifier:
                    # zero order phase correction
                    ph0a = (self.phCorr.maxPh0 * (ydata - self.phCorr.start) / self.phCorr.spcMax) % 360.0
                    ph1a = 0.0
                    ph = self.phasesRemovePivot(ph0a, ph1a, self.phCorr.pivotPoints2d[self.phCorr.dim], len(self.phCorr.spc[0]))
                    ph0 = ((self.phCorr.ph0_2d[self.phCorr.dim] + ph[0] + 180.0) % 360.0) - 180.0
                    ph1 = self.phCorr.ph1_2d[self.phCorr.dim] + ph[1]
                    self.phCorr.ph0_2d[self.phCorr.dim] = ph0
                    self.phCorr.ph1_2d[self.phCorr.dim] = ph1

            else:
                if event.button == 2:
                    # set pivot for phase correction
                    self.phCorr.pivot2d[self.phCorr.dim] = xdata
                    self.phCorr.pivotPoints2d[self.phCorr.dim] = self.nd.nmrdat[s][e].ppm2points(self.phCorr.pivot2d[self.phCorr.dim], self.phCorr.dim)

                else:
                    # first order phase correction
                    ph1 = (self.phCorr.maxPh1 * (ydata - self.phCorr.start) / self.phCorr.spcMax)
                    ph = self.phasesRemovePivot(0.0, ph1, self.phCorr.pivotPoints2d[self.phCorr.dim],len(self.phCorr.spc[0]))
                    ph0 = ((self.phCorr.ph0_2d[self.phCorr.dim] + ph[0] + 180.0) % 360.0) - 180.0
                    ph1 = self.phCorr.ph1_2d[self.phCorr.dim] + ph[1]
                    self.phCorr.ph0_2d[self.phCorr.dim] = ph0
                    self.phCorr.ph1_2d[self.phCorr.dim] = ph1

            cid3 = self.w.MplWidget.canvas.mpl_connect('motion_notify_event', self.onPhCorrDraw2d)
            cid3 = self.w.MplWidget.canvas.mpl_disconnect(cid3)
            self.phCorr.spc2 = np.copy(self.phCorr.spc)
            self.phCorrPlotSpc2d()
        else:
            # zoom mode activated
            if (event.button > 1):
                # Right MB click will unzoom the plot
                try:
                    self.w.MplWidget.canvas.figure.canvas.toolbar.home()
                except:
                    pass

        # end onPhCorrRelease2d

    def openScript(self, fName=""):
        if (fName == False):
            fName = ""

        if (len(fName) == 0):
            fName = QFileDialog.getOpenFileName(None, 'Open Script File', '', 'Python files (*.py)')
            fName = fName[0]

        if (len(fName) > 0):
            f = open(fName, 'r')
            scriptText = f.read()
            self.w.script.setText(scriptText)

        self.w.nmrSpectrum.setCurrentIndex(10)
        # end openScript

    def phCorrPlotSpc(self):
        xlim = self.w.MplWidget.canvas.axes.get_xlim()
        ylim = self.w.MplWidget.canvas.axes.get_ylim()
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        if (d.posCol == "RGB"):
            posCol = d.posColRGB
        else:
            posCol = d.posCol

        if (d.negCol == "RGB"):
            negCol = d.negColRGB
        else:
            negCol = d.negCol

        refCol = d.phRefCol
        posCol = matplotlib.colors.to_hex(posCol)
        negCol = matplotlib.colors.to_hex(negCol)
        refCol = matplotlib.colors.to_hex(refCol)
        xlabel = d.xLabel + " [" + d.axisType1 + "]"
        ylabel = d.yLabel + " [" + d.axisType2 + "]"
        if (self.nd.nmrdat[self.nd.s][self.nd.e].dim == 1):
            self.w.MplWidget.canvas.axes.clear()
            if ((d.phRefDS > 0) & (d.phRefExp > 0) & (
                    ((d.phRefDS - 1 == self.nd.s) & (d.phRefExp - 1 == self.nd.e)) is False)):
                self.w.MplWidget.canvas.axes.plot(self.nd.nmrdat[d.phRefDS - 1][d.phRefExp - 1].ppm1,
                                                  self.nd.nmrdat[d.phRefDS - 1][d.phRefExp - 1].spc[0].real,
                                                  color=refCol)

            self.w.MplWidget.canvas.axes.plot(self.nd.nmrdat[self.nd.s][self.nd.e].ppm1, self.phCorr.spc[0].real,
                                              color=posCol)
            self.w.MplWidget.canvas.axes.plot([self.phCorr.pivot, self.phCorr.pivot],
                                              [2.0 * self.phCorr.spcMax, -2.0 * self.phCorr.spcMax], color='r')
            self.w.MplWidget.canvas.axes.set_xlabel(xlabel)
            self.w.MplWidget.canvas.axes.invert_xaxis()
            self.w.MplWidget.canvas.axes.set_xlim(xlim)
            self.w.MplWidget.canvas.axes.set_ylim(ylim)

        self.setProcPars()
        self.w.MplWidget.canvas.draw()
        # This is a messy solution to force the matplotlib widget to update the plot by introducing an error (calling
        # a figure object and redirecting the error output
        codeErr = io.StringIO()
        sys.stderr = codeErr
        try:
            self.w.MplWidget.canvas.figure()
        except:
            pass

        sys.stderr = sys.__stderr__
        # end phCorrPlotSpc

    def phCorrPlotSpc2d(self, keepZoom = True):
        if keepZoom:
            xlim = self.w.MplWidget.canvas.axes.get_xlim()
            ylim = self.w.MplWidget.canvas.axes.get_ylim()

        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        self.w.MplWidget.canvas.axes.set_prop_cycle(None)
        if self.phCorr.dim == 0:
            xlabel = d.xLabel + " [" + d.axisType1 + "]"
        else:
            xlabel = d.yLabel + " [" + d.axisType2 + "]"

        self.w.MplWidget.canvas.axes.clear()
        self.phCorr.spcMax = 0.0
        for k in range(len(self.phCorr.spc)):
            self.w.MplWidget.canvas.axes.plot(self.phCorr.ppm, self.phCorr.spc[k].real)
            self.phCorr.spcMax = max(self.phCorr.spcMax, np.max(np.max(np.abs(self.phCorr.spc[k].real))))

        self.w.MplWidget.canvas.axes.invert_xaxis()
        if not keepZoom:
            xlim = self.w.MplWidget.canvas.axes.get_xlim()
            ylim = self.w.MplWidget.canvas.axes.get_ylim()

        self.w.MplWidget.canvas.axes.plot([self.phCorr.pivot2d[self.phCorr.dim], self.phCorr.pivot2d[self.phCorr.dim]],
                                          [2.0 * self.phCorr.spcMax, -2.0 * self.phCorr.spcMax], color='r')
        self.w.MplWidget.canvas.axes.set_xlabel(xlabel)
        self.w.MplWidget.canvas.axes.invert_xaxis()
        self.w.MplWidget.canvas.axes.set_xlim(xlim)
        self.w.MplWidget.canvas.axes.set_ylim(ylim)
        self.w.MplWidget.canvas.draw()
        # This is a messy solution to force the matplotlib widget to update the plot by introducing an error (calling
        # a figure object and redirecting the error output
        codeErr = io.StringIO()
        sys.stderr = codeErr
        try:
            self.w.MplWidget.canvas.figure()
        except:
            pass

        sys.stderr = sys.__stderr__
        # end phCorrPlotSpc2d

    def phase1d(self, mat, ph0, ph1, piv):
        npts = len(mat[0])
        ph0 = -ph0 * math.pi / 180.0
        ph1 = -ph1 * math.pi / 180.0
        frac = np.linspace(0, 1, npts) - float(npts - piv) / float(npts)
        ph = ph0 + frac * ph1
        mat = np.cos(ph) * mat.real + np.sin(ph) * mat.imag + 1j * (-np.sin(ph) * mat.real + np.cos(ph) * mat.imag)
        return mat
        # end phase1d

    def phasesRemovePivot(self, phc0, phc1, piv, npts):
        phases = np.array([0.0, 0.0])
        frac = np.linspace(0, 1, npts) - float(npts - piv) / float(npts)
        ph = -phc0 - frac * phc1
        phases[0] = -ph[0]
        phases[1] = ph[0] - ph[len(ph) - 1]
        return phases
        # end phasesRemovePivot

    def pickColRow(self):
        self.w.statusBar().clearMessage()
        self.w.statusBar().showMessage("Click to add row/col")
        self.showAcquisitionParameters()
        self.showNMRSpectrum()
        self.ginput2d()
        # end pickColRow

    def plot2dColRow(self):
        while len(self.w.MplWidget.canvas.axes.lines) > 0:
            self.w.MplWidget.canvas.axes.lines[0].remove()

        self.w.MplWidget.canvas.axes.set_prop_cycle(None)
        ppm1 = self.nd.nmrdat[self.nd.s][self.nd.e].ppm1
        ppm2 = self.nd.nmrdat[self.nd.s][self.nd.e].ppm2
        for k in range(len(self.phCorr.spcRow)):
            pid = self.w.MplWidget.canvas.axes.plot([self.phCorr.spcCol[k], self.phCorr.spcCol[k]],
                                                    [np.min(ppm2), np.max(ppm2)])
            self.w.MplWidget.canvas.axes.plot([np.min(ppm1), np.max(ppm1)],
                                              [self.phCorr.spcRow[k], self.phCorr.spcRow[k]], color=pid[0].get_color())

        self.w.MplWidget.canvas.draw()

    def plotSpc(self, hidePreProcessing = False):
        self.keepZoom = self.w.keepZoom.isChecked()
        xlim = self.w.MplWidget.canvas.axes.get_xlim()
        ylim = self.w.MplWidget.canvas.axes.get_ylim()
        self.w.nmrSpectrum.setCurrentIndex(0)
        if (len(self.nd.nmrdat[self.nd.s]) == 0):
            return

        if (len(self.nd.nmrdat[self.nd.s][self.nd.e].spc) == 0):
            return

        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        if (d.posCol == "RGB"):
            posCol = d.posColRGB
        else:
            posCol = d.posCol

        if (d.negCol == "RGB"):
            negCol = d.negColRGB
        else:
            negCol = d.negCol

        posCol = matplotlib.colors.to_hex(posCol)
        negCol = matplotlib.colors.to_hex(negCol)
        xlabel = d.xLabel + " [" + d.axisType1 + "]"
        ylabel = d.yLabel + " [" + d.axisType2 + "]"
        if (self.nd.nmrdat[self.nd.s][self.nd.e].dim == 1):
            self.w.MplWidget.canvas.axes.clear()
            for k in range(len(self.nd.nmrdat[self.nd.s])):
                if ((k != self.nd.e) and (self.nd.nmrdat[self.nd.s][k].display.displaySpc == True)):
                    d = self.nd.nmrdat[self.nd.s][k].display
                    if (d.posCol == "RGB"):
                        posCol = d.posColRGB
                    else:
                        posCol = d.posCol

                    if (d.negCol == "RGB"):
                        negCol = d.negColRGB
                    else:
                        negCol = d.negCol

                    posCol = matplotlib.colors.to_hex(posCol)
                    negCol = matplotlib.colors.to_hex(negCol)
                    self.w.MplWidget.canvas.axes.plot(self.nd.nmrdat[self.nd.s][k].ppm1,
                                                      self.nd.nmrdat[self.nd.s][k].spc[0].real, color=posCol)

            d = self.nd.nmrdat[self.nd.s][self.nd.e].display
            if (d.posCol == "RGB"):
                posCol = d.posColRGB
            else:
                posCol = d.posCol

            if (d.negCol == "RGB"):
                negCol = d.negColRGB
            else:
                negCol = d.negCol

            posCol = matplotlib.colors.to_hex(posCol)
            negCol = matplotlib.colors.to_hex(negCol)
            xlabel = d.xLabel + " [" + d.axisType1 + "]"
            ylabel = d.yLabel + " [" + d.axisType2 + "]"
            self.w.MplWidget.canvas.axes.plot(self.nd.nmrdat[self.nd.s][self.nd.e].ppm1,
                                              self.nd.nmrdat[self.nd.s][self.nd.e].spc[0].real, color=posCol)
            self.w.MplWidget.canvas.axes.set_xlabel(xlabel)
            self.w.MplWidget.canvas.axes.autoscale()
            self.w.MplWidget.canvas.axes.invert_xaxis()
            if (self.keepZoom == True):
                self.w.MplWidget.canvas.axes.set_xlim(xlim)
                self.w.MplWidget.canvas.axes.set_ylim(ylim)

            # self.w.MplWidget.canvas.toolbar.update()
            self.w.MplWidget.canvas.draw()
            if (self.keepXZoom == True):
                self.w.MplWidget.canvas.axes.set_xlim(xlim)
                self.verticalAutoScale()
                self.keepXZoom = False


        else:
            mm = np.max(np.abs(self.nd.nmrdat[self.nd.s][self.nd.e].spc.real))
            posLev = np.linspace(d.minLevel * mm, d.maxLevel * mm, d.nLevels)
            negLev = np.linspace(-d.maxLevel * mm, -d.minLevel * mm, d.nLevels)
            self.w.MplWidget.canvas.axes.clear()
            self.w.MplWidget.canvas.axes.contour(self.nd.nmrdat[self.nd.s][self.nd.e].ppm1,
                                                 self.nd.nmrdat[self.nd.s][self.nd.e].ppm2,
                                                 self.nd.nmrdat[self.nd.s][self.nd.e].spc.real, posLev, colors=posCol,
                                                 linestyles='solid', antialiased=True)
            self.w.MplWidget.canvas.axes.contour(self.nd.nmrdat[self.nd.s][self.nd.e].ppm1,
                                                 self.nd.nmrdat[self.nd.s][self.nd.e].ppm2,
                                                 self.nd.nmrdat[self.nd.s][self.nd.e].spc.real, negLev, colors=negCol,
                                                 linestyles='solid', antialiased=True)
            self.w.MplWidget.canvas.axes.set_xlabel(xlabel)
            self.w.MplWidget.canvas.axes.set_ylabel(ylabel)
            self.w.MplWidget.canvas.axes.autoscale()
            self.w.MplWidget.canvas.axes.invert_xaxis()
            self.w.MplWidget.canvas.axes.invert_yaxis()
            if (self.keepZoom == True):
                self.w.MplWidget.canvas.axes.set_xlim(xlim)
                self.w.MplWidget.canvas.axes.set_ylim(ylim)
            else:
                if (self.keepXZoom == True):
                    self.w.MplWidget.canvas.axes.set_xlim(xlim)
                    self.keepXZoom = False

            # self.w.MplWidget.canvas.toolbar.update()
            self.w.MplWidget.canvas.draw()

        self.keepZoom = False
        if hidePreProcessing == False:
            pyautogui.click(clicks=1)

        # end plotSpc

    def plotSpcDisp(self):
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.changeDataSetExp()
        if (self.phCorrActive == False):
            self.plotSpc()
        else:
            self.phCorrPlotSpc()

        # end plotSpcDisp

    def plotSpcPreProc(self):
        if (len(self.nd.pp.classSelect) == 0):
            self.nd.preProcInit()

        # self.w.rDolphinExport.setChecked(self.nd.pp.rDolphinExport)
        self.fillPreProcessingNumbers()
        sel = self.w.selectClassTW.selectedIndexes()
        cls = np.array([])
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            cls = np.append(cls, self.w.selectClassTW.item(k, 1).text())

        self.nd.pp.classSelect = cls
        cls2 = np.unique(cls)
        sel2 = np.array([], dtype='int')
        for k in range(len(sel)):
            if (sel[k].column() == 0):
                sel2 = np.append(sel2, int(sel[k].row()))

        self.nd.pp.plotSelect = sel2
        self.keepZoom = self.w.keepZoom.isChecked()
        xlim = self.w.MplWidget.canvas.axes.get_xlim()
        ylim = self.w.MplWidget.canvas.axes.get_ylim()
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.w.MplWidget.canvas.axes.clear()
        if (self.w.preProcessingWidget.currentIndex() == 1):
            for k in range(len(self.nd.pp.excludeStart)):
                self.w.MplWidget.canvas.axes.axvspan(self.nd.pp.excludeStart[k], self.nd.pp.excludeEnd[k],
                                                     alpha=self.nd.pp.alpha, color=self.nd.pp.colour)

        if (self.w.preProcessingWidget.currentIndex() == 2):
            for k in range(len(self.nd.pp.segStart)):
                self.w.MplWidget.canvas.axes.axvspan(self.nd.pp.segStart[k], self.nd.pp.segEnd[k],
                                                     alpha=self.nd.pp.alpha, color=self.nd.pp.colour)

        for k in range(len(self.nd.pp.plotSelect)):
            colIdx = np.where(cls2 == cls[self.nd.pp.plotSelect[k]])[0][0]
            plotCol = matplotlib.colors.to_hex(self.nd.pp.plotColours[colIdx])
            self.w.MplWidget.canvas.axes.plot(self.nd.nmrdat[self.nd.s][self.nd.pp.plotSelect[k]].ppm1,
                                              self.nd.nmrdat[self.nd.s][self.nd.pp.plotSelect[k]].spc[0].real,
                                              color=plotCol)

        if (self.w.preProcessingWidget.currentIndex() == 3):
            self.w.MplWidget.canvas.axes.axvspan(self.nd.pp.noiseStart, self.nd.pp.noiseEnd, alpha=self.nd.pp.alpha,
                                                 color=self.nd.pp.colour)
            val = self.nd.pp.noiseThreshold * self.nd.pp.stdVal
            x = [self.nd.nmrdat[self.nd.s][0].ppm1[0], self.nd.nmrdat[self.nd.s][0].ppm1[-1]]
            y = [val, val]
            self.w.MplWidget.canvas.axes.plot(x, y, color=self.nd.pp.thColour, linewidth=self.nd.pp.thLineWidth)

        if (self.w.preProcessingWidget.currentIndex() == 5):
            for k in range(len(self.nd.pp.compressStart)):
                self.w.MplWidget.canvas.axes.axvspan(self.nd.pp.compressStart[k], self.nd.pp.compressEnd[k],
                                                     alpha=self.nd.pp.alpha, color=self.nd.pp.colour)

        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        xlabel = d.xLabel + " [" + d.axisType1 + "]"
        self.w.MplWidget.canvas.axes.set_xlabel(xlabel)
        self.w.MplWidget.canvas.axes.autoscale()
        self.w.MplWidget.canvas.axes.invert_xaxis()
        if (self.keepZoom == True):
            self.w.MplWidget.canvas.axes.set_xlim(xlim)
            self.w.MplWidget.canvas.axes.set_ylim(ylim)

        # self.w.MplWidget.canvas.toolbar.update()
        self.w.MplWidget.canvas.draw()

    def previousCommand(self):
        if (self.w.cmdLine.hasFocus() == True):
            if (self.nd.cmdIdx > 0):
                self.nd.cmdIdx -= 1
                self.w.cmdLine.setText(self.nd.cmdBuffer[self.nd.cmdIdx])

        # end previousCommand

    def quit_app(self):
        # some actions to perform before actually quitting:
        self.w.close()
        # end quit_app

    def readNMRSpc(self):
        kz = self.w.keepZoom.isChecked()
        if (len(self.nd.nmrdat[0]) == 0):
            self.w.keepZoom.setChecked(False)

        selected_directory = QFileDialog.getExistingDirectory()
        if (len(selected_directory) > 0):
            # Use the selected directory...
            idx = selected_directory.rfind('/')
            dsName = selected_directory[:idx]
            expName = selected_directory[idx + 1:]
            self.nd.readSpc(dsName, expName)
            self.setJres()
            self.nd.ft()
            if self.nd.nmrdat[self.nd.s][self.nd.e].dim == 0:
                self.nd.autoref(True)
            else:
                self.nd.autoref(False)

            self.nd.e = len(self.nd.nmrdat[self.nd.s]) - 1
            self.plotSpc()
            self.w.keepZoom.setChecked(kz)
            self.setProcPars()
            self.setAcqPars()
            self.setTitleFile()
            self.setPulseProgram()
            self.w.expBox.setValue(self.nd.e + 1)
            self.setDispPars()
            self.updateGUI()

        # end readNMRSpc

    def readNMRPipeSpc(self, sfile=False):
        if sfile == False:
            selectedFile = QFileDialog.getOpenFileName()
            if len(selectedFile[0]) == 0:
                return

        else:
            selectedFile = (sfile, '')

        #print(selectedFile)
        fName = os.path.split(selectedFile[0])[1]
        dataPath = os.path.split(os.path.split(selectedFile[0])[0])[0]
        expNum = os.path.split(os.path.split(selectedFile[0])[0])[1]
        if expNum.find('.') > -1:
            expNum = expNum[:expNum.find('.')]

        self.readNMRPipeSpcs([dataPath], [expNum], fName)
        self.setStandardColours()
        self.updateGUI()
        self.resetPlot()
        # end readNMRPipeSpc

    def readNMRPipeSpcs(self, dataPath, dataSets, procDataName='test.dat'):
        zFill = 25
        if (dataPath[0] == 'interactive'):
            dataPath = [QFileDialog.getExistingDirectory()]

        if (len(dataPath) > 0):
            if (str(dataSets) == 'all'):
                folders = []
                for r, d, f in os.walk(dataPath):
                    for folder in d:
                        if (os.path.isfile(os.path.join(r, folder, procDataName))):
                            folders.append(folder.zfill(zFill).rstrip('.proc'))

                folders.sort()
                dataSets = []
                for k in range(len(folders)):
                    dataSets.append(int(folders[k]))

            self.nd.readNMRPipeSpcs(dataPath, dataSets, procDataName)
        # end readNMRPipeSpcs

    def readSpcs(self, dataPath, dataSets):
        zFill = 25
        if (dataPath[0] == 'interactive'):
            dataPath = [QFileDialog.getExistingDirectory()]

        if (len(dataPath) > 0):
            if (str(dataSets[0]) == 'all'):
                folders = []
                for r, d, f in os.walk(dataPath[0]):
                    for folder in d:
                        if (os.path.isfile(os.path.join(r, folder, 'fid'))):
                            if (folder != '99999'):
                                folders.append(folder.zfill(zFill))

                        if (os.path.isfile(os.path.join(r, folder, 'ser'))):
                            folders.append(folder.zfill(zFill))

                folders.sort()
                dataSets = []
                for k in range(len(folders)):
                    dataSets.append(int(folders[k]))

            if (str(dataSets[0]) == 'all1d'):
                folders = []
                for r, d, f in os.walk(dataPath[0]):
                    for folder in d:
                        if (os.path.isfile(os.path.join(r, folder, 'fid'))):
                            if (folder != '99999'):
                                folders.append(folder.zfill(zFill))

                folders.sort()
                dataSets = []
                for k in range(len(folders)):
                    dataSets.append(int(folders[k]))

            if (str(dataSets[0]) == 'all2d'):
                folders = []
                for r, d, f in os.walk(dataPath[0]):
                    for folder in d:
                        if (os.path.isfile(os.path.join(r, folder, 'ser'))):
                            folders.append(folder.zfill(zFill))

                folders.sort()
                dataSets = []
                for k in range(len(folders)):
                    dataSets.append(int(folders[k]))

            if len(dataPath) > 1:
                dp = []
                for d in dataPath:
                    if os.path.isfile(os.path.join(d, dataSets[0], 'fid')) or os.path.isfile(
                            os.path.join(d, dataSets[0], 'ser')):
                        dp.append(d)

                dataPath = dp

            else:
                ds = []
                for d in dataSets:
                    if os.path.isfile(os.path.join(dataPath[0], str(d), 'fid')) or os.path.isfile(
                            os.path.join(dataPath[0], str(d), 'ser')):
                        ds.append(d)

                dataSets = ds

            if len(dataPath) > 0 and len(dataSets) > 0:
                self.nd.readSpcs(dataPath, dataSets)

        # end readSpcs

    def reference1d(self, refShift=0.0):
        self.w.MplWidget.canvas.setFocus()
        self.showNMRSpectrum()
        xy = self.w.MplWidget.canvas.axes.figure.ginput(1)
        self.nd.nmrdat[self.nd.s][self.nd.e].refPoint[0] = self.nd.nmrdat[self.nd.s][self.nd.e].ppm2points(xy[0][0], 0)
        self.nd.nmrdat[self.nd.s][self.nd.e].refShift[0] = refShift
        self.nd.nmrdat[self.nd.s][self.nd.e].calcPPM()
        self.resetPlot()
        # end reference1d

    def reference2d(self, refShift=[0.0, 0.0]):
        self.w.MplWidget.canvas.setFocus()
        self.showNMRSpectrum()
        xy = self.w.MplWidget.canvas.axes.figure.ginput(1)
        self.nd.nmrdat[self.nd.s][self.nd.e].refPoint[0] = self.nd.nmrdat[self.nd.s][self.nd.e].ppm2points(xy[0][0], 0)
        self.nd.nmrdat[self.nd.s][self.nd.e].refShift[0] = refShift[0]
        self.nd.nmrdat[self.nd.s][self.nd.e].refPoint[1] = self.nd.nmrdat[self.nd.s][self.nd.e].ppm2points(xy[0][1], 1)
        self.nd.nmrdat[self.nd.s][self.nd.e].refShift[1] = refShift[1]
        self.nd.nmrdat[self.nd.s][self.nd.e].proc.refPoint[0] = self.nd.nmrdat[self.nd.s][self.nd.e].refPoint[0]*self.nd.nmrdat[self.nd.s][self.nd.e].proc.nPoints[0]/(len(self.nd.nmrdat[self.nd.s][self.nd.e].fid[0])*self.nd.nmrdat[self.nd.s][self.nd.e].proc.multFactor[0])
        self.nd.nmrdat[self.nd.s][self.nd.e].proc.refPoint[1] = self.nd.nmrdat[self.nd.s][self.nd.e].refPoint[1]*self.nd.nmrdat[self.nd.s][self.nd.e].proc.nPoints[1]/(len(self.nd.nmrdat[self.nd.s][self.nd.e].fid)*self.nd.nmrdat[self.nd.s][self.nd.e].proc.multFactor[1])
        self.nd.nmrdat[self.nd.s][self.nd.e].calcPPM()
        self.resetPlot()
        # end reference1d

    def removeLastColRow(self):
        nLines = len(self.w.MplWidget.canvas.axes.lines)
        if nLines > 0:
            self.w.MplWidget.canvas.axes.lines[nLines - 1].remove()
            self.phCorr.spcRow = self.phCorr.spcRow[:-1]
            self.phCorr.spcCol = self.phCorr.spcCol[:-1]
            self.phCorr.spcRowPts = self.phCorr.spcRowPts[:-1]
            self.phCorr.spcColPts = self.phCorr.spcColPts[:-1]
            self.plot2dColRow()
            self.showAcquisitionParameters()
            self.showNMRSpectrum()

        # end removeLastColRow

    def resetConfig(self):
        self.cf = nmrConfig.NmrConfig()
        self.cf.saveConfig()
        self.loadConfig()
        # end resetConfig

    def resetDataPreProcessing(self):
        self.nd.resetDataPreProcessing()
        self.plotSpcPreProc()
        self.verticalAutoScale()
        self.w.MplWidget.canvas.flush_events()
        self.w.MplWidget.canvas.draw()
        # end dataPreProcessing

    def resetHelp(self):
        if self.cf.mode == 'dark':
            fName = os.path.join(os.path.dirname(__file__), "nmr", "web", "introductionDark", "index.html")
        else:
            fName = os.path.join(os.path.dirname(__file__), "nmr", "web", "introduction", "index.html")

        url = "file:///" + fName.replace('\\', '/')
        self.w.helpView.setUrl(url)
        self.w.nmrSpectrum.setCurrentIndex(12)
        # end resetHelp

    def resetPlot(self):
        zoomChecked = self.w.keepZoom.isChecked()
        self.w.keepZoom.setChecked(False)
        self.plotSpc()
        if (zoomChecked == True):
            self.w.keepZoom.setChecked(True)

        # end resetPlot

    def saveButton(self):
        pfName = QFileDialog.getSaveFileName(None, "Save MetaboLabPy DataSet", "", "*.mlpy", "*.mlpy")
        if (os.path.isfile(pfName[0])):
            os.remove(pfName[0])

        if (os.path.isdir(pfName[0])):
            shutil.rmtree(pfName[0])

        self.nd.script = self.w.script.toHtml()
        self.nd.console = self.w.console.toHtml()
        self.nd.save(pfName[0])
        # end saveButton

    def saveConfig(self):
        self.cf.autoPlot = self.w.autoPlot.isChecked()
        self.cf.keepZoom = self.w.keepZoom.isChecked()
        self.cf.fontSize = self.w.fontSize.value()
        self.cf.phaseReferenceColour = self.nd.nmrdat[0][0].display.phRefCol
        self.cf.posCol10 = self.stdPosCol1[0]
        self.cf.posCol11 = self.stdPosCol1[1]
        self.cf.posCol12 = self.stdPosCol1[2]
        self.cf.negCol10 = self.stdNegCol1[0]
        self.cf.negCol11 = self.stdNegCol1[1]
        self.cf.negCol12 = self.stdNegCol1[2]
        self.cf.posCol20 = self.stdPosCol2[0]
        self.cf.posCol21 = self.stdPosCol2[1]
        self.cf.posCol22 = self.stdPosCol2[2]
        self.cf.negCol20 = self.stdNegCol2[0]
        self.cf.negCol21 = self.stdNegCol2[1]
        self.cf.negCol22 = self.stdNegCol2[2]
        self.cf.saveConfig()
        # end saveConfig

    def saveMat(self):
        scipy.io.savemat('/Users/ludwigc/metabolabpy.mat',
                         {'spc': self.nd.nmrdat[0][0].spc, 'fid': self.nd.nmrdat[0][0].fid})
        # end saveMat

    def saveScript(self, fName=""):
        if (fName == False):
            fName = ""

        if (len(fName) == 0):
            fName = QFileDialog.getSaveFileName(None, 'Save Script File', '', 'Python files (*.py)')
            fName = fName[0]

        if (len(fName) > 0):
            scriptText = self.w.script.toPlainText()
            f = open(fName, 'w')
            f.write(scriptText)

        # end openScript

    def scale2DSpectrumUp(self):
        self.nd.nmrdat[self.nd.s][self.nd.e].display.minLevel /= 1.1
        self.nd.nmrdat[self.nd.s][self.nd.e].display.maxLevel /= 1.1
        self.setDispPars()
        self.plotSpc()
        # end scale2DSpectrumUp

    def scale2DSpectrumDown(self):
        self.nd.nmrdat[self.nd.s][self.nd.e].display.minLevel *= 1.1
        self.nd.nmrdat[self.nd.s][self.nd.e].display.maxLevel *= 1.1
        self.setDispPars()
        self.plotSpc()
        # end scale2DSpectrumDown

    def scaleAll2DSpectraUp(self):
        self.nd.nmrdat[self.nd.s][self.nd.e].display.minLevel /= 1.1
        self.nd.nmrdat[self.nd.s][self.nd.e].display.maxLevel /= 1.1
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            self.nd.nmrdat[self.nd.s][k].display.minLevel = self.nd.nmrdat[self.nd.s][self.nd.e].display.minLevel
            self.nd.nmrdat[self.nd.s][k].display.maxLevel = self.nd.nmrdat[self.nd.s][self.nd.e].display.maxLevel

        self.setDispPars()
        self.plotSpc()
        # end scaleAll2DSpectraUp

    def scaleAll2DSpectraDown(self):
        self.nd.nmrdat[self.nd.s][self.nd.e].display.minLevel *= 1.1
        self.nd.nmrdat[self.nd.s][self.nd.e].display.maxLevel *= 1.1
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            self.nd.nmrdat[self.nd.s][k].display.minLevel = self.nd.nmrdat[self.nd.s][self.nd.e].display.minLevel
            self.nd.nmrdat[self.nd.s][k].display.maxLevel = self.nd.nmrdat[self.nd.s][self.nd.e].display.maxLevel

        self.setDispPars()
        self.plotSpc()
        # end scaleAll2DSpectraDown

    def scriptEditor(self):
        self.w.nmrSpectrum.setCurrentIndex(10)
        # end scriptEditor

    def selectAddCompressPreProc(self):
        xy = self.w.MplWidget.canvas.axes.figure.ginput(2)
        t = np.round(1e4 * np.array([xy[0][0], xy[1][0]])) / 1e4
        self.nd.pp.compressStart = np.append(self.nd.pp.compressStart, min(t))
        self.nd.pp.compressEnd = np.append(self.nd.pp.compressEnd, max(t))
        self.fillPreProcessingNumbers()
        self.w.compressBucketsTW.setFocus()
        self.setPlotPreProc()
        self.w.compressBucketsTW.setFocus()
        self.plotSpcPreProc()
        self.setCompressPreProc()
        # end selectAddExcludePreProc

    def selectAddExcludePreProc(self):
        xy = self.w.MplWidget.canvas.axes.figure.ginput(2)
        t = np.round(1e4 * np.array([xy[0][0], xy[1][0]])) / 1e4
        self.nd.pp.excludeStart = np.append(self.nd.pp.excludeStart, min(t))
        self.nd.pp.excludeEnd = np.append(self.nd.pp.excludeEnd, max(t))
        self.fillPreProcessingNumbers()
        self.w.excludeRegionTW.setFocus()
        self.setPlotPreProc()
        self.w.excludeRegionTW.setFocus()
        self.plotSpcPreProc()
        self.setExcludePreProc()
        # end selectAddExcludePreProc

    def selectAddSegAlignPreProc(self):
        xy = self.w.MplWidget.canvas.axes.figure.ginput(2)
        t = np.round(1e4 * np.array([xy[0][0], xy[1][0]])) / 1e4
        self.nd.pp.segStart = np.append(self.nd.pp.segStart, min(t))
        self.nd.pp.segEnd = np.append(self.nd.pp.segEnd, max(t))
        self.fillPreProcessingNumbers()
        self.w.segAlignTW.setFocus()
        self.setPlotPreProc()
        self.w.segAlignTW.setFocus()
        self.plotSpcPreProc()
        self.setSegAlignPreProc()
        # end selectAddExcludePreProc

    def selectAllPreProc(self):
        nSpc = len(self.nd.pp.classSelect)
        self.nd.pp.plotSelect = np.arange(nSpc)
        self.fillPreProcessingNumbers()
        self.setPlotPreProc()
        self.plotSpcPreProc()
        self.w.selectClassTW.setFocus()
        # end selectAllPreProc

    def selectClassPreProc(self):
        cls = self.w.selectClassLE.text()
        cls2 = self.nd.pp.classSelect
        sel = np.array([])
        for k in range(len(cls2)):
            if (cls2[k] == cls):
                sel = np.append(sel, k)

        if (len(sel) == 0):
            sel = np.arange(len(cls2))

        self.nd.pp.plotSelect = sel
        self.fillPreProcessingNumbers()
        self.setPlotPreProc()
        self.plotSpcPreProc()
        self.w.selectClassTW.setFocus()
        # end selectClassPreProc

    def selectClearCompressPreProc(self):
        self.nd.pp.preProcFill = True
        for k in range(len(self.nd.pp.compressStart)):
            self.w.compressBucketsTW.item(k, 0).setText("")
            self.w.compressBucketsTW.setFocus()
            self.w.compressBucketsTW.item(k, 1).setText("")
            self.w.compressBucketsTW.setFocus()

        self.nd.pp.preProcFill = False
        self.nd.pp.compressStart = np.array([])
        self.nd.pp.compressEnd = np.array([])
        self.w.compressBucketsTW.setFocus()
        self.fillPreProcessingNumbers()
        self.w.compressBucketsTW.setFocus()
        self.setPlotPreProc()
        self.w.compressBucketsTW.setFocus()
        self.plotSpcPreProc()
        self.setCompressPreProc()
        self.w.MplWidget.canvas.flush_events()
        self.w.MplWidget.canvas.draw()
        # end selectClearExcludePreProc

    def selectClearExcludePreProc(self):
        self.nd.pp.preProcFill = True
        for k in range(len(self.nd.pp.excludeStart)):
            self.w.excludeRegionTW.item(k, 0).setText("")
            self.w.excludeRegionTW.setFocus()
            self.w.excludeRegionTW.item(k, 1).setText("")
            self.w.excludeRegionTW.setFocus()

        self.nd.pp.preProcFill = False
        self.nd.pp.excludeStart = np.array([])
        self.nd.pp.excludeEnd = np.array([])
        self.w.excludeRegionTW.setFocus()
        self.fillPreProcessingNumbers()
        self.w.excludeRegionTW.setFocus()
        self.setPlotPreProc()
        self.w.excludeRegionTW.setFocus()
        self.plotSpcPreProc()
        self.setExcludePreProc()
        self.w.MplWidget.canvas.flush_events()
        self.w.MplWidget.canvas.draw()
        # end selectClearExcludePreProc

    def selectClearSegAlignPreProc(self):
        self.nd.pp.preProcFill = True
        for k in range(len(self.nd.pp.segStart)):
            self.w.segAlignTW.item(k, 0).setText("")
            self.w.segAlignTW.setFocus()
            self.w.segAlignTW.item(k, 1).setText("")
            self.w.segAlignTW.setFocus()

        self.nd.pp.preProcFill = False
        self.nd.pp.segStart = np.array([])
        self.nd.pp.segEnd = np.array([])
        self.w.segAlignTW.setFocus()
        self.fillPreProcessingNumbers()
        self.w.segAlignTW.setFocus()
        self.setPlotPreProc()
        self.w.segAlignTW.setFocus()
        self.plotSpcPreProc()
        self.setSegAlignPreProc()
        self.w.MplWidget.canvas.flush_events()
        self.w.MplWidget.canvas.draw()
        # end selectClearExcludePreProc

    def selectEvenPreProc(self):
        nSpc = len(self.nd.pp.classSelect)
        self.nd.pp.plotSelect = np.arange(nSpc)
        self.nd.pp.plotSelect = self.nd.pp.plotSelect[1::2]
        self.fillPreProcessingNumbers()
        self.setPlotPreProc()
        self.plotSpcPreProc()
        self.w.selectClassTW.setFocus()
        # end selectEvenPreProc

    def selectOddPreProc(self):
        nSpc = len(self.nd.pp.classSelect)
        self.nd.pp.plotSelect = np.arange(nSpc)
        self.nd.pp.plotSelect = self.nd.pp.plotSelect[0::2]
        self.fillPreProcessingNumbers()
        self.setPlotPreProc()
        self.plotSpcPreProc()
        self.w.selectClassTW.setFocus()
        # end selectOddPreProc

    def selectPlotAll(self):
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            self.nd.nmrdat[self.nd.s][k].display.displaySpc = True

        # self.plotSpc()
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.changeDataSetExp()
        self.plotSpc()
        return "selectPlotAll"
        # end selectPlotAll

    def selectPlotClear(self):
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            self.nd.nmrdat[self.nd.s][k].display.displaySpc = False

        # self.plotSpc()
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.changeDataSetExp()
        self.plotSpc()
        return "selectPlotClear"
        # end selectPlotClear

    def selectPlotList(self, plotSelect):
        plotSelect = np.array(plotSelect)
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            self.nd.nmrdat[self.nd.s][k].display.displaySpc = False

        plotSelect -= 1
        for k in range(len(plotSelect)):
            if ((plotSelect[k] > -1) and (plotSelect[k] < len(self.nd.nmrdat[self.nd.s]))):
                self.nd.nmrdat[self.nd.s][plotSelect[k]].display.displaySpc = True

        # self.plotSpc()
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.changeDataSetExp()
        return "selectPlotList"
        # end selectPlotList

    def setAcqPars(self):
        s = self.nd.s
        e = self.nd.e
        a = self.nd.nmrdat[s][e].acq
        acqStr = "originalDataset      " + self.nd.nmrdat[s][e].origDataSet + "\n"
        acqStr += "___________________________________________________________________________________________________\n"
        acqStr += "\n"
        acqStr += "metaInfo             "
        for k in range(len(a.title)):
            acqStr += a.title[k] + " "

        acqStr += "\n                    "
        acqStr += " Origin\t" + a.origin + "\n                    "
        acqStr += " Owner\t" + a.owner + "\n"
        acqStr += "___________________________________________________________________________________________________\n"
        acqStr += "\n"
        acqStr += "probe                          " + a.probe + "\n"
        pp = a.pulProgName
        pp = pp[1:]
        pp = pp[:len(pp) - 1]
        acqStr += "pulseProgram                   " + pp + "\n\n"
        acqStr += "sw                   [ppm]    " + "% 9.2f"%a.sw[0]           + "        |    % 9.2f"%a.sw[1]           + "        |    % 9.2f\n"%a.sw[2]
        acqStr += "sw_h                 [Hz]     " + "% 9.2f"%a.sw_h[0]         + "        |    % 9.2f"%a.sw_h[1]         + "        |    % 9.2f\n"%a.sw_h[2]
        acqStr += "bf1/2/3              [MHz]    " + "% 9.2f"%a.bf1             + "        |    % 9.2f"%a.bf2             + "        |    % 9.2f\n"%a.bf3
        acqStr += "sfo1/2/3             [MHz]    " + "% 9.2f"%a.sfo1            + "        |    % 9.2f"%a.sfo2            + "        |    % 9.2f\n"%a.sfo3
        acqStr += "o1/2/3               [Hz]     " + "% 9.2f"%a.o1              + "        |    % 9.2f"%a.o2              + "        |    % 9.2f\n"%a.o3
        acqStr += "nPoints                       " + "% 6d"%a.nDataPoints[0] + "           |    % 6d"%a.nDataPoints[1] + "           |    % 6d\n"%a.nDataPoints[2]
        acqStr += "transients                    " + "% 6d\n"%a.transients
        acqStr += "steadyStateScans              " + "% 6d\n\n"%a.steadyStateScans
        acqStr += "groupDelay           [us]     " + "% 9.2f\n"%a.groupDelay
        acqStr += "decim                         " + "% 6d\n"%a.decim
        acqStr += "dspfvs                        " + "% 6d\n"%a.dspfvs
        acqStr += "temperature          [K]      " + "% 9.2f\n"%a.temperature
        self.w.acqPars.setText(acqStr)
        # end setAcqPars

    def setAvoidNegValues(self):
        if (self.nd.pp.preProcFill == False):
            if (self.w.avoidNegValues.isChecked() == True):
                self.nd.pp.avoidNegativeValues = True
            else:
                self.nd.pp.avoidNegativeValues = False

        # end setAvoidNegValues

    def setBucketPPMPreProc(self):
        try:
            bucketPPM = float(self.w.bucketPpmLE.text())
        except:
            bucketPPM = self.nd.pp.bucketPPM

        ppmPerPoint = abs(self.nd.nmrdat[self.nd.s][0].ppm1[0] - self.nd.nmrdat[self.nd.s][0].ppm1[1])
        bucketPoints = round(bucketPPM / ppmPerPoint)
        bucketPPM = np.round(1e4 * bucketPoints * ppmPerPoint) / 1e4
        self.w.bucketPpmLE.setText(str(bucketPPM))
        self.w.bucketDataPointsLE.setText(str(int(bucketPoints)))
        self.nd.pp.bucketPoints = bucketPoints
        self.nd.pp.bucketPPM = bucketPPM
        # end setBucketPPMPreProc

    def setBucketPointsPreProc(self):
        try:
            bucketPoints = float(self.w.bucketDataPointsLE.text())
        except:
            bucketPoints = self.nd.pp.bucketPoints

        ppmPerPoint = abs(self.nd.nmrdat[self.nd.s][0].ppm1[0] - self.nd.nmrdat[self.nd.s][0].ppm1[1])
        bucketPoints = round(bucketPoints)
        bucketPPM = np.round(1e4 * bucketPoints * ppmPerPoint) / 1e4
        self.w.bucketPpmLE.setText(str(bucketPPM))
        self.w.bucketDataPointsLE.setText(str(int(bucketPoints)))
        self.nd.pp.bucketPoints = bucketPoints
        self.nd.pp.bucketPPM = bucketPPM
        # end setBucketPointsPreProc

    def setBucketSpectra(self):
        if (self.nd.pp.preProcFill == False):
            if (self.w.bucketSpectra.isChecked() == True):
                self.nd.pp.flagBucketSpectra = True
                self.w.preProcessingSelect.setCurrentIndex(4)
            else:
                self.nd.pp.flagBucketSpectra = False

        # end setBucketSpectra

    def setChangePreProc(self):
        if (self.nd.pp.preProcFill == False):
            cls = np.array([])
            for k in range(len(self.nd.pp.classSelect)):
                cls = np.append(cls, self.w.selectClassTW.item(k, 1).text())

            self.nd.pp.classSelect = cls

        # end setChangePreProc

    def setCompressBuckets(self):
        if (self.nd.pp.preProcFill == False):
            if (self.w.compressBuckets.isChecked() == True):
                self.nd.pp.flagCompressBuckets = True
                self.w.preProcessingSelect.setCurrentIndex(5)
            else:
                self.nd.pp.flagCompressBuckets = False

        # end setExcludeRegion

    def setCompressPreProc(self):
        if (self.nd.pp.preProcFill == False):
            nRows = self.w.compressBucketsTW.rowCount()
            coStart = np.array([])
            coEnd = np.array([])
            tStart = np.array([])
            tEnd = np.array([])
            for k in range(nRows):
                # tStart = np.array([])
                # tEnd   = np.array([])
                try:
                    tStart = np.append(tStart, float(self.w.compressBucketsTW.item(k, 0).text()))
                    # self.w.compressBucketsTW.item(k,0).clearContents()
                except:
                    tStart = np.append(tStart, -10000.0)

                try:
                    tEnd = np.append(tEnd, float(self.w.compressBucketsTW.item(k, 1).text()))
                    # self.w.compressBucketsTW.item(k,1).clearContents()
                except:
                    tEnd = np.append(tEnd, -10000.0)

            # self.w.compressBucketsTW.clearContents()
            self.w.compressBucketsTW.setRowCount(0)
            self.w.compressBucketsTW.setRowCount(nRows)
            self.nd.pp.preProcFill = True
            for k in np.arange(len(tStart) - 1, -1, -1):  # range(len(tStart)):
                compNumber1 = QTableWidgetItem(2 * k)
                compNumber1.setTextAlignment(QtCore.Qt.AlignHCenter)
                self.w.compressBucketsTW.setItem(k, 0, compNumber1)
                compNumber2 = QTableWidgetItem(2 * k + 1)
                compNumber2.setTextAlignment(QtCore.Qt.AlignHCenter)
                self.w.compressBucketsTW.setItem(k, 1, compNumber2)
                if ((tStart[k] > -10000.0) & (tEnd[k] > -10000.0)):
                    tMin = min(tStart[k], tEnd[k])
                    tEnd[k] = max(tStart[k], tEnd[k])
                    tStart[k] = tMin
                    coStart = np.append(coStart, tStart[k])
                    coEnd = np.append(coEnd, tEnd[k])
                    tStart = np.delete(tStart, k)
                    tEnd = np.delete(tEnd, k)

                if (tStart[k] > -10000.0):
                    self.w.compressBucketsTW.item(k, 0).setText(str(tStart[k]))
                    self.w.compressBucketsTW.setFocus()
                else:
                    self.w.compressBucketsTW.item(k, 0).setText("")
                    self.w.compressBucketsTW.setFocus()

                if (tEnd[k] > -10000.0):
                    self.w.compressBucketsTW.item(k, 1).setText(str(tEnd[k]))
                    self.w.compressBucketsTW.setFocus()
                else:
                    self.w.compressBucketsTW.item(k, 1).setText("")
                    self.w.compressBucketsTW.setFocus()

            self.nd.pp.preProcFill = False
            sortIdx = np.argsort(coStart)
            self.nd.pp.compressStart = coStart[sortIdx]
            self.nd.pp.compressEnd = coEnd[sortIdx]
            self.plotSpcPreProc()

        # end setCompressPreProc

    def setDarkMode(self):
        self.cf.readConfig()
        self.cf.mode = 'dark'
        self.cf.saveConfig()
        # end saveConfig

    def setDispPars(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].display
        self.w.posColR.setText(str(d.posColRGB[0]))
        self.w.posColG.setText(str(d.posColRGB[1]))
        self.w.posColB.setText(str(d.posColRGB[2]))
        self.w.negColR.setText(str(d.negColRGB[0]))
        self.w.negColG.setText(str(d.negColRGB[1]))
        self.w.negColB.setText(str(d.negColRGB[2]))
        self.w.nLevels.setText(str(d.nLevels))
        self.w.minLevel.setText(str(d.minLevel))
        self.w.maxLevel.setText(str(d.maxLevel))
        self.w.spcOffset.setText(str(d.spcOffset))
        self.w.spcScale.setText(str(d.spcScale))
        self.w.xLabel.setText(d.xLabel)
        self.w.yLabel.setText(d.yLabel)
        self.w.spcLabel.setText(d.spcLabel)
        self.w.posCol.setCurrentIndex(d.colours.get(d.posCol))
        self.w.negCol.setCurrentIndex(d.colours.get(d.negCol))
        self.w.axisType1.setCurrentIndex(d.axes.get(d.axisType1))
        self.w.axisType2.setCurrentIndex(d.axes.get(d.axisType2))
        self.w.displaySpc.setCurrentIndex(d.falseTrue.get(d.displaySpc))
        self.w.phRefColour.setCurrentIndex(d.colours2.get(d.phRefCol))
        self.w.phRefDS.setValue(d.phRefDS)
        self.w.phRefExp.setValue(d.phRefExp)
        # end setDispPars

    def setExcludePreProc(self):
        if (self.nd.pp.preProcFill == False):
            nRows = self.w.excludeRegionTW.rowCount()
            exStart = np.array([])
            exEnd = np.array([])
            tStart = np.array([])
            tEnd = np.array([])
            for k in range(nRows):
                # tStart = np.array([])
                # tEnd   = np.array([])
                try:
                    tStart = np.append(tStart, float(self.w.excludeRegionTW.item(k, 0).text()))
                    # self.w.excludeRegionTW.item(k,0).clearContents()
                except:
                    tStart = np.append(tStart, -10000.0)

                try:
                    tEnd = np.append(tEnd, float(self.w.excludeRegionTW.item(k, 1).text()))
                    # self.w.excludeRegionTW.item(k,1).clearContents()
                except:
                    tEnd = np.append(tEnd, -10000.0)

            # self.w.excludeRegionTW.clearContents()
            self.w.excludeRegionTW.setRowCount(0)
            self.w.excludeRegionTW.setRowCount(nRows)
            self.nd.pp.preProcFill = True
            for k in np.arange(len(tStart) - 1, -1, -1):  # range(len(tStart)):
                exclNumber1 = QTableWidgetItem(2 * k)
                exclNumber1.setTextAlignment(QtCore.Qt.AlignHCenter)
                self.w.excludeRegionTW.setItem(k, 0, exclNumber1)
                exclNumber2 = QTableWidgetItem(2 * k + 1)
                exclNumber2.setTextAlignment(QtCore.Qt.AlignHCenter)
                self.w.excludeRegionTW.setItem(k, 1, exclNumber2)
                if ((tStart[k] > -10000.0) & (tEnd[k] > -10000.0)):
                    tMin = min(tStart[k], tEnd[k])
                    tEnd[k] = max(tStart[k], tEnd[k])
                    tStart[k] = tMin
                    exStart = np.append(exStart, tStart[k])
                    exEnd = np.append(exEnd, tEnd[k])
                    tStart = np.delete(tStart, k)
                    tEnd = np.delete(tEnd, k)

                if (tStart[k] > -10000.0):
                    self.w.excludeRegionTW.item(k, 0).setText(str(tStart[k]))
                    self.w.excludeRegionTW.setFocus()
                else:
                    self.w.excludeRegionTW.item(k, 0).setText("")
                    self.w.excludeRegionTW.setFocus()

                if (tEnd[k] > -10000.0):
                    self.w.excludeRegionTW.item(k, 1).setText(str(tEnd[k]))
                    self.w.excludeRegionTW.setFocus()
                else:
                    self.w.excludeRegionTW.item(k, 1).setText("")
                    self.w.excludeRegionTW.setFocus()

            self.nd.pp.preProcFill = False
            sortIdx = np.argsort(exStart)
            self.nd.pp.excludeStart = exStart[sortIdx]
            self.nd.pp.excludeEnd = exEnd[sortIdx]
            self.plotSpcPreProc()

        # end setExcludePreProc

    def setExcludeRegion(self):
        if (self.nd.pp.preProcFill == False):
            if (self.w.excludeRegion.isChecked() == True):
                self.nd.pp.flagExcludeRegion = True
                self.w.preProcessingSelect.setCurrentIndex(1)
            else:
                self.nd.pp.flagExcludeRegion = False

        # end setExcludeRegion

    def setExportCharacter(self):
        tt = self.w.exportCharacter.text()
        if (len(tt) > 0):
            self.nd.pp.exportCharacter = tt[0]
            self.w.exportCharacter.setText(tt[0])

        # end setExportCharacter

    def setExportDelimiterTab(self):
        self.nd.pp.exportDelimiterTab = self.w.exportDelimiterTab.isChecked()
        # end setExportDelimiterTab

    def setExportFileName(self):
        if self.nd.pp.exportMethod == 0:
            self.nd.pp.exportExcel = self.w.exportFileName.text()

        if self.nd.pp.exportMethod == 1:
            self.nd.pp.exportFileName = self.w.exportFileName.text()

        if self.nd.pp.exportMethod == 2:
            self.nd.pp.exportMetaboAnalyst = self.w.exportFileName.text()

        if self.nd.pp.exportMethod == 3:
            self.nd.pp.exportrDolphin = self.w.exportFileName.text()

        if self.nd.pp.exportMethod == 4:
            self.nd.pp.exportBatman = self.w.exportFileName.text()

        if self.nd.pp.exportMethod == 5:
            self.nd.pp.exportBruker = self.w.exportFileName.text()

        # end setExportFileName

    def setExportMethod(self):
        if self.nd.pp.exportMethod == 0:
            self.w.delimiterLabel.setHidden(True)
            self.w.exportDelimiterTab.setHidden(True)
            self.w.exportDelimiterCharacter.setHidden(True)
            self.w.exportCharacter.setHidden(True)
            self.w.samplesInRowsLabel.setHidden(False)
            self.w.samplesInComboBox.setHidden(False)
            self.w.exportPath.setText(self.nd.pp.exportExcelPath)
            self.w.exportFileName.setText(self.nd.pp.exportExcel)

        if self.nd.pp.exportMethod == 1:
            self.w.delimiterLabel.setHidden(False)
            self.w.exportDelimiterTab.setHidden(False)
            self.w.exportDelimiterCharacter.setHidden(False)
            self.w.exportCharacter.setHidden(False)
            self.w.samplesInRowsLabel.setHidden(False)
            self.w.samplesInComboBox.setHidden(False)
            self.w.exportPath.setText(self.nd.pp.exportPathName)
            self.w.exportFileName.setText(self.nd.pp.exportFileName)

        if self.nd.pp.exportMethod == 2:
            self.w.delimiterLabel.setHidden(True)
            self.w.exportDelimiterTab.setHidden(True)
            self.w.exportDelimiterCharacter.setHidden(True)
            self.w.exportCharacter.setHidden(True)
            self.w.samplesInRowsLabel.setHidden(True)
            self.w.samplesInComboBox.setHidden(True)
            self.w.exportPath.setText(self.nd.pp.exportMetaboAnalystPath)
            self.w.exportFileName.setText(self.nd.pp.exportMetaboAnalyst)

        if self.nd.pp.exportMethod == 3:
            self.w.delimiterLabel.setHidden(True)
            self.w.exportDelimiterTab.setHidden(True)
            self.w.exportDelimiterCharacter.setHidden(True)
            self.w.exportCharacter.setHidden(True)
            self.w.samplesInRowsLabel.setHidden(True)
            self.w.samplesInComboBox.setHidden(True)
            self.w.exportPath.setText(self.nd.pp.exportrDolphinPath)
            self.w.exportFileName.setText(self.nd.pp.exportrDolphin)

        if self.nd.pp.exportMethod == 4:
            self.w.delimiterLabel.setHidden(True)
            self.w.exportDelimiterTab.setHidden(True)
            self.w.exportDelimiterCharacter.setHidden(True)
            self.w.exportCharacter.setHidden(True)
            self.w.samplesInRowsLabel.setHidden(True)
            self.w.samplesInComboBox.setHidden(True)
            self.w.exportPath.setText(self.nd.pp.exportBatmanPath)
            self.w.exportFileName.setText(self.nd.pp.exportBatman)

        if self.nd.pp.exportMethod == 5:
            self.w.delimiterLabel.setHidden(True)
            self.w.exportDelimiterTab.setHidden(True)
            self.w.exportDelimiterCharacter.setHidden(True)
            self.w.exportCharacter.setHidden(True)
            self.w.samplesInRowsLabel.setHidden(True)
            self.w.samplesInComboBox.setHidden(True)
            self.w.exportPath.setText(self.nd.pp.exportBrukerPath)
            self.w.exportFileName.setText(self.nd.pp.exportBruker)

        # end setExportMethod

    def setExportMethodOptions(self):
        self.nd.pp.exportMethod = self.w.exportMethod.currentIndex()
        self.setExportMethod()
        # end setExportMethodOptions

    def setExportPath(self):
        if self.nd.pp.exportMethod == 0:
            self.nd.pp.exportExcelPath = self.w.exportPath.text()

        if self.nd.pp.exportMethod == 1:
            self.nd.pp.exportPathName = self.w.exportPath.text()

        if self.nd.pp.exportMethod == 2:
            self.nd.pp.exportMetaboAnalystPath = self.w.exportPath.text()

        if self.nd.pp.exportMethod == 3:
            self.nd.pp.exportrDolphinPath = self.w.exportPath.text()

        if self.nd.pp.exportMethod == 4:
            self.nd.pp.exportBatmanPath = self.w.exportPath.text()

        if self.nd.pp.exportMethod == 5:
            self.nd.pp.exportBrukerPath = self.w.exportPath.text()

        # end setExportPath

    def setExportDataSet(self):
        if (self.nd.pp.preProcFill == False):
            if (self.w.exportDataSet.isChecked() == True):
                self.nd.pp.flagExportDataSet = True
                self.w.preProcessingSelect.setCurrentIndex(8)
            else:
                self.nd.pp.flagExportDataSet = False

        # end setExportDataSet

    def setExportTable(self):
        pName = QFileDialog.getExistingDirectory()
        # pName = pName[0]
        if (len(pName) > 0):
            if self.nd.pp.exportMethod == 0:
                self.w.exportPath.setText(pName)
                self.nd.pp.exportExcelPath = pName

            if self.nd.pp.exportMethod == 1:
                self.w.exportPath.setText(pName)
                self.nd.pp.exportPathName = pName

            if self.nd.pp.exportMethod == 2:
                self.w.exportPath.setText(pName)
                self.nd.pp.exportMetaboAnalystPath = pName

            if self.nd.pp.exportMethod == 3:
                self.w.exportPath.setText(pName)
                self.nd.pp.exportrDolphinPath = pName

            if self.nd.pp.exportMethod == 4:
                self.w.exportPath.setText(pName)
                self.nd.pp.exportBatmanPath = pName

            if self.nd.pp.exportMethod == 5:
                self.w.exportPath.setText(pName)
                self.nd.pp.exportBrukerPath = pName

        # end setExportTable

    def setFontSize(self):
        fontSize = self.w.fontSize.value()
        f = self.w.acqPars.font()
        f.setPointSize(fontSize)
        self.w.acqPars.setFont(f)
        self.w.titleFile.setFont(f)
        cursor = self.w.script.textCursor()
        self.w.script.selectAll()
        self.w.script.setFontPointSize(fontSize)
        self.w.script.setTextCursor(cursor)
        self.w.script.setCurrentFont(f)
        # self.w.script.setFont(f)
        cursor = self.w.console.textCursor()
        self.w.console.selectAll()
        self.w.console.setFontPointSize(fontSize)
        self.w.console.setTextCursor(cursor)
        self.w.console.setCurrentFont(f)
        # self.w.console.setFont(f)
        self.w.pulseProgram.setFont(f)
        self.w.cmdLine.setFont(f)
        self.w.setStyleSheet("font-size: " + str(fontSize) + "pt")
        # end setFontSize

    def setHelp(self):
        url = []
        idx = self.w.helpComboBox.currentIndex()
        if self.cf.mode == 'dark':
            fName = os.path.join(os.path.dirname(__file__), "nmr", "web", "introductionDark", "index.html")
        else:
            fName = os.path.join(os.path.dirname(__file__), "nmr", "web", "introduction", "index.html")

        url.append("file:///" + fName.replace('\\', '/'))
        url.append("https://www.hmdb.ca")
        url.append("https://www.smpdb.ca")
        url.append("https://bmrb.io/metabolomics/")
        url.append("https://www.genome.jp/kegg/pathway.html#metabolism")
        url.append("https://nmrshiftdb.nmr.uni-koeln.de")
        url.append("https://sdbs.db.aist.go.jp/sdbs/cgi-bin/cre_index.cgi")
        url.append("http://dmar.riken.jp/spincouple/")
        self.w.helpView.setUrl(url[idx])
        # end setHelp

    def setJres(self):
        if (self.nd.nmrdat[self.nd.s][self.nd.e].acq.fnMode == 1):
            self.nd.nmrdat[self.nd.s][self.nd.e].display.yLabel = '1H'
            self.nd.nmrdat[self.nd.s][self.nd.e].display.axisType2 = 'Hz'
            self.nd.nmrdat[self.nd.s][self.nd.e].proc.windowType = np.array([5, 3, 0])
            self.nd.nmrdat[self.nd.s][self.nd.e].proc.lb[0] = 0.5

        # end setJres


    def setLightMode(self):
        self.cf.readConfig()
        self.cf.mode = 'light'
        self.cf.saveConfig()
        # end saveConfig

    def setNoiseFiltering(self):
        if (self.nd.pp.preProcFill == False):
            if (self.w.noiseFiltering.isChecked() == True):
                self.nd.pp.flagNoiseFiltering = True
                self.w.preProcessingSelect.setCurrentIndex(3)
            else:
                self.nd.pp.flagNoiseFiltering = False

        # end setNoiseFiltering

    def setnoiseRegPreProc(self):
        try:
            th = float(self.w.noiseThresholdLE.text())
        except:
            th = self.nd.pp.noiseThreshold

        try:
            ns = float(self.w.noiseRegionStartLE.text())
        except:
            ns = self.nd.pp.noiseStart

        try:
            ne = float(self.w.noiseRegionEndLE.text())
        except:
            ne = self.nd.pp.noiseEnd

        try:
            lw = float(self.w.thLineWidthLE.text())
        except:
            lw = self.nd.pp.thLineWidth

        tm = min(ns, ne)
        ne = max(ns, ne)
        ns = tm
        self.nd.pp.noiseThreshold = th
        self.nd.pp.noiseStart = ns
        self.nd.pp.noiseEnd = ne
        self.nd.pp.thLineWidth = lw
        self.w.noiseThresholdLE.setText(str(th))
        self.w.noiseRegionStartLE.setText(str(ns))
        self.w.noiseRegionEndLE.setText(str(ne))
        self.w.thLineWidthLE.setText(str(lw))
        self.plotSpcPreProc()
        # end setnoiseRegPreProc

    def setPhRefExp(self, phRefExp, phRefDS=1):
        self.w.phRefDS.setValue(phRefDS)
        self.w.phRefExp.setValue(phRefExp)
        # end setPhRefExp

    def setPlotPreProc(self):
        if (self.nd.pp.preProcFill == False):
            sel = np.array([])
            sel = self.w.selectClassTW.selectedIndexes()
            sel2 = np.array([])
            for k in range(len(sel)):
                if (sel[k].column() == 0):
                    sel2 = np.append(sel2, sel[k].row())

            self.nd.pp.plotSelect = sel2
            self.plotSpcPreProc()

        # end setPlotPreProc

    def setPqnTsaScaling(self):
        if self.w.pqnButton.isChecked() is True:
            self.nd.pp.scalePQN = True
        else:
            self.nd.pp.scalePQN = False

        self.w.preserveOverallScale.setDisabled(self.nd.pp.scalePQN)

        # end setPqnTsaScaling

    def setPreProcessing(self):
        if (self.w.preprocessing.isChecked() == True):
            if len(self.nd.nmrdat[self.nd.s]) != len(self.nd.pp.classSelect):
                self.nd.preProcInit()

            self.showPreProcessing()
            self.fillPreProcessingNumbers()
            self.nd.noiseFilteringInit()
        else:
            self.hidePreProcessing()

        # end setPreProcessing

    def setStandardColours(self):
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            if self.cf.mode == 'dark':
                self.nd.nmrdat[self.nd.s][k].display.posColRGB = self.stdPosCol2
                self.nd.nmrdat[self.nd.s][k].display.negColRGB = self.stdNegCol2
            else:
                self.nd.nmrdat[self.nd.s][k].display.posColRGB = self.stdPosCol1
                self.nd.nmrdat[self.nd.s][k].display.negColRGB = self.stdNegCol1


        self.plotSpc()

    def setHsqcAnalysis(self):
        if (self.w.hsqcAnalysis.isChecked() == True):
            self.w.multipletAnalysis.setVisible(True)
            self.w.isotopomerAnalysis.setVisible(True)
            self.w.nmrSpectrum.setTabEnabled(1, True)
            self.w.nmrSpectrum.setTabEnabled(2, True)
            self.w.nmrSpectrum.setStyleSheet("QTabBar::tab::disabled {width: 0; height: 0; margin: 0; padding: 0; border: none;} ")
            self.w.nmrSpectrum.setCurrentIndex(1)
            self.activateCommandLine()
            self.activateCommandLine()
        else:
            self.w.multipletAnalysis.setChecked(False)
            self.w.isotopomerAnalysis.setChecked(False)
            self.w.multipletAnalysis.setVisible(False)
            self.w.isotopomerAnalysis.setVisible(False)
            self.w.nmrSpectrum.setTabEnabled(1, False)
            self.w.nmrSpectrum.setTabEnabled(2, False)
            self.w.nmrSpectrum.setStyleSheet("QTabBar::tab::disabled {width: 0; height: 0; margin: 0; padding: 0; border: none;} ")
            self.w.nmrSpectrum.setCurrentIndex(0)

        # end setHsqcAnalysis

    def setMultipletAnalysis(self):
        if (self.w.multipletAnalysis.isChecked() == True):
            self.w.nmrSpectrum.setTabEnabled(3, True)
            self.w.nmrSpectrum.setStyleSheet("QTabBar::tab::disabled {width: 0; height: 0; margin: 0; padding: 0; border: none;} ")
            self.w.nmrSpectrum.setCurrentIndex(3)
        else:
            self.w.nmrSpectrum.setTabEnabled(3, False)
            self.w.nmrSpectrum.setStyleSheet("QTabBar::tab::disabled {width: 0; height: 0; margin: 0; padding: 0; border: none;} ")
            self.w.nmrSpectrum.setCurrentIndex(1)

        # end setMultipletAnalysis

    def setIsotopomerAnalysis(self):
        if (self.w.isotopomerAnalysis.isChecked() == True):
            self.w.nmrSpectrum.setTabEnabled(4, True)
            self.w.nmrSpectrum.setStyleSheet("QTabBar::tab::disabled {width: 0; height: 0; margin: 0; padding: 0; border: none;} ")
            self.w.nmrSpectrum.setCurrentIndex(4)
        else:
            self.w.nmrSpectrum.setTabEnabled(4, False)
            self.w.nmrSpectrum.setStyleSheet("QTabBar::tab::disabled {width: 0; height: 0; margin: 0; padding: 0; border: none;} ")
            self.w.nmrSpectrum.setCurrentIndex(1)

        # end setIsotopomerAnalysis

    def setPreProcessingOptions(self):
        curIdx = self.w.preProcessingSelect.currentIndex()
        self.w.preProcessingWidget.setCurrentIndex(curIdx)
        if self.nd.nmrdat[self.nd.s][0].acq.manufacturer == 'Bruker':
            if self.w.exportMethod.count() == 5:
                self.w.exportMethod.addItem('Bruker Dataset')

        else:
            if self.w.exportMethod.count() == 6:
                self.w.exportMethod.removeItem(5)

        self.plotSpcPreProc()
        # end setPreProcessingOption

    def setPreserveOverallScale(self):
        self.nd.pp.preserveOverallScale = self.w.preserveOverallScale.isChecked()
        # end setPreserveOverallScale

    def setProcPars(self):
        p = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        a = self.nd.nmrdat[self.nd.s][self.nd.e].apc
        self.w.zeroFilling.setText(str(p.nPoints[0]))
        self.w.zeroFilling_2.setText(str(p.nPoints[1]))
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
        self.w.polyOrder.setText(str(p.polyOrder))
        self.w.extrapolationSize.setText(str(p.convExtrapolationSize[0]))
        self.w.windowSize.setText(str(p.convWindowSize[0]))
        self.w.fidOffsetCorrection.setText(str(p.fidOffsetCorrection))
        self.w.windowFunction.setCurrentIndex(p.windowType[0])
        self.w.windowFunction_2.setCurrentIndex(p.windowType[1])
        self.w.phaseCorrection.setCurrentIndex(p.phCorr[0])
        self.w.phaseCorrection_2.setCurrentIndex(p.phCorr[1])
        self.w.waterSuppression.setCurrentIndex(p.waterSuppression)
        self.w.stripTransformStart.setText(str(p.stripStart))
        self.w.stripTransformEnd.setText(str(p.stripEnd))
        self.w.winType.setCurrentIndex(p.convWindowType[0])
        self.w.gibbs.setCurrentIndex(p.gibbsP.get(p.gibbs[0]))
        self.w.gibbs_2.setCurrentIndex(p.gibbsP.get(p.gibbs[1]))
        self.w.rSpc_p0.setText(str(a.rSpc[0]))
        self.w.rSpc_p1.setText(str(a.rSpc[1]))
        self.w.rSpc_p2.setText(str(a.rSpc[2]))
        self.w.rSpc_p3.setText(str(a.rSpc[3]))
        self.w.rSpc_p4.setText(str(a.rSpc[4]))
        self.w.rSpc_p5.setText(str(a.rSpc[5]))
        self.w.rSpc_p6.setText(str(a.rSpc[6]))
        self.w.iSpc_p0.setText(str(a.iSpc[0]))
        self.w.iSpc_p1.setText(str(a.iSpc[1]))
        self.w.iSpc_p2.setText(str(a.iSpc[2]))
        self.w.iSpc_p3.setText(str(a.iSpc[3]))
        self.w.iSpc_p4.setText(str(a.iSpc[4]))
        self.w.iSpc_p5.setText(str(a.iSpc[5]))
        self.w.iSpc_p6.setText(str(a.iSpc[6]))
        self.w.baselineOrder.setCurrentIndex(a.nOrder)
        self.w.baselineCorrection.setCurrentIndex(a.correctBaseline)
        # end setProcPars

    def setPulseProgram(self):
        self.w.pulseProgram.setText(self.nd.nmrdat[self.nd.s][self.nd.e].pulseProgram)
        # end setPulseProgram

    # def setrDolphinExport(self):
    #    self.nd.pp.rDolphinExport = self.w.rDolphinExport.isChecked()
    #
    def setSamplesInComboBox(self):
        self.nd.pp.exportSamplesInRowsCols = self.w.samplesInComboBox.currentIndex()
        # end setSamplesInComboBox

    def setScaleSpectra(self):
        if (self.nd.pp.preProcFill == False):
            if (self.w.scaleSpectra.isChecked() == True):
                self.nd.pp.flagScaleSpectra = True
                self.w.preProcessingSelect.setCurrentIndex(6)
            else:
                self.nd.pp.flagScaleSpectra = False

        # end setScaleSpectra

    def setSegAlignPreProc(self):
        if (self.nd.pp.preProcFill == False):
            nRows = self.w.segAlignTW.rowCount()
            segStart = np.array([])
            segEnd = np.array([])
            tStart = np.array([])
            tEnd = np.array([])
            for k in range(nRows):
                # tStart = np.array([])
                # tEnd   = np.array([])
                try:
                    tStart = np.append(tStart, float(self.w.segAlignTW.item(k, 0).text()))
                    # self.w.segAlignTW.item(k,0).clearContents()
                except:
                    tStart = np.append(tStart, -10000.0)

                try:
                    tEnd = np.append(tEnd, float(self.w.segAlignTW.item(k, 1).text()))
                    # self.w.segAlignTW.item(k,1).clearContents()
                except:
                    tEnd = np.append(tEnd, -10000.0)

            # self.w.segAlignTW.clearContents()
            self.w.segAlignTW.setRowCount(0)
            self.w.segAlignTW.setRowCount(nRows)
            self.nd.pp.preProcFill = True
            for k in np.arange(len(tStart) - 1, -1, -1):  # range(len(tStart)):
                segNumber1 = QTableWidgetItem(2 * k)
                segNumber1.setTextAlignment(QtCore.Qt.AlignHCenter)
                self.w.segAlignTW.setItem(k, 0, segNumber1)
                segNumber2 = QTableWidgetItem(2 * k + 1)
                segNumber2.setTextAlignment(QtCore.Qt.AlignHCenter)
                self.w.segAlignTW.setItem(k, 1, segNumber2)
                if ((tStart[k] > -10000.0) & (tEnd[k] > -10000.0)):
                    tMin = min(tStart[k], tEnd[k])
                    tEnd[k] = max(tStart[k], tEnd[k])
                    tStart[k] = tMin
                    segStart = np.append(segStart, tStart[k])
                    segEnd = np.append(segEnd, tEnd[k])
                    tStart = np.delete(tStart, k)
                    tEnd = np.delete(tEnd, k)

                if (tStart[k] > -10000.0):
                    self.w.segAlignTW.item(k, 0).setText(str(tStart[k]))
                    self.w.segAlignTW.setFocus()
                else:
                    self.w.segAlignTW.item(k, 0).setText("")
                    self.w.segAlignTW.setFocus()

                if (tEnd[k] > -10000.0):
                    self.w.segAlignTW.item(k, 1).setText(str(tEnd[k]))
                    self.w.segAlignTW.setFocus()
                else:
                    self.w.segAlignTW.item(k, 1).setText("")
                    self.w.segAlignTW.setFocus()

            self.nd.pp.preProcFill = False
            sortIdx = np.argsort(segStart)
            self.nd.pp.segStart = segStart[sortIdx]
            self.nd.pp.segEnd = segEnd[sortIdx]
            self.plotSpcPreProc()

        # end setSegAlignPreProc

    def setSegmentalAlignment(self):
        if (self.nd.pp.preProcFill == False):
            if (self.w.segmentalAlignment.isChecked() == True):
                self.nd.pp.flagSegmentalAlignment = True
                self.w.preProcessingSelect.setCurrentIndex(2)
            else:
                self.nd.pp.flagSegmentalAlignment = False

        # end setSegmentalAlignment

    def setSelectClass(self):
        for k in range(len(self.nd.pp.classSelect)):
            self.w.selectClassTW.item(k, 1).setText(self.nd.pp.classSelect[k])

        # end setSelectClass

    def setSymJ(self):
        curIdx = self.w.symJ.currentIndex()
        if (curIdx == 0):
            self.nd.nmrdat[self.nd.s][self.nd.e].proc.symj = True
            self.nd.nmrdat[self.nd.s][self.nd.e].proc.tilt = True
            self.w.tilt.setCurrentIndex(0)
        else:
            self.nd.nmrdat[self.nd.s][self.nd.e].proc.symj = False

        # end setTilt

    def setTilt(self):
        curIdx = self.w.tilt.currentIndex()
        if (curIdx == 0):
            self.nd.nmrdat[self.nd.s][self.nd.e].proc.tilt = True
        else:
            self.nd.nmrdat[self.nd.s][self.nd.e].proc.tilt = False
            self.nd.nmrdat[self.nd.s][self.nd.e].proc.symj = False
            self.w.symJ.setCurrentIndex(1)

        # end setTilt

    def setTitleFile(self):
        self.w.titleFile.setText(self.nd.nmrdat[self.nd.s][self.nd.e].title)
        # end setTitleFile

    def setupProcessingParameters(self):
        self.w.nmrSpectrum.setCurrentIndex(5)
        # end setupProcessingParameters

    def setVarianceStabilisation(self):
        if (self.nd.pp.preProcFill == False):
            if (self.w.varianceStabilisation.isChecked() == True):
                self.nd.pp.flagVarianceStabilisation = True
                self.w.preProcessingSelect.setCurrentIndex(7)
            else:
                self.nd.pp.flagVarianceStabilisation = False

        # end setVarianceStabilisation

    def setVarianceStabilisationOptions(self):
        if self.w.autoScaling.isChecked():
            self.nd.pp.autoScaling = True
            self.nd.pp.paretoScaling = False
            self.nd.pp.gLogTransform = False
        elif self.w.paretoScaling.isChecked():
            self.nd.pp.autoScaling = False
            self.nd.pp.paretoScaling = True
            self.nd.pp.gLogTransform = False
        else:
            self.nd.pp.autoScaling = False
            self.nd.pp.paretoScaling = False
            self.nd.pp.gLogTransform = True

        self.w.lambdaText.setEnabled(self.nd.pp.gLogTransform)
        self.w.y0Text.setEnabled(self.nd.pp.gLogTransform)
        self.w.lambdaLE.setEnabled(self.nd.pp.gLogTransform)
        self.w.y0LE.setEnabled(self.nd.pp.gLogTransform)
        # end setVarianceStabilisationOptions

    def setVarLambda(self):
        self.nd.pp.varLambda = float(self.w.lambdaLE.text())
        # end setVarLambda

    def setVary0(self):
        self.nd.pp.varY0 = float(self.w.y0LE.text())
        # end setVarLambda

    def setPan(self):
        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.setZoomRelease)
        cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
        try:
            self.w.MplWidget.canvas.figure.canvas.toolbar.pan()
        except:
            pass

        self.zoomWasOn = False
        self.panWasOn = True

    def setZoom(self):
        try:
            self.w.MplWidget.canvas.figure.canvas.toolbar.zoom()
        except:
            pass

        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.setZoomRelease)
        self.zoomWasOn = True
        self.panWasOn = False

    def setZoomOff(self):
        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.setZoomRelease)
        cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)

    def setZoomRelease(self, event):
        if (event.button > 1):
            # Right MB click will unzoom the plot
            try:
                self.w.MplWidget.canvas.figure.canvas.toolbar.home()
            except:
                pass

            pyautogui.click(clicks=1)


    def show(self):
        self.w.show()
        # end show

    def showAcquisitionParameters(self):
        self.w.nmrSpectrum.setCurrentIndex(7)
        # end showAcquisitionParameters

    def showAutoBaseline(self):
        self.w.statusBar().clearMessage()
        self.w.statusBar().showMessage("Automatic baseline correction in progress...")
        self.showAcquisitionParameters()
        self.showNMRSpectrum()
        # end showAutoBaseline

    def showAutoPhase(self):
        self.w.statusBar().clearMessage()
        self.w.statusBar().showMessage("Automatic phase correction in progress...")
        self.showAcquisitionParameters()
        self.showNMRSpectrum()
        # end showAutoPhase

    def showConsole(self):
        self.w.nmrSpectrum.setCurrentIndex(11)
        # end showConsole

    def showDisplayParameters(self):
        self.w.nmrSpectrum.setCurrentIndex(6)
        # end showDisplayParameters

    def showHelp(self):
        self.w.nmrSpectrum.setCurrentIndex(12)
        # end showHelp

    def showMainWindow(self):
        if (self.w.isFullScreen() == True):
            self.w.showNormal()
        else:
            self.w.showFullScreen()

        # end showMainWindow

    def showNMRSpectrum(self):
        self.w.nmrSpectrum.setCurrentIndex(0)
        # if (self.w.preprocessing.isChecked() == False):
        #    self.plotSpc()
        # end showNMRSpectrum

    def showPhCorr(self):
        self.w.statusBar().clearMessage()
        self.w.statusBar().showMessage(
            "Left Mouse Button (MB) for ph0, Right MB or Left MB + shift for ph1, Middle MB or Left MB + Cmd to set pivot")
        #    #"Left Mouse Button (MB) for ph0, Right MB or Left MB + shift for ph1, Middle MB or Left MB + Cmd to set pivot        |        Press Alt+p to exit    |   Press Alt+z to zoom")
        self.showAcquisitionParameters()
        self.showNMRSpectrum()
        # end showPhCorr

    def showPhCorr2d(self):
        self.w.statusBar().clearMessage()
        self.w.statusBar().showMessage(
            "2D Interactive Phase Correction")
        #    "Press: Alt+k to pick row/col | Alt+e to empty selection | Alt+r to remove last row/col | Alt+1 for horizontal phase correction | Alt+2 for vertical phase correction | Alt+x to eXit")
        self.showAcquisitionParameters()
        self.showNMRSpectrum()
        # end showPhCorr2d

    def showPhCorr2d_1d(self, dim=0):
        self.w.statusBar().clearMessage()
        self.w.statusBar().showMessage(
            "Left Mouse Button (MB) for ph0, Right MB or Left MB + shift for ph1, Middle MB or Left MB + Cmd to set pivot")
        #    "Left Mouse Button (MB) for ph0, Right MB or Left MB + shift for ph1, Middle MB or Left MB + Cmd to set pivot | Press: Alt+Shift+p to apply phCorr | Alt+Shift+x to cancel | Alt+z to zoom")
        self.showAcquisitionParameters()
        self.showNMRSpectrum()
        # end showPhCorr2d

    def showPhZoom(self):
        self.w.statusBar().clearMessage()
        self.w.statusBar().showMessage(
            "Left Mouse Button (MB) for rectangular zoom, Right MB to unzoom")
        #    "Left Mouse Button (MB) for rectangular zoom, Right MB to unzoom        |        Press Alt+z to exit to phase correction")
        self.showAcquisitionParameters()
        self.showNMRSpectrum()
        # end showPhZoom

    def showPreProcessing(self):
        self.w.preProcessingGroupBox.setHidden(False)
        self.w.preProcessingSelect.setHidden(False)
        self.w.preProcessingWidget.setHidden(False)
        self.w.runPreProcessingButton.setHidden(False)
        self.w.resetPreProcessingButton.setHidden(False)
        self.w.writeScriptButton.setHidden(False)
        self.setExportMethod()
        # self.setSelectClass()
        self.plotSpcPreProc()
        # end showPreProcessing

    def showPulseProgram(self):
        self.w.nmrSpectrum.setCurrentIndex(9)
        # end showPulseProgram

    def showTitleFileInformation(self):
        self.w.nmrSpectrum.setCurrentIndex(8)
        # end showTitleFileInformation

    def showVersion(self):
        self.w.statusBar().clearMessage()
        self.w.statusBar().showMessage("MetaboLabPy " + self.__version__)
        self.showAcquisitionParameters()
        self.showNMRSpectrum()
        # end showVersion

    def startStopPhCorr(self):
        s = self.nd.s
        e = self.nd.e
        if self.zoomWasOn == True:
            try:
                self.setZoom()
                #self.w.MplWidget.canvas.figure.canvas.toolbar.zoom()
            except:
                pass

            self.setZoomOff()

        if self.panWasOn == True:
            try:
                self.setPan()
                self.setZoomOff()
                #self.w.MplWidget.canvas.figure.canvas.toolbar.pan()
            except:
                pass

        if (self.nd.nmrdat[s][e].dim == 1):
            if (self.phCorrActive == False):
                self.phCorr.spc = self.nd.nmrdat[s][e].spc
                self.phCorr.spcMax = max(max(abs(self.phCorr.spc)))
                self.phCorr.pivPoints = self.nd.nmrdat[s][e].ppm2points(self.phCorr.pivot, 0)
                cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.onPhCorrClick)
                cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.onPhCorrRelease)
                self.phCorrActive = True
                self.showPhCorr()
                # self.w.MplWidget.canvas.figure.canvas.toolbar.setEnabled(False)
                self.w.exitPhCorr1d.setVisible(True)
                self.w.zoomPhCorr1d.setVisible(True)
                self.w.exitZoomPhCorr1d.setVisible(False)
                self.updateGUI()
                self.phCorrPlotSpc()
            else:
                cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.onPhCorrClick)
                cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.onPhCorrRelease)
                cid = self.w.MplWidget.canvas.mpl_disconnect(cid)
                cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
                self.phCorrActive = False
                # self.w.MplWidget.canvas.figure.canvas.toolbar.setEnabled(True)
                self.showVersion()
                self.w.exitPhCorr1d.setVisible(False)
                self.w.zoomPhCorr1d.setVisible(False)
                self.w.exitZoomPhCorr1d.setVisible(False)
                self.updateGUI()
                self.plotSpc()
                self.setZoom()

        else:  # dim == 2
            if (self.phCorrActive == False):
                if (self.zoomWasOn == True):
                    self.setZoom()

                if (self.panWasOn == True):
                    self.setPan()

                self.w.pickRowColPhCorr2d.setVisible(True)
                self.w.emptyRowColPhCorr2d.setVisible(True)
                self.w.removeRowColPhCorr2d.setVisible(True)
                self.w.horzPhCorr2d.setVisible(True)
                self.w.vertPhCorr2d.setVisible(True)
                self.w.exitPhCorr2d.setVisible(True)
                self.phCorrActive = True
                self.showPhCorr2d()
            else:
                self.emptyColRow()
                self.w.pickRowColPhCorr2d.setVisible(False)
                self.w.emptyRowColPhCorr2d.setVisible(False)
                self.w.removeRowColPhCorr2d.setVisible(False)
                self.w.horzPhCorr2d.setVisible(False)
                self.w.vertPhCorr2d.setVisible(False)
                self.w.exitPhCorr2d.setVisible(False)
                self.phCorrActive = False
                self.showVersion()
                self.setZoomOff()
                self.setZoom()
                cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.setZoomRelease)

        self.showAcquisitionParameters()
        self.showNMRSpectrum()
        # end startStopPhCorr

    @contextlib.contextmanager
    def stdoutIO(self, stdout=None):
        old = sys.stdout
        if stdout is None:
            stdout = StringIO()
        sys.stdout = stdout
        yield stdout
        sys.stdout = old
        #end stdoutIO

    def tutorials(self):
        #url = "http://beregond.bham.ac.uk/~ludwigc/tutorials"
        fName = os.path.join(os.path.dirname(__file__), "nmr", "web", "tutorials", "index.html")
        url = "file:///" + fName.replace('\\', '/')
        self.w.helpView.setUrl(url)
        self.w.nmrSpectrum.setCurrentIndex(12)
        # end tutorials

    def updateGUI(self):
        s = self.nd.s
        e = self.nd.e
        self.w.setBox.valueChanged.disconnect()
        self.w.expBox.valueChanged.disconnect()
        self.w.expBox.setValue(e + 1)
        self.w.setBox.setValue(s + 1)
        self.w.setBox.valueChanged.connect(lambda: self.changeDataSetExp())
        self.w.expBox.valueChanged.connect(lambda: self.changeDataSetExp())
        self.setDispPars()
        self.setProcPars()
        self.setAcqPars()
        self.setTitleFile()
        self.setPulseProgram()
        self.w.expBox.setValue(e + 1)
        if (self.nd.nmrdat[s][e].dim == 1):
            self.w.preprocessing.setVisible(True)
        else:
            self.w.preprocessing.setVisible(False)

        if self.nd.nmrdat[s][e].dim > 1:
            if self.nd.nmrdat[self.nd.s][self.nd.e].acq.pulProgName.find("hsqc") > 0 or self.nd.nmrdat[self.nd.s][self.nd.e].acq.pulProgName.find("hmqc") > 0:
                self.w.hsqcAnalysis.setVisible(False) # develop set true
            else:
                self.w.hsqcAnalysis.setVisible(False)

        else:
            self.w.hsqcAnalysis.setVisible(False)

        self.w.multipletAnalysis.setVisible(False)
        self.w.isotopomerAnalysis.setVisible(False)
        return "updated GUI"
        # end updateGUI

    def verticalAutoScale(self):
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

        # end verticalAutoScale

    def vertPhCorr2d(self):
        s = self.nd.s
        e = self.nd.e
        self.phCorr.nDims = 2
        self.phCorr.dim = 1
        nLines = len(self.phCorr.spcColPts)
        if nLines > 0:
            npts0 = len(self.nd.nmrdat[s][e].spc)
            npts = len(self.nd.nmrdat[s][e].spc[0])
            self.phCorr.spc = np.zeros((nLines, npts0), dtype='complex')
            spc1 = np.copy(self.nd.nmrdat[s][e].spc)
            spc1 = np.ndarray.transpose(spc1)
            for k in range(nLines):
                spc = np.array([spc1[npts - self.phCorr.spcColPts[k]]])
                spc = self.hilbert(spc)
                self.phCorr.spc[k] = spc[0]

            self.phCorr.ppm = self.nd.nmrdat[s][e].ppm2
            if self.phCorr.pivotPoints2d[1] < 0:
                self.phCorr.pivotPoints2d[1] = int(len(self.phCorr.ppm) / 2)
                self.phCorr.pivot2d[1] = self.nd.nmrdat[s][e].points2ppm(self.phCorr.pivotPoints2d[1], 1)

        self.showPhCorr2d_1d(self.phCorr.dim)
        self.phCorr.spcMax = np.max(np.max(np.abs(self.phCorr.spc)))
        try:
            zwo = True
            self.w.MplWidget.canvas.figure.canvas.toolbar.zoom()
        except:
            pass

        self.setZoomOff()
        self.phCorr.maxPh0 = 90.0
        self.phCorr.maxPh1 = 90.0
        cid = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.onPhCorrClick2d)
        cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.onPhCorrRelease2d)
        #self.w.actionApplyPhCorr.triggered.connect(self.apply2dPhCorr)
        #self.w.actionCancelPhCorr.triggered.connect(self.cancel2dPhCorr)
        self.phCorrPlotSpc2d(False)
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
        self.phCorrPlotSpc2d(False)
        self.showAcquisitionParameters()
        self.showNMRSpectrum()
        # end vertPhCorr2d

    def zeroAcqPars(self):
        self.w.acqPars.setText("")
        # end zeroAcqPars

    def zeroConsole(self):
        self.w.console.setText("")
        # end zeroConsole

    def zeroDispPars(self):
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
        # end zeroDispPars

    def zeroProcPars(self):
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
        # end zeroProcPars

    def zeroPulseProgram(self):
        self.w.pulseProgram.setText("")
        # end zeroPulseProgram

    def zeroScript(self):
        self.w.script.setText("")
        # end zeroConsole

    def zeroTitleFile(self):
        self.w.titleFile.setText("")
        # end zeroTitleFile

    def zoomPhCorr(self):
        print(self.zoomWasOn)
        zwo = False
        pwo = False
        if (self.phCorrActive == True):
            #if self.zoomWasOn == True:
            try:
            #    print("setZoom")
                self.setZoom()
            #    #self.w.MplWidget.canvas.figure.canvas.toolbar.zoom()
            #    zwo = True
            except:
                pass

            print(self.zoom)
            if (self.zoom == False):
                # Enable zoom
                self.zoom = True
                self.showPhZoom()
                if self.phCorr.nDims == 1:
                    self.w.exitPhCorr1d.setVisible(False)
                    self.w.zoomPhCorr1d.setVisible(False)
                    self.w.exitZoomPhCorr1d.setVisible(True)
                else:
                    self.w.zoomPhCorr2d.setVisible(False)
                    self.w.applyPhCorr2d.setVisible(False)
                    self.w.cancelPhCorr2d.setVisible(False)
                    self.w.exitZoomPhCorr2d.setVisible(True)

                #try:
                #    self.setZoom()
                #    #self.w.MplWidget.canvas.figure.canvas.toolbar.zoom()
                #except:
                #    pass
                #
            else:
                # Disable zoom
                self.zoom = False
                if self.phCorr.nDims == 1:
                    self.showPhCorr()
                    self.zoomWasOn = False
                    self.w.exitPhCorr1d.setVisible(True)
                    self.w.zoomPhCorr1d.setVisible(True)
                    self.w.exitZoomPhCorr1d.setVisible(False)
                else:
                    self.showPhCorr2d_1d()
                    self.w.zoomPhCorr2d.setVisible(True)
                    self.w.applyPhCorr2d.setVisible(True)
                    self.w.cancelPhCorr2d.setVisible(True)
                    self.w.exitZoomPhCorr2d.setVisible(False)
                    self.setZoomOff()



        self.showAcquisitionParameters()
        self.showNMRSpectrum()
        # end zoomPhCorr


    def _download_requested(self, download_item) -> None:
        codeOut = io.StringIO()
        codeErr = io.StringIO()
        sys.stdout = codeOut
        sys.stderr = codeErr
        dialog = QFileDialog(self.w)
        path = dialog.getSaveFileName(dialog, "Save File", download_item.path())

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
        self.w.console.setTextColor('Black')
        self.w.console.append(codeOut.getvalue())
        self.w.console.setTextColor('Red')
        self.w.console.append(codeErr.getvalue())
        codeOut.close()
        codeErr.close()
        self.w.nmrSpectrum.setCurrentIndex(11)

    def _download_finished(self) -> None:
        codeOut = io.StringIO()
        codeErr = io.StringIO()
        sys.stdout = codeOut
        sys.stderr = codeErr
        print("Download complete")
        if self.w.autoUnzip.isChecked() == True:
            fName, fExt = os.path.splitext(self.download_item.downloadFileName())
            if fExt == '.zip':
                print('Extracting .zip-file')
                with zipfile.ZipFile(os.path.join(self.download_item.downloadDirectory(), self.download_item.downloadFileName()), 'r') as zip_ref:
                    zip_ref.extractall(self.download_item.downloadDirectory())

                print('.zip-file extraction finished')


        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        self.w.console.setTextColor('Black')
        self.w.console.append(codeOut.getvalue())
        self.w.console.setTextColor('Red')
        self.w.console.append(codeErr.getvalue())
        codeOut.close()
        codeErr.close()
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
    if (len(dd[1]) > 0):
        sys.argv.pop()

    args = vars(ap.parse_args())
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_ShareOpenGLContexts)
    app = QApplication(['pyMetaboLab'])  # sys.argv)
    icon = QIcon()
    pName = os.path.join(os.path.dirname(__file__), "icon")
    icon.addFile(os.path.join(pName, "icon-16.png"), QtCore.QSize(16, 16))
    icon.addFile(os.path.join(pName, "icon-24.png"), QtCore.QSize(24, 24))
    icon.addFile(os.path.join(pName, "icon-32.png"), QtCore.QSize(32, 32))
    icon.addFile(os.path.join(pName, "icon-48.png"), QtCore.QSize(48, 48))
    icon.addFile(os.path.join(pName, "icon-256.png"), QtCore.QSize(256, 256))
    app.setWindowIcon(icon)
    app.setApplicationDisplayName("MetaboLabPy")
    w = main_w()
    if (args["FullScreen"] == True):
        w.w.showFullScreen()

    if (args["noSplash"] == False):
        ##
        # Create and display the splash screen
        pName = os.path.join(os.path.dirname(__file__), "png")
        splash_pix = QPixmap(os.path.join(pName, "metabolabpy.png"))
        splash = QSplashScreen(splash_pix)
        splash.setMask(splash_pix.mask())
        # adding progress bar
        progressBar = QProgressBar(splash)
        splash.show()
        progressBar.show()
        app.processEvents()
        maxTime = 0.5
        maxRange = 30
        timeInc = maxRange
        for i in range(maxRange):
            progressBar.setValue(1.0 * float(i + 1) / float(maxRange))
            # Simulate something that takes time
            time.sleep(maxTime / float(maxRange))
            progressBar.repaint()

        splash.close()
        ## End of splash screen

    if (args["fileName"] != "None"):
        try:
            w.loadFile(args["fileName"])
        except:
            if (args["script"] != None):
                w.openScript(args["script"])
                w.scriptEditor()
                w.execScript()

    else:
        if (args["script"] != None):
            w.openScript(args["script"])
            w.scriptEditor()
            w.execScript()

    cf = nmrConfig.NmrConfig()
    cf.readConfig()
    if cf.mode == 'light':
        qtmodern.styles.light(app)
    else:
        qtmodern.styles.dark(app)

    w.show()
    sys.exit(app.exec_())


if __name__ == "__main__":  # pragma: no cover
    main()
