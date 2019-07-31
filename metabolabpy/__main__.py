#!/usr/bin/env python
import argparse
from   PySide2.QtUiTools                  import QUiLoader
from   PySide2.QtCore                     import QFile
from   PySide2.QtWidgets                  import *
from   PySide2.QtGui                      import *
from   PySide2                            import QtGui
from   PySide2                            import QtCore
import matplotlib
matplotlib.use('Qt5Agg')
from   matplotlib.backends.backend_qt5agg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from   matplotlib.figure                  import Figure
import matplotlib.pyplot                  as     pl
import numpy                              as     np
import io
import sys
from metabolabpy.nmr import nmrDataSet 
from metabolabpy.GUI import phCorr 
import matplotlib
import time
import platform
import math
from metabolabpy.nmr import nmrConfig
import os
import traceback


# ------------------ MplWidget ------------------
class MplWidget(QWidget):
    
    def __init__(self, parent = None):
        
        QWidget.__init__(self, parent)
        
        self.canvas = FigureCanvas(Figure())
        
        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)
        self.toolbar    = NavigationToolbar(self.canvas, self)
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
class main_w(object):
    def __init__(self):
        self.cmdBuffer   = np.array([])
        self.cmdIdx      = -1
        self.__version__ = '2019.07101730'
        self.nd          = nmrDataSet.NmrDataSet()
        self.phCorr      = phCorr.PhCorr()
        # load ui; create w
        if(platform.system()=='Darwin'):
            fName = os.path.join(os.path.dirname(__file__),"ui","metabolabpy_mainwindow_mac.ui")
        else:
            fName = os.path.join(os.path.dirname(__file__),"ui","metabolabpy_mainwindow.ui")
            
        self.file = QFile(fName)
        self.file.open(QFile.ReadOnly)
        self.loader = QUiLoader()
        self.loader.registerCustomWidget(MplWidget)
        self.w = self.loader.load(self.file)
        self.zoom = False
        
        self.hidePreProcessing()        
        # connections
        self.w.excludeRegionTW.cellChanged.connect(self.setExcludePreProc)
        self.w.selectClassTW.itemSelectionChanged.connect(self.setPlotPreProc)
        self.w.selectClassTW.cellChanged.connect(self.setChangePreProc)
        self.w.excludeClearButton.clicked.connect(self.selectClearExcludePreProc)
        #self.w.excludeAddButton.clicked.connect(self.selectAddExcludePreProc)
        self.w.selectAllButton.clicked.connect(self.selectAllPreProc)
        self.w.selectEvenButton.clicked.connect(self.selectEvenPreProc)
        self.w.selectOddButton.clicked.connect(self.selectOddPreProc)
        self.w.selectClassButton.clicked.connect(self.selectClassPreProc)
        self.w.selectClassLE.returnPressed.connect(self.selectClassPreProc)
        self.w.cmdLine.returnPressed.connect(self.execCmd)
        self.w.actionActivate_Command_Line.triggered.connect(self.activateCommandLine)
        self.w.actionPrevious_command.triggered.connect(self.previousCommand)
        self.w.actionNext_command.triggered.connect(self.nextCommand)
        self.w.actionCorrect_Phase.triggered.connect(self.startStopPhCorr)
        self.w.actionZoomCorrect_Phase.triggered.connect(self.zoomPhCorr)
        self.w.actionClear.triggered.connect(self.clear)
        self.w.actionRead_Bruker_Spectrum.triggered.connect(lambda: self.readBrukerSpc())
        self.w.preprocessing.stateChanged.connect(lambda: self.setPreProcessing())
        self.w.actionReset.triggered.connect(lambda: self.resetPlot())
        self.w.actionShow_NMR_Spectrum.triggered.connect(lambda: self.showNMRSpectrum())
        self.w.actionSetup_Processing_Parameters.triggered.connect(lambda: self.setupProcessingParameters())
        self.w.actionShow_Display_Parameters.triggered.connect(lambda: self.showDisplayParameters())
        self.w.actionShow_Acquisition_Parameters.triggered.connect(lambda: self.showAcquisitionParameters())
        self.w.actionShow_Title_File_Information.triggered.connect(lambda: self.showTitleFileInformation())
        self.w.actionShow_pulseProgram.triggered.connect(lambda: self.showPulseProgram())
        self.w.actionFourier_Transform.triggered.connect(lambda: self.ft())
        self.w.actionScript_Editor.triggered.connect(lambda: self.scriptEditor())
        self.w.actionChange_to_next_Exp.triggered.connect(lambda: self.changeToNextExp())
        self.w.actionChange_to_previous_Exp.triggered.connect(lambda: self.changeToPreviousExp())
        self.w.actionChange_to_next_DS.triggered.connect(lambda: self.changeToNextDS())
        self.w.actionChange_to_previous_DS.triggered.connect(lambda: self.changeToPreviousDS())
        self.w.actionAutomatic_Phase_Correction.triggered.connect(lambda: self.autophase1d())
        self.w.actionAutomatic_Baseline_Correction.triggered.connect(lambda: self.autobaseline1d())
        self.w.actionSelect_All.triggered.connect(lambda: self.selectPlotAll())
        self.w.actionClear_All.triggered.connect(lambda: self.selectPlotClear())
        self.w.actionConsole.triggered.connect(lambda: self.showConsole())
        self.w.actionToggle_FullScreen.triggered.connect(lambda: self.showMainWindow())
        self.w.setBox.valueChanged.connect(lambda: self.changeDataSetExp())
        self.w.expBox.valueChanged.connect(lambda: self.changeDataSetExp())
        self.w.posCol.currentIndexChanged.connect(lambda: self.getDispPars1())
        self.w.negCol.currentIndexChanged.connect(lambda: self.getDispPars2())
        self.w.posColR.returnPressed.connect(lambda: self.getDispPars3())
        self.w.posColG.returnPressed.connect(lambda: self.getDispPars3())
        self.w.posColB.returnPressed.connect(lambda: self.getDispPars3())
        self.w.negColR.returnPressed.connect(lambda: self.getDispPars3())
        self.w.negColG.returnPressed.connect(lambda: self.getDispPars3())
        self.w.negColB.returnPressed.connect(lambda: self.getDispPars3())
        self.w.nLevels.returnPressed.connect(lambda: self.getDispPars4())
        self.w.minLevel.returnPressed.connect(lambda: self.getDispPars5())
        self.w.maxLevel.returnPressed.connect(lambda: self.getDispPars6())
        self.w.axisType1.currentIndexChanged.connect(lambda: self.getDispPars7())
        self.w.axisType2.currentIndexChanged.connect(lambda: self.getDispPars8())
        self.w.displaySpc.currentIndexChanged.connect(lambda: self.getDispPars9())
        self.w.baselineCorrection.currentIndexChanged.connect(lambda: self.checkBaselineCorrection())
        self.w.baselineOrder.currentIndexChanged.connect(lambda: self. checkBaselineOrder())
        self.w.spcOffset.returnPressed.connect(lambda: self.getDispPars10())
        self.w.spcScale.returnPressed.connect(lambda: self.getDispPars11())
        self.w.fontSize.valueChanged.connect(lambda: self.setFontSize())
        self.w.xLabel.returnPressed.connect(lambda: self.getDispPars12())
        self.w.yLabel.returnPressed.connect(lambda: self.getDispPars13())
        self.w.spcLabel.returnPressed.connect(lambda: self.getDispPars14())
        self.w.preProcessingSelect.currentIndexChanged.connect(lambda: self.setPreProcessingOptions())
        self.w.windowFunction.currentIndexChanged.connect(lambda: self.getProcPars1())
        self.w.windowFunction_2.currentIndexChanged.connect(lambda: self.getProcPars2())
        self.w.phaseCorrection.currentIndexChanged.connect(lambda: self.getProcPars3())
        self.w.phaseCorrection_2.currentIndexChanged.connect(lambda: self.getProcPars4())
        self.w.waterSuppression.currentIndexChanged.connect(lambda: self.getProcPars5())
        self.w.winType.currentIndexChanged.connect(lambda: self.getProcPars6())
        self.w.gibbs.currentIndexChanged.connect(lambda: self.getProcPars7())
        self.w.gibbs_2.currentIndexChanged.connect(lambda: self.getProcPars8())
        self.w.zeroFilling.returnPressed.connect(lambda: self.getProcPars9())
        self.w.zeroFilling_2.returnPressed.connect(lambda: self.getProcPars10())
        self.w.lb.returnPressed.connect(lambda: self.getProcPars11())
        self.w.gb.returnPressed.connect(lambda: self.getProcPars12())
        self.w.ssb.returnPressed.connect(lambda: self.getProcPars13())
        self.w.lb_2.returnPressed.connect(lambda: self.getProcPars14())
        self.w.gb_2.returnPressed.connect(lambda: self.getProcPars15())
        self.w.ssb_2.returnPressed.connect(lambda: self.getProcPars16())
        self.w.ph0.returnPressed.connect(lambda: self.getProcPars17())
        self.w.ph1.returnPressed.connect(lambda: self.getProcPars18())
        self.w.ph0_2.returnPressed.connect(lambda: self.getProcPars19())
        self.w.ph1_2.returnPressed.connect(lambda: self.getProcPars20())
        self.w.polyOrder.returnPressed.connect(lambda: self.getProcPars21())
        self.w.extrapolationSize.returnPressed.connect(lambda: self.getProcPars22())
        self.w.windowSize.returnPressed.connect(lambda: self.getProcPars23())
        self.w.fidOffsetCorrection.returnPressed.connect(lambda: self.getProcPars24())
        self.w.phRefDS.valueChanged.connect(lambda: self.changeDataSetExpPhRef())
        self.w.phRefExp.valueChanged.connect(lambda: self.changeDataSetExpPhRef())
        self.w.phRefColour.currentIndexChanged.connect(lambda: self.getDispPars15())
        self.w.fourierTransformButton.clicked.connect(self.ft)
        self.w.executeScript.clicked.connect(self.execScript)
        # Quit Button
        self.w.quitButton.clicked.connect(self.quit_app)
        self.w.saveButton.clicked.connect(self.saveButton)
        self.w.exportPathSelectButton.clicked.connect(lambda: self.setExportTable())
        self.w.actionQuit.triggered.connect(lambda: self.quit_app())
        self.w.dispPlotButton.clicked.connect(self.plotSpcDisp)
        self.w.openScript.clicked.connect(self.openScript)
        self.w.nmrSpectrum.currentChanged.connect(lambda: self.tabIndexChanged())
        self.showVersion()
        self.keepZoom = False
        self.phCorrActive = False
        self.setFontSize()
        self.cf = nmrConfig.NmrConfig()
        self.cf.readConfig()
        self.w.autoPlot.setChecked(self.cf.autoPlot)
        self.w.keepZoom.setChecked(self.cf.keepZoom)
        self.w.fontSize.setValue(self.cf.fontSize)
        self.w.actionSave_as_Default.triggered.connect(self.saveConfig)
        self.w.actionLoad_Default.triggered.connect(self.loadConfig)
        self.w.actionReset_Config.triggered.connect(self.resetConfig)
        self.w.rSpc_p0.returnPressed.connect(lambda: self.get_rSpc_p0())
        self.w.rSpc_p1.returnPressed.connect(lambda: self.get_rSpc_p1())
        self.w.rSpc_p2.returnPressed.connect(lambda: self.get_rSpc_p2())
        self.w.rSpc_p3.returnPressed.connect(lambda: self.get_rSpc_p3())
        self.w.rSpc_p4.returnPressed.connect(lambda: self.get_rSpc_p4())
        self.w.rSpc_p5.returnPressed.connect(lambda: self.get_rSpc_p5())
        self.w.rSpc_p6.returnPressed.connect(lambda: self.get_rSpc_p6())
        self.w.iSpc_p0.returnPressed.connect(lambda: self.get_iSpc_p0())
        self.w.iSpc_p1.returnPressed.connect(lambda: self.get_iSpc_p1())
        self.w.iSpc_p2.returnPressed.connect(lambda: self.get_iSpc_p2())
        self.w.iSpc_p3.returnPressed.connect(lambda: self.get_iSpc_p3())
        self.w.iSpc_p4.returnPressed.connect(lambda: self.get_iSpc_p4())
        self.w.iSpc_p5.returnPressed.connect(lambda: self.get_iSpc_p5())
        self.w.iSpc_p6.returnPressed.connect(lambda: self.get_iSpc_p6())
        # end __init__


    def activateCommandLine(self):
        self.w.cmdLine.setFocus()
        # end activateCommandLine
        
    def autobaseline1d(self):
        self.showAutoBaseline()
        self.nd.ft()
        self.nd.autoref()
        self.nd.autobaseline1d()
        self.w.baselineCorrection.setCurrentIndex(1)
        self.nd.ft()
        self.nd.baseline1d()
        #self.w.baselineCorrection.setCurrentIndex(1)
        self.setProcPars()
        self.showVersion()
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.changeDataSetExp()
        self.plotSpc()
        # end autobaseline1d
        
    def autobaseline1dAll(self):
        self.showAutoBaseline()
        self.nd.ft()
        self.nd.autoref()
        self.nd.autobaseline1dAll()
        self.w.baselineCorrection.setCurrentIndex(1)
        self.nd.ft()
        self.nd.baseline1d()
        #self.w.baselineCorrection.setCurrentIndex(1)
        self.setProcPars()
        self.showVersion()
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.changeDataSetExp()
        self.plotSpc()
        # end autobaseline1dAll
        
    def autophase1d(self):
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
        # end autophase1d
        
    def autophase1dAll(self):
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
        # end autophase1dAll
        
    def autoref(self):
        self.nd.autoref()
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
        
    def changeDataSetExp(self):
        if(len(self.nd.nmrdat)>0):
            if(len(self.nd.nmrdat[self.nd.s])>0):
                self.keepZoom = self.w.keepZoom.isChecked()
                oldSet = self.nd.s
                oldExp = self.nd.e
                if(self.w.setBox.value()<1):
                    self.w.setBox.setValue(1)
                    
                if(self.w.expBox.value()<1):
                    self.w.expBox.setValue(1)
                    
                if(self.w.setBox.value()>len(self.nd.nmrdat)):
                    self.w.setBox.setValue(len(self.nd.nmrdat))
                    
                self.nd.s = self.w.setBox.value() - 1
                if(self.w.expBox.value()>len(self.nd.nmrdat[self.nd.s])):
                    self.w.expBox.setValue(len(self.nd.nmrdat[self.nd.s]))
                    
                self.nd.e = self.w.expBox.value() - 1
                if(not((oldSet==self.nd.s) and (oldExp==self.nd.e))):
                    self.setDispPars()
                    self.setProcPars()
                    self.setAcqPars()
                    self.setTitleFile()
                    self.setPulseProgram()
                    if(self.phCorrActive==False):
                        if(self.w.autoPlot.isChecked()):
                            self.plotSpc()
                        elif(self.w.nmrSpectrum.currentIndex()==0):
                            print("a")
                            self.plotSpc()
                        
                    else:
                        self.phCorr.spc = self.nd.nmrdat[self.nd.s][self.nd.e].spc
                        self.phCorrPlotSpc()
                    
                                    
                self.keepZoom = False
                
            else:
                self.w.setBox.valueChanged.disconnect()
                self.w.expBox.valueChanged.disconnect()
                self.w.expBox.setValue(0)
                self.w.setBox.setValue(0)
                self.w.setBox.valueChanged.connect(lambda: self.changeDataSetExp())
                self.w.expBox.valueChanged.connect(lambda: self.changeDataSetExp())
                
        else:
            self.w.setBox.valueChanged.disconnect()
            self.w.expBox.valueChanged.disconnect()
            self.w.expBox.setValue(0)
            self.w.setBox.setValue(0)
            self.w.setBox.valueChanged.connect(lambda: self.changeDataSetExp())
            self.w.expBox.valueChanged.connect(lambda: self.changeDataSetExp())
        
        # end changeDataSetExp

    def changeDataSetExpPhRef(self):
        if(len(self.nd.nmrdat)>0):
            s = self.nd.s
            e = self.nd.e
            if(len(self.nd.nmrdat[self.nd.s])>0):
                if(self.w.phRefDS.value()<0):
                    self.w.phRefDS.setValue(0)
                    
                if(self.w.phRefExp.value()<0):
                    self.w.phRefExp.setValue(0)
                    
                if(self.w.phRefDS.value()>len(self.nd.nmrdat)):
                    self.w.phRefExp.setValue(len(self.nd.nmrdat))
                    
                if(self.w.expBox.value()>len(self.nd.nmrdat[self.nd.s])):
                    self.w.expBox.setValue(len(self.nd.nmrdat[self.nd.s]))
                    
                for k in range(len(self.nd.nmrdat)):
                    for l in range(len(self.nd.nmrdat[k])):
                        self.nd.nmrdat[k][l].disp.phRefDS  = self.w.phRefDS.value()
                        self.nd.nmrdat[k][l].disp.phRefExp = self.w.phRefExp.value()
                    
                
            
        
        # end changeDataSetExpPhRef
                        
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
        cbl                                                      = self.w.baselineCorrection.currentIndex()
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.correctBaseline = cbl
        if(cbl == 1):
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
        if(self.w.baselineOrder.isEnabled() == True):
            self.w.rSpc_p0.setEnabled(True)
            self.w.iSpc_p0.setEnabled(True)
            if(blo>0):
                self.w.rSpc_p1.setEnabled(True)
                self.w.iSpc_p1.setEnabled(True)
                
            if(blo>1):
                self.w.rSpc_p2.setEnabled(True)
                self.w.iSpc_p2.setEnabled(True)
                
            if(blo>2):
                self.w.rSpc_p3.setEnabled(True)
                self.w.iSpc_p3.setEnabled(True)
                
            if(blo>3):
                self.w.rSpc_p4.setEnabled(True)
                self.w.iSpc_p4.setEnabled(True)
                
            if(blo>4):
                self.w.rSpc_p5.setEnabled(True)
                self.w.iSpc_p5.setEnabled(True)
                
            if(blo>5):
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
        #end clear
    
    def enableBaseline(self):
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            self.nd.nmrdat[self.nd.s][k].apc.correctBaseline = 1

        self.w.baselineOrder.setCurrentIndex(self.nd.nmrdat[self.nd.s][0].apc.nOrder)
        self.w.baselineCorrection.setCurrentIndex(1)
        return "baselineCorrection enabled"
        # end enableBaseline
        
    def execCmd(self):
        cmdText        = self.w.cmdLine.text()
        if(len(cmdText) > 0):
            self.w.nmrSpectrum.setCurrentIndex(7)
            self.w.cmdLine.setText("")
            self.cmdBuffer = np.append(self.cmdBuffer, cmdText)
            self.cmdIdx    = len(self.cmdBuffer)
            codeOut        = io.StringIO()
            codeErr        = io.StringIO()
            sys.stdout     = codeOut
            sys.stderr     = codeErr
            print(">>> " + cmdText)
            try:
                output = eval(cmdText)
                print(output)
                self.w.console.setTextColor('Black')
                self.w.console.append(codeOut.getvalue())
            except: #(SyntaxError, NameError, TypeError, ZeroDivisionError, AttributeError, ArithmeticError, BufferError, LookupError):
                cmdText2 = "self." + cmdText
                try:
                    output = eval(cmdText2)
                    print(output)
                    self.w.console.setTextColor('Black')
                    self.w.console.append(codeOut.getvalue())
                except:
                    traceback.print_exc()
                    self.w.console.setTextColor('Red')
                    self.w.console.append(codeOut.getvalue())
                    self.w.console.append(codeErr.getvalue())
                        
            
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__
            codeOut.close()
            codeErr.close()
            self.w.console.verticalScrollBar().setValue(self.w.console.verticalScrollBar().maximum())
            
        # end execCmd
        
    def execScript(self):
        codeOut = io.StringIO()
        codeErr = io.StringIO()
        code = self.w.script.toPlainText()
        try:
            exec(code)
        except: # (SyntaxError, NameError, TypeError, ZeroDivisionError, AttributeError):
            traceback.print_exc()
            self.w.console.setTextColor('Red')
            self.w.console.append(codeOut.getvalue())
            self.w.console.append(codeErr.getvalue())
                        
        sys.stdout = codeOut
        sys.stderr = codeErr
        # restore stdout and stderr
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        self.w.console.setTextColor('Blue')
        self.w.console.append('Executing script...\n')
        self.w.console.append(code)
        self.w.console.setTextColor('Black')
        self.w.console.append(codeOut.getvalue())
        self.w.console.setTextColor('Red')
        self.w.console.append(codeErr.getvalue())
        self.w.console.setTextColor('Black')
        codeOut.close()
        codeErr.close()
        self.updateGUI()
        # end execScript
        
    def fillPreProcessingNumbers(self):
        self.nd.pp.preProcFill = True
        nSpc = len(self.nd.pp.classSelect)
        self.w.selectClassTW.setRowCount(nSpc)
        for k in range(nSpc):
            spcNumber = QTableWidgetItem(str(k))
            spcNumber.setTextAlignment(QtCore.Qt.AlignHCenter)
            self.w.selectClassTW.setItem(k, 0, spcNumber)
            #self.w.selectClassTW.setItemSelected(spcNumber, False)
            classNumber = QTableWidgetItem(self.nd.pp.classSelect[k])
            classNumber.setTextAlignment(QtCore.Qt.AlignHCenter)
            self.w.selectClassTW.setItem(k, 1, classNumber)
        
        self.w.selectClassTW.selectAll()
        selIt = self.w.selectClassTW.selectedItems()
        for k in np.arange(len(selIt)-1,-1,-1):
            if(self.w.selectClassTW.selectedItems()[k].column() == 1):
                self.w.selectClassTW.selectedItems()[k].setSelected(False)
                
                
        selIt = self.w.selectClassTW.selectedItems()
        for k in np.arange(len(selIt)-1,-1,-1):
            if(np.isin(self.w.selectClassTW.selectedItems()[k].row(),self.nd.pp.plotSelect)):
                self.w.selectClassTW.selectedItems()[k].setSelected(True)
            else:
                self.w.selectClassTW.selectedItems()[k].setSelected(False)
                
                
        for k in range(len(self.nd.pp.excludeStart)):
            exclNumber1 = QTableWidgetItem(str(2*k))
            exclNumber1.setTextAlignment(QtCore.Qt.AlignHCenter)
            exclNumber2 = QTableWidgetItem(str(2*k+1))
            exclNumber2.setTextAlignment(QtCore.Qt.AlignHCenter)
            self.w.excludeRegionTW.setItem(k, 0, exclNumber1)
            self.w.excludeRegionTW.setItem(k, 1, exclNumber2)
            self.w.excludeRegionTW.item(k,0).setText(str(self.nd.pp.excludeStart[k]))
            self.w.excludeRegionTW.item(k,1).setText(str(self.nd.pp.excludeEnd[k]))
            
        self.nd.pp.preProcFill = False
        #    
        #d  = np.arange(nSpc)
        #dd = np.isin(d, self.nd.pp.plotSelect, invert = True)
        #e  = d[dd]
        #for k in range(len(e)):
        #    self.w.selectClassTW.selectedItems()[2*e[k]].setSelected(False)
        #    
        #self.w.selectClassTW.setRangeSelected(QTableWidgetSelectionRange(0,0,nSpc-1,0),True)
        #for k in range(len(self.nd.pp.plotSelect)):
        #    self.w.selectClassTW.selectionModel().select((self.nd.pp.plotSelect[k],0), True)
        #    #self.w.selectClassTW.setItemSelected(spcNumber, True)
        #    
        # end fillPreProcessingNumbers
        
    def ft(self):
        self.nd.ft()
        if(self.w.baselineCorrection.currentIndex() > 0):
            self.baseline1d()
            
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.changeDataSetExp()
        self.plotSpc()
        # end ft
        
    def getDispPars1(self):
        d            = self.nd.nmrdat[self.nd.s][self.nd.e].disp
        d.posCol     = d.colours.get(self.w.posCol.currentIndex())
        self.nd.nmrdat[self.nd.s][self.nd.e].disp = d
        # end getDispPars1

    def getDispPars2(self):
        d            = self.nd.nmrdat[self.nd.s][self.nd.e].disp
        d.negCol     = d.colours.get(self.w.negCol.currentIndex())
        self.nd.nmrdat[self.nd.s][self.nd.e].disp = d
        # end getDispPars2

    def getDispPars3(self):
        d            = self.nd.nmrdat[self.nd.s][self.nd.e].disp
        posR         = float(self.w.posColR.text())
        posG         = float(self.w.posColG.text())
        posB         = float(self.w.posColB.text())
        negR         = float(self.w.negColR.text())
        negG         = float(self.w.negColG.text())
        negB         = float(self.w.negColB.text())
        d.posColRGB  = (posR, posG, posB)
        d.negColRGB  = (negR, negG, negB)
        self.nd.nmrdat[self.nd.s][self.nd.e].disp = d
        # end getDispPars3

    def getDispPars4(self):
        d            = self.nd.nmrdat[self.nd.s][self.nd.e].disp
        d.nLevels    = round(float(self.w.nLevels.text()))
        # end getDispPars4

        self.nd.nmrdat[self.nd.s][self.nd.e].disp = d
    def getDispPars5(self):
        d            = self.nd.nmrdat[self.nd.s][self.nd.e].disp
        d.minLevel   = float(self.w.minLevel.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].disp = d
        # end getDispPars5

    def getDispPars6(self):
        d            = self.nd.nmrdat[self.nd.s][self.nd.e].disp
        d.maxLevel   = float(self.w.maxLevel.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].disp = d
        # end getDispPars6

    def getDispPars7(self):
        d            = self.nd.nmrdat[self.nd.s][self.nd.e].disp
        d.axisType1  = d.axes.get(self.w.axisType1.currentIndex())
        self.nd.nmrdat[self.nd.s][self.nd.e].disp = d
        # end getDispPars7

    def getDispPars8(self):
        d            = self.nd.nmrdat[self.nd.s][self.nd.e].disp
        d.axisType2  = d.axes.get(self.w.axisType2.currentIndex())
        self.nd.nmrdat[self.nd.s][self.nd.e].disp = d
        # end getDispPars8

    def getDispPars9(self):
        d            = self.nd.nmrdat[self.nd.s][self.nd.e].disp
        d.displaySpc = d.falseTrue.get(self.w.displaySpc.currentIndex())
        self.nd.nmrdat[self.nd.s][self.nd.e].disp = d
        # end getDispPars9

    def getDispPars10(self):
        d            = self.nd.nmrdat[self.nd.s][self.nd.e].disp
        d.spcOffset  = float(self.w.spcOffset.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].disp = d
        # end getDispPars10

    def getDispPars11(self):
        d            = self.nd.nmrdat[self.nd.s][self.nd.e].disp
        d.spcScale   = float(self.w.spcScale.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].disp = d
        # end getDispPars11

    def getDispPars12(self):
        d            = self.nd.nmrdat[self.nd.s][self.nd.e].disp
        d.xLabel     = self.w.xLabel.text()
        self.nd.nmrdat[self.nd.s][self.nd.e].disp = d
        # end getDispPars12

    def getDispPars13(self):
        d            = self.nd.nmrdat[self.nd.s][self.nd.e].disp
        d.yLabel     = self.w.yLabel.text()
        self.nd.nmrdat[self.nd.s][self.nd.e].disp = d
        # end getDispPars13

    def getDispPars14(self):
        d            = self.nd.nmrdat[self.nd.s][self.nd.e].disp
        d.spcLabel   = self.w.spcLabel.text()
        self.nd.nmrdat[self.nd.s][self.nd.e].disp = d
        # end getDispPars14
       
    def getDispPars15(self):
        d            = self.nd.nmrdat[self.nd.s][self.nd.e].disp
        d.phRefCol   = d.colours2.get(self.w.phRefColour.currentIndex())
        for k in range(len(self.nd.nmrdat)):
            for l in range(len(self.nd.nmrdat[k])):
                self.nd.nmrdat[k][l].disp.phRefCol = d.phRefCol
            
        
        # end getDispPars15

    def getProcPars1(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.windowType[0]                           = self.w.windowFunction.currentIndex()
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars1

    def getProcPars2(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.windowType[1]                           = self.w.windowFunction_2.currentIndex()
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars2

    def getProcPars3(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.phCorr[0]                               = self.w.phaseCorrection.currentIndex()
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars3

    def getProcPars4(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.phCorr[1]                               = self.w.phaseCorrection_2.currentIndex()
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars4

    def getProcPars5(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.waterSuppression                        = self.w.waterSuppression.currentIndex()
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars5

    def getProcPars6(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.convWindowType[0]                       = self.w.winType.currentIndex()
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars6

    def getProcPars7(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.gibbs[0]                                = p.gibbsP.get(self.w.gibbs.currentIndex())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars7

    def getProcPars8(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.gibbs[1]                                = p.gibbsP.get(self.w.gibbs_2.currentIndex())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars8

    def getProcPars9(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.nPoints[0]                              = int(self.w.zeroFilling.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars9

    def getProcPars10(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.nPoints[1]                              = int(self.w.zeroFilling_2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars10

    def getProcPars11(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.lb[0]                                   = float(self.w.lb.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars11

    def getProcPars12(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.gb[0]                                   = float(self.w.gb.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars12

    def getProcPars13(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.ssb[0]                                  = float(self.w.ssb.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars13

    def getProcPars14(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.lb[1]                                   = float(self.w.lb_2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars14

    def getProcPars15(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.gb[1]                                   = float(self.w.gb_2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars15

    def getProcPars16(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.ssb[1]                                  = float(self.w.ssb_2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars16

    def getProcPars17(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.ph0[0]                                  = float(self.w.ph0.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars17

    def getProcPars18(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.ph1[0]                                  = float(self.w.ph1.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars18

    def getProcPars19(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.ph0[1]                                  = float(self.w.ph0_2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars19

    def getProcPars20(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.ph1[1]                                  = float(self.w.ph1_2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars20

    def getProcPars21(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.polyOrder                               = int(self.w.polyOrder.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars21

    def getProcPars22(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.convExtrapolationSize[0]                = int(self.w.extrapolationSize.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars22

    def getProcPars23(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.convWindowSize[0]                       = int(self.w.windowSize.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars23

    def getProcPars24(self):
        p                                         = self.nd.nmrdat[self.nd.s][self.nd.e].proc
        p.fidOffsetCorrection                     = int(self.w.fidOffsetCorrection.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].proc = p
        # end getProcPars24
        
    def get_rSpc_p0(self):
        r                                             = self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc
        r[0]                                          = float(self.w.rSpc_p0.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = r
        # end get_rSpc_p0
        
    def get_rSpc_p1(self):
        r                                             = self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc
        r[1]                                          = float(self.w.rSpc_p1.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = r
        # end get_rSpc_p1
        
    def get_rSpc_p2(self):
        r                                             = self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc
        r[2]                                          = float(self.w.rSpc_p2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = r
        # end get_rSpc_p2
        
    def get_rSpc_p3(self):
        r                                             = self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc
        r[3]                                          = float(self.w.rSpc_p3.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = r
        # end get_rSpc_p3
        
    def get_rSpc_p4(self):
        r                                             = self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc
        r[4]                                          = float(self.w.rSpc_p4.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = r
        # end get_rSpc_p4
        
    def get_rSpc_p5(self):
        r                                             = self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc
        r[5]                                          = float(self.w.rSpc_p5.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = r
        # end get_rSpc_p5
        
    def get_rSpc_p6(self):
        r                                             = self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc
        r[6]                                          = float(self.w.rSpc_p6.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = r
        # end get_rSpc_p6
        
    def get_iSpc_p0(self):
        i                                             = self.nd.nmrdat[self.nd.s][self.nd.e].apc.iSpc
        i[0]                                          = float(self.w.iSpc_p0.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = i
        # end get_iSpc_p0
        
    def get_iSpc_p1(self):
        i                                             = self.nd.nmrdat[self.nd.s][self.nd.e].apc.iSpc
        i[1]                                          = float(self.w.iSpc_p1.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = i
        # end get_iSpc_p1
        
    def get_iSpc_p2(self):
        i                                             = self.nd.nmrdat[self.nd.s][self.nd.e].apc.iSpc
        i[2]                                          = float(self.w.iSpc_p2.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = i
        # end get_iSpc_p2
        
    def get_iSpc_p3(self):
        i                                             = self.nd.nmrdat[self.nd.s][self.nd.e].apc.iSpc
        i[3]                                          = float(self.w.iSpc_p3.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = i
        # end get_iSpc_p3
        
    def get_iSpc_p4(self):
        i                                             = self.nd.nmrdat[self.nd.s][self.nd.e].apc.iSpc
        i[4]                                          = float(self.w.iSpc_p4.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = i
        # end get_iSpc_p4
        
    def get_iSpc_p5(self):
        i                                             = self.nd.nmrdat[self.nd.s][self.nd.e].apc.iSpc
        i[5]                                          = float(self.w.iSpc_p5.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = i
        # end get_iSpc_p5
        
    def get_iSpc_p6(self):
        i                                             = self.nd.nmrdat[self.nd.s][self.nd.e].apc.iSpc
        i[6]                                          = float(self.w.iSpc_p6.text())
        self.nd.nmrdat[self.nd.s][self.nd.e].apc.rSpc = i
        # end get_iSpc_p6
        
    def h(self):
        print("Command history: ")
        print(">>><<<")
        for k in range(len(self.cmdBuffer)):
            print(self.cmdBuffer[k])
            
        return(">>><<<")
        # end h
        
    def hidePreProcessing(self):
        self.w.preProcessingGroupBox.setHidden(True)
        self.w.preProcessingSelect.setHidden(True)
        self.w.preProcessingWidget.setHidden(True)
        self.w.runPreProcessingButton.setHidden(True)
        self.w.writeScriptButton.setHidden(True)
        self.plotSpc()
        # end hidePreProcessing

    def loadConfig(self):
        self.cf.readConfig()
        self.w.phRefColour.setCurrentIndex(self.nd.nmrdat[0][0].disp.colours2.get(self.cf.phaseReferenceColour))
        self.w.autoPlot.setChecked(self.cf.autoPlot)
        self.w.keepZoom.setChecked(self.cf.keepZoom)
        self.w.fontSize.setValue(self.cf.fontSize)
        # end loadConfig
        
    def nextCommand(self):
        if(self.w.cmdLine.hasFocus() == True):
            if(self.cmdIdx<len(self.cmdBuffer)):
                self.cmdIdx += 1
                if(self.cmdIdx == len(self.cmdBuffer)):
                    self.w.cmdLine.setText("")
                else:
                    self.w.cmdLine.setText(self.cmdBuffer[self.cmdIdx])
                
            
        
        # end nextCommand

    def onPhCorrClick(self,event):
        s = self.nd.s
        e = self.nd.e
        if(self.zoom==False):
            self.phCorr.spc    = self.nd.nmrdat[s][e].spc
            self.phCorr.spcMax = max(max(abs(self.phCorr.spc)))
            self.w.MplWidget.canvas.toolbar._zoom_mode.__init__()
            if(event.button==1):
                mods = QApplication.queryKeyboardModifiers()
                if(mods==QtCore.Qt.ControlModifier):
                    # set pivot for phase correction
                    self.phCorr.start = event.xdata
                    self.phCorr.pivot     = event.xdata
                    self.phCorr.pivPoints = self.nd.nmrdat[s][e].ppm2points(self.phCorr.pivot,0)
                    
                if(mods==QtCore.Qt.ShiftModifier):
                    # first order phase correction
                    self.phCorr.start = event.ydata
                    
                if(mods==QtCore.Qt.NoModifier):
                    # zero order phase correction
                    self.phCorr.start = event.ydata
    
                if(mods==QtCore.Qt.AltModifier):
                    f.canvas.manager.toolbar.zoom()
                                
            else:
                if(event.button==2):
                    # set pivot for phase correction
                    self.phCorr.start = event.xdata
                    self.phCorr.pivot     = event.xdata
                    self.phCorr.pivPoints = self.nd.nmrdat[s][e].ppm2points(self.phCorr.pivot,0)
                else:
                    # first order phase correction
                    self.phCorr.start = event.ydata
            
    
            cid3 = self.w.MplWidget.canvas.mpl_connect('motion_notify_event', self.onPhCorrDraw)
            
        # end onPhCorrClick
        
    def onPhCorrDraw(self,event):
        if(self.zoom==False):
            s = self.nd.s
            e = self.nd.e
            if((event.xdata!=None) & (event.ydata!=None)):
                self.phCorr.xData = event.xdata
                self.phCorr.yData = event.ydata
                if(event.button==1):
                    mods = QApplication.queryKeyboardModifiers()
                    if(mods==QtCore.Qt.ControlModifier):
                        # set pivot for phase correction
                        self.phCorr.pivot     = event.xdata
                        self.phCorr.pivPoints = self.nd.nmrdat[s][e].ppm2points(self.phCorr.pivot,0)
                    
                    if(mods==QtCore.Qt.ShiftModifier):
                        # first order phase correction
                        ph0 = 0
                        ph1 = self.phCorr.maxPh1*(event.ydata - self.phCorr.start)/self.phCorr.spcMax
                        self.phCorr.spc = self.phase1d(self.nd.nmrdat[s][e].spc,ph0,ph1,self.phCorr.pivPoints)
                        
                    if(mods==QtCore.Qt.NoModifier):
                        # zero order phase correction
                        ph0 = self.phCorr.maxPh0*(event.ydata - self.phCorr.start)/self.phCorr.spcMax
                        ph1 = 0
                        self.phCorr.spc = self.phase1d(self.nd.nmrdat[s][e].spc,ph0,ph1,self.phCorr.pivPoints)
                        
                else:
                    if(event.button==2):
                        # set pivot for phase correction
                        self.phCorr.xData = event.xdata
                        self.phCorr.yData = event.ydata
                        self.phCorr.pivot     = event.xdata
                        self.phCorr.pivPoints = self.nd.nmrdat[s][e].ppm2points(self.phCorr.pivot,0)
                    else:
                        # first order phase correction
                        self.phCorr.xData = event.xdata
                        self.phCorr.yData = event.ydata
                        ph0 = 0
                        ph1 = self.phCorr.maxPh1*(event.ydata - self.phCorr.start)/self.phCorr.spcMax
                        self.phCorr.spc = self.phase1d(self.nd.nmrdat[s][e].spc,ph0,ph1,self.phCorr.pivPoints)
                
            
            self.phCorrPlotSpc()
        
        # end onPhCorrDraw
        
    def onPhCorrRelease(self,event):
        s = self.nd.s
        e = self.nd.e
        if((event.xdata!=None) & (event.ydata!=None)):
            xdata = event.xdata
            ydata = event.ydata
        else:
            xdata = self.phCorr.xData
            ydata = self.phCorr.yData
            
        if(self.zoom==False):
            if(event.button==1):
                mods = QApplication.queryKeyboardModifiers()
                if(mods==QtCore.Qt.ControlModifier):
                    # set pivot for phase correction
                    self.phCorr.pivot     = xdata
                    self.phCorr.pivPoints = self.nd.nmrdat[s][e].ppm2points(self.phCorr.pivot,0)
                    
                if(mods==QtCore.Qt.ShiftModifier):
                    # first order phase correction
                    ph1 = (self.phCorr.maxPh1*(ydata - self.phCorr.start)/self.phCorr.spcMax)
                    ph  =  self.phasesRemovePivot(0.0, ph1, self.phCorr.pivPoints, len(self.phCorr.spc[0]))
                    ph0 = ((self.nd.nmrdat[s][e].proc.ph0[0] + ph[0] + 180.0) % 360.0) - 180.0
                    ph1 =  self.nd.nmrdat[s][e].proc.ph1[0] + ph[1]
                    self.nd.nmrdat[s][e].proc.ph0[0] = ph0
                    self.nd.nmrdat[s][e].proc.ph1[0] = ph1
                    
                if(mods==QtCore.Qt.NoModifier):
                    # zero order phase correction
                    ph0a = (self.phCorr.maxPh0*(ydata - self.phCorr.start)/self.phCorr.spcMax) % 360.0
                    ph1a = 0.0
                    ph  =  self.phasesRemovePivot(ph0a, ph1a, self.phCorr.pivPoints, len(self.phCorr.spc[0]))
                    ph0 = ((self.nd.nmrdat[s][e].proc.ph0[0] + ph[0] + 180.0) % 360.0) - 180.0
                    ph1 =  self.nd.nmrdat[s][e].proc.ph1[0] + ph[1]
                    self.nd.nmrdat[s][e].proc.ph0[0] = ph0
                    self.nd.nmrdat[s][e].proc.ph1[0] = ph1
                        
            else:
                if(event.button==2):
                    # set pivot for phase correction
                    self.phCorr.pivot     = xdata
                    self.phCorr.pivPoints = self.nd.nmrdat[s][e].ppm2points(self.phCorr.pivot,0)
                else:
                    # first order phase correction
                    ph1 = (self.phCorr.maxPh1*(ydata - self.phCorr.start)/self.phCorr.spcMax)
                    ph  =  self.phasesRemovePivot(0.0, ph1, self.phCorr.pivPoints, len(self.phCorr.spc[0]))
                    ph0 = ((self.nd.nmrdat[s][e].proc.ph0[0] + ph[0] + 180.0) % 360.0) - 180.0
                    ph1 =  self.nd.nmrdat[s][e].proc.ph1[0] + ph[1]
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
            if(event.button>1):
                # Right MB click will unzoom the plot
                try:
                    self.w.MplWidget.canvas.figure.canvas.toolbar.home()
                except:
                    pass
                
            
        
        # end onPhCorrRelease
    
    def openScript(self,fName=''):
        if(len(fName)==0):
            fName      = QFileDialog.getOpenFileName()
            fName      = fName[0]
            
        if(len(fName)>0):
            f          = open(fName,'r')
            scriptText = f.read()
            self.w.script.setText(scriptText)
        
        # end openScript
        
    def phCorrPlotSpc(self):
        xlim = self.w.MplWidget.canvas.axes.get_xlim()
        ylim = self.w.MplWidget.canvas.axes.get_ylim()
        d    = self.nd.nmrdat[self.nd.s][self.nd.e].disp
        if(d.posCol=="RGB"):
            posCol = d.posColRGB
        else:
            posCol = d.posCol
            
        if(d.negCol=="RGB"):
            negCol = d.negColRGB
        else:
            negCol = d.negCol
            
        refCol                   = d.phRefCol
        posCol                   = matplotlib.colors.to_hex(posCol)
        negCol                   = matplotlib.colors.to_hex(negCol)
        refCol                   = matplotlib.colors.to_hex(refCol)
        xlabel                   = d.xLabel + " [" + d.axisType1 + "]"
        ylabel                   = d.yLabel + " [" + d.axisType2 + "]"
        if(self.nd.nmrdat[self.nd.s][self.nd.e].dim==1):
            self.w.MplWidget.canvas.axes.clear()
            if((d.phRefDS>0) & (d.phRefExp>0) & (((d.phRefDS-1==self.nd.s) & (d.phRefExp-1==self.nd.e))==False)):
                self.w.MplWidget.canvas.axes.plot(self.nd.nmrdat[d.phRefDS-1][d.phRefExp-1].ppm1, self.nd.nmrdat[d.phRefDS-1][d.phRefExp-1].spc[0].real, color = refCol)
                
            self.w.MplWidget.canvas.axes.plot(self.nd.nmrdat[self.nd.s][self.nd.e].ppm1, self.phCorr.spc[0].real, color = posCol)
            self.w.MplWidget.canvas.axes.plot([self.phCorr.pivot, self.phCorr.pivot], [2.0*self.phCorr.spcMax, -2.0*self.phCorr.spcMax],color = 'r')
            self.w.MplWidget.canvas.axes.set_xlabel(xlabel)
            self.w.MplWidget.canvas.axes.invert_xaxis()
            self.w.MplWidget.canvas.axes.set_xlim(xlim)
            self.w.MplWidget.canvas.axes.set_ylim(ylim)
            self.w.MplWidget.canvas.toolbar.update()
            self.w.MplWidget.canvas.draw()
        
        # end phCorrPlotSpc
        
    def phase1d(self, mat, ph0, ph1, piv):
        npts  = len(mat[0])
        ph0   = -ph0*math.pi/180.0
        ph1   = -ph1*math.pi/180.0
        frac  = np.linspace(0, 1, npts) - float(npts - piv)/float(npts)
        ph    = ph0 + frac*ph1
        mat   = np.cos(ph)*mat.real + np.sin(ph)*mat.imag + 1j*(-np.sin(ph)*mat.real + np.cos(ph)*mat.imag)
        return mat
        #end phase1d
    
    def phasesRemovePivot(self, phc0, phc1, piv, npts):
        phases    = np.array([0.0, 0.0])
        frac      = np.linspace(0, 1, npts) - float(npts - piv)/float(npts)
        ph        = -phc0 - frac*phc1
        phases[0] = -ph[0]
        phases[1] = ph[0] - ph[len(ph)-1]
        return phases
        # end phasesRemovePivot
    
    def plotSpc(self):
        self.keepZoom = self.w.keepZoom.isChecked()
        xlim = self.w.MplWidget.canvas.axes.get_xlim()
        ylim = self.w.MplWidget.canvas.axes.get_ylim()
        self.w.nmrSpectrum.setCurrentIndex(0)
        if(len(self.nd.nmrdat[self.nd.s])==0):
            return
        
        if(len(self.nd.nmrdat[self.nd.s][self.nd.e].spc)==0):
            return
        
        d = self.nd.nmrdat[self.nd.s][self.nd.e].disp
        if(d.posCol=="RGB"):
            posCol = d.posColRGB
        else:
            posCol = d.posCol
            
        if(d.negCol=="RGB"):
            negCol = d.negColRGB
        else:
            negCol = d.negCol
            
        posCol                   = matplotlib.colors.to_hex(posCol)
        negCol                   = matplotlib.colors.to_hex(negCol)
        xlabel                   = d.xLabel + " [" + d.axisType1 + "]"
        ylabel                   = d.yLabel + " [" + d.axisType2 + "]"
        if(self.nd.nmrdat[self.nd.s][self.nd.e].dim==1):
            self.w.MplWidget.canvas.axes.clear()
            for k in range(len(self.nd.nmrdat[self.nd.s])):
                if((k != self.nd.e) and (self.nd.nmrdat[self.nd.s][k].disp.displaySpc == True)):
                    d = self.nd.nmrdat[self.nd.s][k].disp
                    if(d.posCol=="RGB"):
                        posCol = d.posColRGB
                    else:
                        posCol = d.posCol
                
                    if(d.negCol=="RGB"):
                        negCol = d.negColRGB
                    else:
                        negCol = d.negCol
                
                    posCol                   = matplotlib.colors.to_hex(posCol)
                    negCol                   = matplotlib.colors.to_hex(negCol)
                    self.w.MplWidget.canvas.axes.plot(self.nd.nmrdat[self.nd.s][k].ppm1, self.nd.nmrdat[self.nd.s][k].spc[0].real, color = posCol)
                
                
            d = self.nd.nmrdat[self.nd.s][self.nd.e].disp
            if(d.posCol=="RGB"):
                posCol = d.posColRGB
            else:
                posCol = d.posCol
                    
            if(d.negCol=="RGB"):
                negCol = d.negColRGB
            else:
                negCol = d.negCol
                
            posCol                   = matplotlib.colors.to_hex(posCol)
            negCol                   = matplotlib.colors.to_hex(negCol)
            xlabel                   = d.xLabel + " [" + d.axisType1 + "]"
            ylabel                   = d.yLabel + " [" + d.axisType2 + "]"
            self.w.MplWidget.canvas.axes.plot(self.nd.nmrdat[self.nd.s][self.nd.e].ppm1, self.nd.nmrdat[self.nd.s][self.nd.e].spc[0].real, color = posCol)
            self.w.MplWidget.canvas.axes.set_xlabel(xlabel)
            self.w.MplWidget.canvas.axes.autoscale()
            self.w.MplWidget.canvas.axes.invert_xaxis()
            if(self.keepZoom==True):
                self.w.MplWidget.canvas.axes.set_xlim(xlim)
                self.w.MplWidget.canvas.axes.set_ylim(ylim)

            self.w.MplWidget.canvas.toolbar.update()
            self.w.MplWidget.canvas.draw()
            
        else:
            mm                       = self.nd.nmrdat[self.nd.s][self.nd.e].spc.real.max()
            posLev                   = np.linspace( d.minLevel*mm, d.maxLevel*mm,d.nLevels)
            negLev                   = np.linspace(-d.maxLevel*mm,-d.minLevel*mm,d.nLevels)
            self.w.MplWidget.canvas.axes.clear()
            self.w.MplWidget.canvas.axes.contour(self.nd.nmrdat[self.nd.s][self.nd.e].ppm1, self.nd.nmrdat[self.nd.s][self.nd.e].ppm2, self.nd.nmrdat[self.nd.s][self.nd.e].spc.real, posLev, colors = posCol, linestyles = 'solid')
            self.w.MplWidget.canvas.axes.contour(self.nd.nmrdat[self.nd.s][self.nd.e].ppm1, self.nd.nmrdat[self.nd.s][self.nd.e].ppm2, self.nd.nmrdat[self.nd.s][self.nd.e].spc.real, negLev, colors = negCol, linestyles = 'solid')
            self.w.MplWidget.canvas.axes.set_xlabel(xlabel)
            self.w.MplWidget.canvas.axes.set_ylabel(ylabel)
            self.w.MplWidget.canvas.axes.autoscale()
            self.w.MplWidget.canvas.axes.invert_xaxis()
            self.w.MplWidget.canvas.axes.invert_yaxis()
            if(self.keepZoom==True):
                self.w.MplWidget.canvas.axes.set_xlim(xlim)
                self.w.MplWidget.canvas.axes.set_ylim(ylim)

            self.w.MplWidget.canvas.toolbar.update()
            self.w.MplWidget.canvas.draw()
        
        self.keepZoom = False
        # end plotSpc
        
    def plotSpcDisp(self):
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.changeDataSetExp()
        if(self.phCorrActive==False):
            print("b")
            self.plotSpc()
        else:
            self.phCorrPlotSpc()
            
        # end plotSpcDisp
        
    def plotSpcPreProc(self):
        if(len(self.nd.pp.classSelect) == 0):
            self.nd.preProcInit()
            
        self.fillPreProcessingNumbers()    
        sel  = self.w.selectClassTW.selectedIndexes()    
        cls  = np.array([])
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            cls = np.append(cls,self.w.selectClassTW.item(k,1).text())
            
        self.nd.pp.classSelect = cls
        cls2 = np.unique(cls)
        sel2 = np.array([], dtype = 'int')
        for k in range(len(sel)):
            if(sel[k].column() == 0):
                sel2 = np.append(sel2, int(sel[k].row()))
            
        self.nd.pp.plotSelect = sel2
        self.keepZoom = self.w.keepZoom.isChecked()
        xlim = self.w.MplWidget.canvas.axes.get_xlim()
        ylim = self.w.MplWidget.canvas.axes.get_ylim()
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.w.MplWidget.canvas.axes.clear()
        if(self.w.preProcessingWidget.currentIndex() == 1):
            for k in range(len(self.nd.pp.excludeStart)):
                self.w.MplWidget.canvas.axes.axvspan(self.nd.pp.excludeStart[k], self.nd.pp.excludeEnd[k], alpha=self.nd.pp.alpha, color=self.nd.pp.colour)
                
            
        for k in range(len(self.nd.pp.plotSelect)):
            colIdx  = np.where(cls2 == cls[self.nd.pp.plotSelect[k]])[0][0]
            plotCol = matplotlib.colors.to_hex(self.nd.pp.plotColours[colIdx])
            self.w.MplWidget.canvas.axes.plot(self.nd.nmrdat[self.nd.s][self.nd.pp.plotSelect[k]].ppm1, self.nd.nmrdat[self.nd.s][self.nd.pp.plotSelect[k]].spc[0].real, color = plotCol)            
                
        d = self.nd.nmrdat[self.nd.s][self.nd.e].disp
        xlabel = d.xLabel + " [" + d.axisType1 + "]"
        self.w.MplWidget.canvas.axes.set_xlabel(xlabel)
        self.w.MplWidget.canvas.axes.autoscale()
        self.w.MplWidget.canvas.axes.invert_xaxis()
        if(self.keepZoom==True):
            self.w.MplWidget.canvas.axes.set_xlim(xlim)
            self.w.MplWidget.canvas.axes.set_ylim(ylim)

        self.w.MplWidget.canvas.toolbar.update()
        self.w.MplWidget.canvas.draw()
            
    def previousCommand(self):
        if(self.w.cmdLine.hasFocus() == True):
            if(self.cmdIdx>0):
                self.cmdIdx -= 1
                self.w.cmdLine.setText(self.cmdBuffer[self.cmdIdx])
            
        
        # end previousCommand
    def quit_app(self):
        # some actions to perform before actually quitting:
        self.w.close()
        # end quit_app

    def readBrukerSpc(self):
        kz = self.w.keepZoom.isChecked()
        if(len(self.nd.nmrdat[0])==0):
            self.w.keepZoom.setChecked(False)
            
        selected_directory = QFileDialog.getExistingDirectory()
        if(len(selected_directory)>0):
            # Use the selected directory...
            idx     = selected_directory.rfind('/')
            dsName  = selected_directory[:idx]
            expName = selected_directory[idx+1:]
            self.nd.readSpc(dsName, expName)
            self.nd.ft()
            self.nd.autoref()
            self.nd.e = len(self.nd.nmrdat[self.nd.s]) - 1
            self.plotSpc()
            self.w.keepZoom.setChecked(kz)
            self.setProcPars()
            self.setAcqPars()
            self.setTitleFile()
            self.setPulseProgram()
            self.w.expBox.setValue(self.nd.e+1)
            self.setDispPars()
            self.updateGUI()
    
        # end readBrukerSpc
        
    def resetConfig(self):
        self.cf = nmrConfig.NmrConfig()
        self.cf.saveConfig()
        self.loadConfig()
        #end resetConfig

    def resetPlot(self):
        zoomChecked = self.w.keepZoom.isChecked()
        self.w.keepZoom.setChecked(False)
        self.plotSpc()
        if(zoomChecked == True):
            self.w.keepZoom.setChecked(True)
        
        # end resetPlot
        
    def saveButton(self):
        print('save')
        # end save
        
    def saveConfig(self):
        self.cf.autoPlot             = self.w.autoPlot.isChecked()
        self.cf.keepZoom             = self.w.keepZoom.isChecked()
        self.cf.fontSize             = self.w.fontSize.value()
        self.cf.phaseReferenceColour = self.nd.nmrdat[0][0].disp.phRefCol
        self.cf.saveConfig()
        # end saveConfig
        
    def scriptEditor(self):
        self.w.nmrSpectrum.setCurrentIndex(6)
        # end scriptEditor
        
    def selectAllPreProc(self):
        nSpc = len(self.nd.pp.classSelect)
        self.nd.pp.plotSelect = np.arange(nSpc)
        self.fillPreProcessingNumbers()
        self.setPlotPreProc()
        self.plotSpcPreProc()
        self.w.selectClassTW.setFocus()
        # end selectAllPreProc
        
    def selectClassPreProc(self):
        cls  = self.w.selectClassLE.text()
        cls2 = self.nd.pp.classSelect
        sel  = np.array([])
        for k in range(len(cls2)):
            if(cls2[k] == cls):
                sel = np.append(sel, k)
            
        
        if(len(sel) == 0):
            sel = np.arange(len(cls2))
            
        self.nd.pp.plotSelect = sel
        self.fillPreProcessingNumbers()
        self.setPlotPreProc()
        self.plotSpcPreProc()
        self.w.selectClassTW.setFocus()
        # end selectClassPreProc
        
    def selectClearExcludePreProc(self):
        self.nd.pp.preProcFill = True
        for k in range(len(self.nd.pp.excludeStart)):
            self.w.excludeRegionTW.item(k,0).setText("")
            self.w.excludeRegionTW.setFocus()
            self.w.excludeRegionTW.item(k,1).setText("")
            self.w.excludeRegionTW.setFocus()
            
        self.nd.pp.preProcFill = False
        self.nd.pp.excludeStart = np.array([])
        self.nd.pp.excludeEnd   = np.array([])
        self.w.excludeRegionTW.setFocus()
        self.fillPreProcessingNumbers()
        self.w.excludeRegionTW.setFocus()
        self.setPlotPreProc()
        self.w.excludeRegionTW.setFocus()
        self.plotSpcPreProc()
        self.setExcludePreProc()
        # end selectAllPreProc
        
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
            self.nd.nmrdat[self.nd.s][k].disp.displaySpc = True
            
        #self.plotSpc()
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.changeDataSetExp()
        self.plotSpc()
        return "selectPlotAll"
        # end selectPlotAll
        
    def selectPlotClear(self):
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            self.nd.nmrdat[self.nd.s][k].disp.displaySpc = False
        
        #self.plotSpc()
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.changeDataSetExp()
        self.plotSpc()
        return "selectPlotClear"
        # end selectPlotClear
        
    def selectPlotList(self, plotSelect):
        plotSelect = np.array(plotSelect)
        for k in range(len(self.nd.nmrdat[self.nd.s])):
            self.nd.nmrdat[self.nd.s][k].disp.displaySpc = False
        
        plotSelect -= 1
        for k in range(len(plotSelect)):
            if((plotSelect[k]>-1) and (plotSelect[k]<len(self.nd.nmrdat[self.nd.s]))):
                self.nd.nmrdat[self.nd.s][plotSelect[k]].disp.displaySpc = True

        #self.plotSpc()
        self.w.nmrSpectrum.setCurrentIndex(0)
        self.changeDataSetExp()
        return "selectPlotList"
        # end selectPlotList
        
    def setAcqPars(self):
        s       = self.nd.s
        e       = self.nd.e
        a       = self.nd.nmrdat[s][e].acq
        acqStr  = "originalDataset\t" + self.nd.nmrdat[s][e].origDataSet + "\n"
        acqStr += "_________________________________________________________________________________________________________________________________________\n"
        acqStr += "\n"
        acqStr += "metaInfo\t\t"
        for k in range(len(a.title)):
            acqStr += a.title[k] + " "
            
        acqStr += "\n\t\t"
        acqStr += "Origin\t" + a.origin + "\n\t\t"
        acqStr += "Owner\t" + a.owner + "\n"
        acqStr += "_________________________________________________________________________________________________________________________________________\n"
        acqStr += "\n"
        acqStr += "probe\t\t\t" + a.probe + "\n"
        pp      = a.pulProgName
        pp      = pp[1:]
        pp      = pp[:len(pp)-1]
        acqStr += "pulseProgram\t\t" + pp + "\n\n"
        acqStr += "sw\t\t[ppm]\t{:4.2f}\t|\t{:4.2f}\t|\t{:4.2f}\n".format(a.sw[0],a.sw[1],a.sw[2])
        acqStr += "sw\t\t[Hz]\t{:4.2f}\t|\t{:4.2f}\t|\t{:4.2f}\n".format(a.sw_h[0],a.sw_h[1],a.sw_h[2])
        acqStr += "bf1/2/3\t\t[MHz]\t{:4.2f}\t|\t{:4.2f}\t|\t{:4.2f}\n".format(a.bf1,a.bf2,a.bf3)
        acqStr += "sfo1/2/3\t\t[MHz]\t{:4.2f}\t|\t{:4.2f}\t|\t{:4.2f}\n".format(a.sfo1,a.sfo2,a.sfo3)
        acqStr += "o1/2/3\t\t[Hz]\t{:4.2f}\t|\t{:4.2f}\t|\t{:4.2f}\n\n".format(a.o1,a.o2,a.o3)
        acqStr += "nPoints\t\t\t{:0.0f}\t|\t{:0.0f}\t|\t{:0.0f}\n\n".format(a.nDataPoints[0],a.nDataPoints[1],a.nDataPoints[2])
        acqStr += "transients\t\t\t{:0.0f}\n".format(a.transients)
        acqStr += "steadyStateScans\t\t{:0.0f}\n\n".format(a.steadyStateScans)
        acqStr += "groupDelay\t\t\t {:0.0f}\n".format(a.groupDelay)
        acqStr += "decim\t\t\t {:0.0f}\n".format(a.decim)
        acqStr += "dspfvs\t\t\t {:0.0f}\n\n".format(a.dspfvs)
        acqStr += "temperature\t\t[K]\t {:4.2f}\n".format(a.temperature)
        self.w.acqPars.setText(acqStr)
        # end setAcqPars
        
    def setChangePreProc(self):
        if(self.nd.pp.preProcFill == False):
            cls = np.array([])
            for k in range(len(self.nd.pp.classSelect)):
                cls = np.append(cls, self.w.selectClassTW.item(k,1).text())
            
            self.nd.pp.classSelect = cls
            
        # end setChangePreProc
        
    def setDispPars(self):
        d = self.nd.nmrdat[self.nd.s][self.nd.e].disp
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
        if(self.nd.pp.preProcFill == False):
            nRows   = self.w.excludeRegionTW.rowCount()
            exStart = np.array([])
            exEnd   = np.array([])
            for k in range(nRows):
                tStart = np.array([])
                tEnd   = np.array([])
                try:
                    tStart = np.append(tStart, float(self.w.excludeRegionTW.item(k,0).text()))
                    tEnd   = np.append(tEnd,   float(self.w.excludeRegionTW.item(k,1).text()))
                except:
                    pass
                
                if((len(tStart)>0)&(len(tEnd)>0)):
                    tMin   = min(tStart, tEnd)
                    tEnd   = max(tStart, tEnd)
                    tStart = tMin
                    print("tStart/End: %4.2f | %4.2f" % (tStart[0], tEnd[0]))
                    self.w.excludeRegionTW.item(k,0).setText(str(tStart[0]))
                    self.w.excludeRegionTW.item(k,1).setText(str(tEnd[0]))
                    exStart = np.append(exStart, tStart[0])
                    exEnd   = np.append(exEnd,   tEnd[0])
                
            
            self.nd.pp.excludeStart = exStart
            self.nd.pp.excludeEnd   = exEnd
            self.plotSpcPreProc()
            
        # end setExcludePreProc
        
    def setExportTable(self):
        pfName = QFileDialog.getSaveFileName()
        pfName = pfName[0]
        pfName = os.path.split(pfName)
        self.w.exportPath.setText(pfName[0])
        self.w.exportFileName.setText(pfName[1])
        # end setExportTable

    def setFontSize(self):
        fontSize = self.w.fontSize.value()
        f        = self.w.acqPars.font()
        f.setPointSize(fontSize)
        self.w.acqPars.setFont(f)
        self.w.titleFile.setFont(f)
        self.w.pulseProgram.setFont(f)
        self.w.script.setFont(f)
        self.w.console.setFont(f)
        self.w.cmdLine.setFont(f)
        # end setFontSize
        
    def setPhRefExp(self, phRefExp, phRefDS = 1):
        self.w.phRefDS.setValue(phRefDS)
        self.w.phRefExp.setValue(phRefExp)
        # end setPhRefExp
        
    def setPlotPreProc(self):
        if(self.nd.pp.preProcFill == False):
            sel  = np.array([])
            sel  = self.w.selectClassTW.selectedIndexes()
            sel2 = np.array([])
            for k in range(len(sel)):
                if(sel[k].column() == 0):
                    sel2 = np.append(sel2, sel[k].row())
            
        
            self.nd.pp.plotSelect = sel2
            self.plotSpcPreProc()
        
        # end setPlotPreProc
        
    def setPreProcessing(self):
        if(self.w.preprocessing.isChecked() == True):
            self.showPreProcessing()
            #self.nd.preProcInit()
            self.fillPreProcessingNumbers()
        else:
            self.hidePreProcessing()
            
        # end setPreProcessing
        
    def setPreProcessingOptions(self):
        curIdx = self.w.preProcessingSelect.currentIndex()
        self.w.preProcessingWidget.setCurrentIndex(curIdx)
        self.plotSpcPreProc()
        # end setPreProcessingOption
        
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
        
    def setSelectClass(self):
        for k in range(len(self.nd.pp.classSelect)):
            self.w.selectClassTW.item(k,1).setText(self.nd.pp.classSelect[k])
            
        # end setSelectClass
        
    def setTitleFile(self):
        self.w.titleFile.setText(self.nd.nmrdat[self.nd.s][self.nd.e].title)
        # end setTitleFile
        
    def setupProcessingParameters(self):
        self.w.nmrSpectrum.setCurrentIndex(1)
        # end setupProcessingParameters

    def show(self):
        self.w.show()
        # end show
        
    def showAcquisitionParameters(self):
        self.w.nmrSpectrum.setCurrentIndex(3)
        # end showAcquisitionParameters

    def showAutoBaseline(self):
        self.w.statusBar().showMessage("Automatic basline correction in progress...")
        # end showAutoBaseline
        
    def showAutoPhase(self):
        self.w.statusBar().showMessage("Automatic phase correction in progress...")
        # end showAutoPhase
        
    def showConsole(self):
        self.w.nmrSpectrum.setCurrentIndex(7)
        # end showConsole

    def showDisplayParameters(self):
        self.w.nmrSpectrum.setCurrentIndex(2)
        # end showDisplayParameters

    def showMainWindow(self):
        if(self.w.isFullScreen() == True):
            self.w.showNormal()
        else:
            self.w.showFullScreen()
            
        # end showMainWindow
        
    def showNMRSpectrum(self):
        self.w.nmrSpectrum.setCurrentIndex(0)
        if(self.w.preprocessing.isChecked() == False):
            self.plotSpc()
        # end showNMRSpectrum

    def showPhCorr(self):
        self.w.statusBar().showMessage("Left Mouse Button (MB) for ph0, Right MB or Left MB + shift for ph1, Middle MB or Left MB + Cmd to set pivot        |        Press Alt+p to exit    |   Press Alt+z to zoom")
        # end showPhCorr
        
    def showPhZoom(self):
        self.w.statusBar().showMessage("Left Mouse Button (MB) for rectangular zoom, Right MB to unzoom        |        Press Alt+z to exit to phase correction")
        # end showPhZoom
        
    def showPreProcessing(self):
        self.w.preProcessingGroupBox.setHidden(False)
        self.w.preProcessingSelect.setHidden(False)
        self.w.preProcessingWidget.setHidden(False)
        self.w.runPreProcessingButton.setHidden(False)
        self.w.writeScriptButton.setHidden(False)
        #self.setSelectClass()
        self.plotSpcPreProc()
        # end showPreProcessing

    def showPulseProgram(self):
        self.w.nmrSpectrum.setCurrentIndex(5)
        # end showPulseProgram

    def showTitleFileInformation(self):
        self.w.nmrSpectrum.setCurrentIndex(4)
        # end showTitleFileInformation

    def showVersion(self):
        self.w.statusBar().showMessage("MetaboLabPy " + self.__version__)
        # end showVersion
        
    def startStopPhCorr(self):
        s = self.nd.s
        e = self.nd.e
        if(self.w.MplWidget.canvas.figure.canvas.toolbar._active=='ZOOM'):
            try:
                self.w.MplWidget.canvas.figure.canvas.toolbar.zoom()
            except:
                pass
            
            
        if(self.w.MplWidget.canvas.figure.canvas.toolbar._active=='PAN'):
            try:
                self.w.MplWidget.canvas.figure.canvas.toolbar.pan()
            except:
                pass
            
            
        if(self.phCorrActive==False):
            self.phCorr.spc       = self.nd.nmrdat[s][e].spc
            self.phCorr.spcMax    = max(max(abs(self.phCorr.spc)))
            self.phCorr.pivPoints = self.nd.nmrdat[s][e].ppm2points(self.phCorr.pivot,0)
            cid                   = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.onPhCorrClick)
            cid2                  = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.onPhCorrRelease)
            self.phCorrActive     = True
            self.showPhCorr()
            self.phCorrPlotSpc()
            self.w.MplWidget.canvas.figure.canvas.toolbar.setEnabled(False)
        else:
            cid  = self.w.MplWidget.canvas.mpl_connect('button_press_event', self.onPhCorrClick)
            cid2 = self.w.MplWidget.canvas.mpl_connect('button_release_event', self.onPhCorrRelease)
            cid  = self.w.MplWidget.canvas.mpl_disconnect(cid)
            cid2 = self.w.MplWidget.canvas.mpl_disconnect(cid2)
            self.phCorrActive = False
            self.showVersion()
            self.plotSpc()
            self.w.MplWidget.canvas.figure.canvas.toolbar.setEnabled(True)
        
        # end startStopPhCorr

    def tabIndexChanged(self):
        if(self.w.nmrSpectrum.currentIndex()==0):
            if(self.w.preprocessing.isChecked() == False):
                self.plotSpc()
            
        
        # end tabIndexChanged
            
    def updateGUI(self):
        self.w.expBox.setValue(self.nd.e+1)
        self.w.setBox.setValue(self.nd.s+1)
        self.plotSpc()
        self.setDispPars()
        self.setProcPars()
        self.setAcqPars()
        self.setTitleFile()
        self.setPulseProgram()
        self.w.expBox.setValue(self.nd.e+1)
        if(self.nd.nmrdat[self.nd.s][self.nd.e].dim == 1):
            self.w.preprocessing.setEnabled(True)
        else:
            self.w.preprocessing.setEnabled(False)
            
        return "updated GUI"
        #end updateGUI
    
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
        
    def zeroTitleFile(self):
        self.w.titleFile.setText("")
        # end zeroTitleFile
        
    def zoomPhCorr(self):
        if(self.phCorrActive==True):
            if(self.w.MplWidget.canvas.figure.canvas.toolbar._active=='ZOOM'):
                try:
                    self.w.MplWidget.canvas.figure.canvas.toolbar.zoom()
                except:
                    pass
                
                
            if(self.zoom==False):
                # Enable zoom
                self.zoom = True
                self.showPhZoom()
                try:
                    self.w.MplWidget.canvas.figure.canvas.toolbar.zoom()
                except:
                    pass
                
            else:
                # Disable zoom
                self.zoom = False
                self.showPhCorr()
            
        
        # end zoomPhCorr        
    
    
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-s", "--script",      required = False, help = "optional script argument")
    ap.add_argument("-noSplash",           required = False, help = "turn splash screen off",              action = "store_true")
    ap.add_argument("-fs", "--FullScreen", required = False, help = "open applicatin in full screen mode", action = "store_true")
    args  = vars(ap.parse_args())
    app   = QApplication(['pyMetaboLab']) #sys.argv)
    icon  = QIcon()
    pName = os.path.join(os.path.dirname(__file__),"icon")
    icon.addFile(os.path.join(pName,"icon-16.png"), QtCore.QSize(16,16))
    icon.addFile(os.path.join(pName,"icon-24.png"), QtCore.QSize(24,24))
    icon.addFile(os.path.join(pName,"icon-32.png"), QtCore.QSize(32,32))
    icon.addFile(os.path.join(pName,"icon-48.png"), QtCore.QSize(48,48))
    icon.addFile(os.path.join(pName,"icon-256.png"), QtCore.QSize(256,256))
    app.setWindowIcon(icon)
    app.setApplicationDisplayName("MetaboLabPy")
    w      = main_w()
    w.show()
    if(args["FullScreen"]==True):
        w.w.showFullScreen()
    
    if(args["noSplash"]==False):
        ##
        # Create and display the splash screen
        pName = os.path.join(os.path.dirname(__file__),"png")
        splash_pix = QPixmap(os.path.join(pName,"metabolabpy.png"))
        splash = QSplashScreen(splash_pix)
        splash.setMask(splash_pix.mask())
        # adding progress bar
        progressBar = QProgressBar(splash)
        splash.show()
        progressBar.show()
        app.processEvents()
        maxTime  = 0.5
        maxRange = 1000
        timeInc  = maxRange
        for i in range(maxRange):
            progressBar.setValue(1.0*float(i+1)/float(maxRange))
            # Simulate something that takes time
            time.sleep(maxTime/float(maxRange))
            progressBar.repaint()
        
        splash.close()
        ## End of splash screen
    
    if(args["script"]!=None):
        w.openScript(args["script"])
        w.scriptEditor()
        w.execScript()
        
    sys.exit(app.exec_())
    

if __name__ == "__main__":
    main()    
