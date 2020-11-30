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

class nmrDataSetTestCase(unittest.TestCase):
    
    def test_readSpc(self):
        pName = os.path.join(os.path.dirname(__file__),"data","nmrData") # directory of test data set
        eName = "1"                                                      # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()                                 # create nmrDataSet object
        nd.readSpc(pName,eName)                                      # check if Bruker data can be read
        self.assertEqual(len(nd.nmrdat[0][0].fid[0]), 32768)         # check number of data points in fid
        self.assertEqual(len(nd.nmrdat[0][0].spc[0]), 65536)         # check number of data points in 1r

    def test_ft(self):
        pName = os.path.join(os.path.dirname(__file__),"data","nmrData") # directory of test data set
        eName = "1"                                                      # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()                                 # create nmrDataSet object
        nd.readSpc(pName,eName)                                      # check if Bruker data can be read
        nd.ft()
        self.assertEqual(len(nd.nmrdat[0][0].spc[0]), 65536)         # check number of data points in Fourier transformed spectrum

    def test_ftAll(self):
        pName = os.path.join(os.path.dirname(__file__),"data","nmrData") # directory of test data set
        eName = "1"                                                      # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()                                 # create nmrDataSet object
        nd.readSpc(pName,eName)                                      # check if Bruker data can be read
        nd.ftAll()
        self.assertEqual(len(nd.nmrdat[0][0].spc[0]), 65536)         # check number of data points in Fourier transformed spectrum

    def test_load(self):
        nd = nmrDataSet.NmrDataSet()
        dataSetPath = os.path.join(os.path.dirname(__file__), "data", "loadData.mlpy")
        nd.load(dataSetPath)
        self.assertEqual(nd.nmrdat[0][0].refShift[0], 0.0)
        self.assertEqual(nd.nmrdat[0][0].refPoint[0], 10342)

    def test_setZeroFill(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.setZeroFill([131072])
        nd.ft()
        self.assertEqual(len(nd.nmrdat[0][0].spc[0]),131072)  # check number of data points in Fourier transformed spectrum

    def test_autobaseline1d(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.autobaseline1d()

    def test_autobaseline1dAll(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.autobaseline1dAll()

    def test_autophase1d(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.autophase1d()

    def test_autophase1dAll(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.autophase1dAll()

    def test_autoref(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.autoref()
        self.assertEqual(nd.nmrdat[0][0].refPoint[0], 10342)
        self.assertEqual(nd.nmrdat[0][0].refShift[0], 0.0)
        nd = [[]]
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "3"  # 2D HSQC NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.ft()
        nd.autoref()
        self.assertEqual(nd.nmrdat[0][0].refPoint[0], 152)
        self.assertEqual(nd.nmrdat[0][0].refShift[0], 0.0)
        self.assertEqual(nd.nmrdat[0][0].refPoint[1], 1792)
        self.assertAlmostEqual(nd.nmrdat[0][0].refShift[1], 80.0, 4)

    def test_clear(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.clear()
        self.assertEqual(len(nd.nmrdat),1)
        self.assertEqual(len(nd.nmrdat[0]),0)

    def test_excludeRegion(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.pp.excludeStart = [8.9]
        nd.pp.excludeEnd   = [11.0]
        nd.pp.flagExcludeRegion = True
        nd.dataPreProcessing()
        pts = len(nd.nmrdat[0][0].spc[0]) - nd.nmrdat[0][0].ppm2points([nd.pp.excludeStart, nd.pp.excludeEnd],0)
        pts2 = range(int(min(pts)), int(max(pts)))
        ssum = np.sum(nd.nmrdat[0][0].spc[0][pts2].real)
        self.assertEqual(ssum,0.0)
        nd = [[]]
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "2"  # 2D Jres NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.ft()
        nd.autoref()
        nd.pjres(2, 'skyline')
        nd.pp.excludeStart = [8.9]
        nd.pp.excludeEnd   = [10.5]
        nd.pp.flagExcludeRegion = True
        nd.s = 1
        nd.e = 0
        nd.dataPreProcessing()
        pts = len(nd.nmrdat[1][0].spc[0]) - nd.nmrdat[1][0].ppm2points([nd.pp.excludeStart, nd.pp.excludeEnd],0)
        pts2 = range(int(min(pts)), int(max(pts)))
        ssum = np.sum(nd.nmrdat[1][0].spc[0][pts2].real)
        self.assertEqual(ssum,0.0)

    def test_exportDataSet(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.pp.excludeStart = [8.9]
        nd.pp.excludeEnd   = [11.0]
        nd.pp.flagExcludeRegion = True
        nd.pp.bucketPoints = 64
        nd.pp.flagBucketSpectra = True
        nd.dataPreProcessing()
        nd.pp.exportExcelPath = os.path.expanduser("~")
        nd.pp.exportPathName = os.path.expanduser("~")
        nd.pp.exportMetaboAnalystPath = os.path.expanduser("~")
        nd.pp.exportrDolphinPath = os.path.expanduser("~")
        nd.pp.exportBatmanPath = os.path.expanduser("~")
        nd.pp.exportBrukerPath = os.path.expanduser("~")
        nd.pp.exportFileName = "testExport.csv"
        nd.pp.exportSamplesInRowsCols = 0
        nd.pp.classSelect = ["1"]
        nd.pp.exportMethod = 1
        nd.exportDataSet()
        self.assertEqual(os.path.isfile(nd.pp.exportPathName + os.sep + nd.pp.exportFileName), True)
        os.remove(nd.pp.exportPathName + os.sep + nd.pp.exportFileName)
        nd.pp.exportFileName = "testExport.csv"
        nd.pp.exportSamplesInRowsCols = 1
        nd.pp.exportMethod = 1
        nd.exportDataSet()
        self.assertEqual(os.path.isfile(nd.pp.exportPathName + os.sep + nd.pp.exportFileName), True)
        os.remove(nd.pp.exportPathName + os.sep + nd.pp.exportFileName)
        nd.pp.exportExcel = "testExport.xlsx"
        nd.pp.exportSamplesInRowsCols = 0
        nd.pp.exportMethod = 0
        nd.exportDataSet()
        self.assertEqual(os.path.isfile(nd.pp.exportExcelPath + os.sep + nd.pp.exportExcel), True)
        os.remove(nd.pp.exportExcelPath + os.sep + nd.pp.exportExcel)
        nd.pp.exportExcel = "testExport.xlsx"
        nd.pp.exportSamplesInRowsCols = 1
        nd.pp.exportMethod = 0
        nd.exportDataSet()
        self.assertEqual(os.path.isfile(nd.pp.exportExcelPath + os.sep + nd.pp.exportExcel), True)
        os.remove(nd.pp.exportExcelPath + os.sep + nd.pp.exportExcel)
        nd.pp.exportMetaboAnalyst = "metaboAnalystExport.csv"
        nd.pp.exportSamplesInRowsCols = 0
        nd.pp.exportMethod = 2
        nd.exportDataSet()
        self.assertEqual(os.path.isfile(nd.pp.exportMetaboAnalystPath + os.sep + nd.pp.exportMetaboAnalyst), True)
        os.remove(nd.pp.exportMetaboAnalystPath + os.sep + nd.pp.exportMetaboAnalyst)
        nd.pp.exportrDolphin = "rDolphinExport.csv"
        nd.pp.exportSamplesInRowsCols = 0
        nd.pp.exportMethod = 3
        nd.exportDataSet()
        self.assertEqual(os.path.isfile(nd.pp.exportrDolphinPath + os.sep + nd.pp.exportrDolphin), True)
        os.remove(nd.pp.exportrDolphinPath + os.sep + nd.pp.exportrDolphin)
        nd.pp.exportBatman = "batmanExport.txt"
        nd.pp.exportSamplesInRowsCols = 0
        nd.pp.exportMethod = 4
        nd.exportDataSet()
        self.assertEqual(os.path.isfile(nd.pp.exportBatmanPath + os.sep + nd.pp.exportBatman), True)
        os.remove(nd.pp.exportBatmanPath + os.sep + nd.pp.exportBatman)
        nd.pp.exportBruker = "brukerExport"
        nd.pp.exportSamplesInRowsCols = 0
        nd.pp.exportMethod = 5
        nd.exportDataSet()
        self.assertEqual(os.path.isfile(nd.pp.exportBrukerPath + os.sep + nd.pp.exportBruker + os.sep + '1' + os.sep + 'acqus'), True)
        shutil.rmtree(nd.pp.exportBrukerPath + os.sep + nd.pp.exportBruker)

    def test_noiseFiltering(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.pp.noiseStart     = 10.0
        nd.pp.noiseEnd       = 10.5
        nd.pp.noiseThreshold = 8.0
        nd.pp.thLineWidth    = 2.0
        nd.pp.flagNoiseFiltering = True
        nd.dataPreProcessing()
        pts = len(nd.nmrdat[0][0].spc[0]) - nd.nmrdat[0][0].ppm2points([9.5, 11.5],0)
        pts2 = range(int(min(pts)), int(max(pts)))
        ssum = np.sum(nd.nmrdat[0][0].spc[0][pts2].real)
        self.assertEqual(ssum,0.0)

    def test_bucketSpectra(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.pp.bucketPoints = 64
        nd.pp.flagBucketSpectra = True
        nd.dataPreProcessing()
        pts = len(nd.nmrdat[0][0].spc[0])
        self.assertEqual(pts,1024)

    def test_compressBuckets(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.nmrdat[0][0].proc.ph0[0] = -18.49143139318872
        nd.nmrdat[0][0].proc.ph1[0] = -0.7082407011136772
        nd.ft()
        nd.pp.compressStart = [6.5448]
        nd.pp.compressEnd = [7.4705]
        nd.pp.flagCompressBuckets = True
        nd.dataPreProcessing()
        pts = nd.nmrdat[0][0].ppm2points([nd.pp.compressEnd[0], nd.pp.compressStart[0]], 0)
        npts = len(nd.nmrdat[0][0].spc[0])
        spc = nd.nmrdat[0][0].spc[0][npts - pts[0]:npts - pts[1]].real
        self.assertEqual(spc.max(), spc.sum())

    def test_avoidNegativeValues(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.nmrdat[0][0].proc.ph0[0] = -18.49143139318872
        nd.nmrdat[0][0].proc.ph1[0] = -0.7082407011136772
        nd.ft()
        nd.pp.avoidNegativeValues = True
        nd.dataPreProcessing()
        self.assertEqual(nd.nmrdat[0][0].spc[0].real.min(), 0.0)

    def test_segmentalAlignment(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "10"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.nmrdat[0][0].proc.ph0[0] = -85.94556
        nd.nmrdat[0][0].proc.ph1[0] =   6.64201
        nd.nmrdat[0][0].proc.lb[0] = 2.0
        nd.ft()
        nd.nmrdat[0][0].autobaseline1d()
        nd.nmrdat[0][0].autoRef()
        eName = "11"  # 1D NMR data in exp 1
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.nmrdat[0][1].proc.ph0[0] = -79.29513383337031
        nd.nmrdat[0][1].proc.ph1[0] =   7.790498486770286
        nd.nmrdat[0][1].proc.lb[0] = 2.0
        nd.e = 1
        nd.ft()
        nd.nmrdat[0][0].autobaseline1d()
        nd.nmrdat[0][1].autobaseline1d()
        nd.nmrdat[0][0].autoRef()
        nd.nmrdat[0][1].autoRef()
        startPpm = 8.1091
        endPpm = 8.1825
        nd.pp.segStart = [startPpm]
        nd.pp.segEnd = [endPpm]
        nd.pp.flagSegmentalAlignment = True
        nd.pp.segAlignRefSpc = 1
        nd.dataPreProcessing()
        pts1 = len(nd.nmrdat[0][0].spc[0]) - nd.nmrdat[0][0].ppm2points([endPpm, startPpm], 0)
        spc1 = nd.nmrdat[0][0].spc[0][pts1[0]:pts1[1]].real
        spc2 = nd.nmrdat[0][1].spc[0][pts1[0]:pts1[1]].real
        self.assertEqual(np.argmax(spc1), np.argmax(spc2))

    def test_scaleSpectra(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "10"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.nmrdat[0][0].proc.ph0[0] = -85.94556
        nd.nmrdat[0][0].proc.ph1[0] =   6.64201
        nd.nmrdat[0][0].proc.lb[0] = 0.3
        nd.ft()
        nd.nmrdat[0][0].autobaseline1d()
        nd.nmrdat[0][0].autoRef()
        eName = "11"  # 1D NMR data in exp 1
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.nmrdat[0][1].proc.ph0[0] = -79.29513383337031
        nd.nmrdat[0][1].proc.ph1[0] =   7.790498486770286
        nd.nmrdat[0][1].proc.lb[0] = 0.3
        nd.e = 1
        nd.ft()
        nd.pp.flagScaleSpectra = True
        nd.pp.scaleSpectraRefSpc = 1
        self.assertAlmostEqual(nd.nmrdat[0][1].spc[0].real.max() / 10331062121.45446, 1.0, 1)
        nd.pp.scalePQN = True
        nd.dataPreProcessing()
        self.assertAlmostEqual(nd.nmrdat[0][1].spc[0].real.max() / 7644944702.898952, 1.0, 1)
        nd.resetDataPreProcessing()
        self.assertAlmostEqual(nd.nmrdat[0][0].spc[0].real.max() / 9406699558.794441, 1.0, 1)
        nd.pp.scalePQN = False
        nd.pp.flagScaleSpectra = True
        nd.pp.preserveOverallScale = True
        nd.dataPreProcessing()
        self.assertAlmostEqual(nd.nmrdat[0][0].spc[0].real.max() / 11455353501.485264, 1.0, 1)
        nd.resetDataPreProcessing()
        nd.pp.preserveOverallScale = False
        nd.dataPreProcessing()
        self.assertAlmostEqual(nd.nmrdat[0][0].spc[0].real.max() / 0.011113980938806603, 1.0, 1)



    def test_pjres(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "2"  # 2D JresNMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.ft()
        nd.pjres(2, 'skyline')
        self.assertEqual(len(nd.nmrdat[1][0].spc[0]), 8192)

    def test_preProcInit(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.pp.init(1)
        self.assertEqual(nd.pp.classSelect, "1")
        # end test_preProcInit

    def test_preProcInitPlotColours(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.pp.initPlotColours(0.4, 0.8, 0.4)
        self.assertEqual(nd.pp.plotColours[0], (0.0, 0.0, 0.4))
        self.assertEqual(nd.pp.plotColours[-1], (0.8, 0.4, 0.8))
        # end test_preProcInitPlotColours

    def test_preProcSetAlpha(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        dd = nd.pp.setAlpha(0.5)
        self.assertEqual(nd.pp.alpha, 0.5)
        # end test_preProcSetAlpha

    def test_readSpc(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        self.assertEqual(len(nd.nmrdat[0][0].spc[0]), 65536)

    def test_readSpcs(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = ["1", "1"]  # 1D data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpcs([pName], eName)  # check if Bruker data can be read
        self.assertEqual(len(nd.nmrdat[0][0].spc[0]), 65536)
        self.assertEqual(len(nd.nmrdat[0][1].spc[0]), 65536)

    def test_readNMRPipeSpc(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "5"  # 2D HSQC NMRPipe data in exp 5
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readNMRPipeSpc(pName, eName)
        self.assertEqual(len(nd.nmrdat[0][0].spc), 16384)
        self.assertEqual(len(nd.nmrdat[0][0].spc[0]), 922)

    def test_readNMRPipeSpcs(self):
        pName = [os.path.join(os.path.dirname(__file__), "data", "nmrData")]  # directory of test data set
        eName = ["5"]  # 2D HSQC NMRPipe data in exp 5
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readNMRPipeSpcs(pName, eName)
        self.assertEqual(len(nd.nmrdat[0][0].spc), 16384)
        self.assertEqual(len(nd.nmrdat[0][0].spc[0]), 922)

    def test_resetDataPreProcessing(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.pp.excludeStart = [8.9]
        nd.pp.excludeEnd   = [11.0]
        nd.pp.flagExcludeRegion = True
        nd.dataPreProcessing()
        nd.resetDataPreProcessing()
        self.assertEqual(len(np.where(nd.nmrdat[0][0].spc[0].real == 0)[0]), 0)

    def test_save(self):
        nd = nmrDataSet.NmrDataSet()
        dataSetPath = os.path.join(os.path.dirname(__file__), "data", "loadData.mlpy")
        nd.load(dataSetPath)
        dataSetPath = os.path.join(os.path.dirname(__file__), "data", "saveData.mlpy")
        nd.save(dataSetPath)
        self.assertEqual(os.path.exists(dataSetPath), True)
        shutil.rmtree(dataSetPath)

    def test_setGb(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.setGb([15.0])
        self.assertEqual(nd.nmrdat[0][0].proc.gb[0],15.0)

    def test_setLb(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.setLb([15.0])
        self.assertEqual(nd.nmrdat[0][0].proc.lb[0],15.0)

    def test_setPhFromExp(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.setPhFromExp()
        self.assertEqual(nd.nmrdat[0][0].proc.ph0[0],22.0)
        self.assertEqual(nd.nmrdat[0][0].proc.ph1[0],22.0)

    def test_ph0(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.setPh0([15.0])
        self.assertEqual(nd.nmrdat[0][0].proc.ph0[0],15.0)

    def test_ph1(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.setPh1([15.0])
        self.assertEqual(nd.nmrdat[0][0].proc.ph1[0],15.0)

    def test_setSsb(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.setSsb([1.0])
        self.assertEqual(nd.nmrdat[0][0].proc.ssb[0],1.0)

    def test_setWindowType(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.setWindowType([2])
        self.assertEqual(nd.nmrdat[0][0].proc.windowType[0],2)

    def test_shiftRef(self):
        pName = os.path.join(os.path.dirname(__file__), "data", "nmrData")  # directory of test data set
        eName = "1"  # 1D NMR data in exp 1
        nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
        nd.readSpc(pName, eName)  # check if Bruker data can be read
        nd.shiftRef()
        self.assertEqual(nd.nmrdat[0][0].refShift[0],0.0)




if __name__ == "__main__":
    unittest.main()