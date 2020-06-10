#!/usr/bin/env python


"""

test_nmrConfig

author(s): C. Ludwig
origin: 10-06-2020


"""


import unittest
import metabolabpy.nmr.nmrConfig as nmrConfig
import os
import configparser


class nmrConfigTestCase(unittest.TestCase):

    def test_makeConfig(self):
        cf = nmrConfig.NmrConfig()
        c = cf.makeConfig()
        self.assertEqual(c['GUI']['autoplot'], 'yes')

    def test_saveConfig(self):
        cf = nmrConfig.NmrConfig()
        pName = os.path.join(os.path.dirname(__file__), "data")  # directory of test data
        cf.homeDir = pName
        cf.configFile = os.path.join(cf.homeDir, cf.fName)
        cf.saveConfig()
        self.assertEqual(os.path.isfile(cf.configFile), True)
        os.remove(cf.configFile)

    def test_readConfig(self):
        cf = nmrConfig.NmrConfig()
        pName = os.path.join(os.path.dirname(__file__), "data")  # directory of test data
        cf.homeDir = pName
        cf.configFile = os.path.join(cf.homeDir, cf.fName)
        cf.saveConfig()
        cf.readConfig()
        os.remove(cf.configFile)
        self.assertEqual(cf.keepZoom, True)

    def test_set_autoplot(self):
        cf = nmrConfig.NmrConfig()
        cf.set_autoplot("no")
        self.assertEqual(cf.autoPlot, False)
        cf.set_autoplot("yes")
        self.assertEqual(cf.autoPlot, True)

    def test_set_keepzoom(self):
        cf = nmrConfig.NmrConfig()
        cf.set_keepzoom("no")
        self.assertEqual(cf.keepZoom, False)
        cf.set_keepzoom("yes")
        self.assertEqual(cf.keepZoom, True)

    def test_set_fontsize(self):
        cf = nmrConfig.NmrConfig()
        cf.set_fontsize("18")
        self.assertEqual(cf.fontSize, 18.0)
        cf.set_fontsize("21")
        self.assertEqual(cf.fontSize, 21.0)

    def test_set_phasereferencecolour(self):
        cf = nmrConfig.NmrConfig()
        cf.set_phasereferencecolour("Black")
        self.assertEqual(cf.phaseReferenceColour, "Black")
        cf.set_phasereferencecolour("Red")
        self.assertEqual(cf.phaseReferenceColour, "Red")

    def test_setValues(self):
        cf = nmrConfig.NmrConfig()
        cf.setValues("autoplot", "no")
        self.assertEqual(cf.autoPlot, False)
        cf.setValues("autoplot", "yes")
        self.assertEqual(cf.autoPlot, True)


if __name__ == "__main__":
    unittest.main()