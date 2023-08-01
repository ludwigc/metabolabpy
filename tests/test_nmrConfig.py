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

    def test_make_config(self):
        cf = nmrConfig.NmrConfig()
        c = cf.make_config()
        self.assertEqual(c['GUI']['auto_plot'], 'no')

    def test_save_config(self):
        cf = nmrConfig.NmrConfig()
        p_name = os.path.join(os.path.dirname(__file__), "data")  # directory of test data
        cf.home_dir = p_name
        cf.config_file = os.path.join(cf.home_dir, cf.f_name)
        cf.save_config()
        self.assertEqual(os.path.isfile(cf.config_file), True)
        os.remove(cf.config_file)

    def test_read_config(self):
        cf = nmrConfig.NmrConfig()
        p_name = os.path.join(os.path.dirname(__file__), "data")  # directory of test data
        cf.home_dir = p_name
        cf.config_file = os.path.join(cf.home_dir, cf.f_name)
        cf.save_config()
        cf.read_config()
        os.remove(cf.config_file)
        self.assertEqual(cf.keep_zoom, True)

    def test_set_auto_plot(self):
        cf = nmrConfig.NmrConfig()
        cf.set_auto_plot("no")
        self.assertEqual(cf.auto_plot, False)
        cf.set_auto_plot("yes")
        self.assertEqual(cf.auto_plot, True)

    def test_set_keep_zoom(self):
        cf = nmrConfig.NmrConfig()
        cf.set_keep_zoom("no")
        self.assertEqual(cf.keep_zoom, False)
        cf.set_keep_zoom("yes")
        self.assertEqual(cf.keep_zoom, True)

    def test_set_font_size(self):
        cf = nmrConfig.NmrConfig()
        cf.set_font_size("18")
        self.assertEqual(cf.font_size, 18.0)
        cf.set_font_size("21")
        self.assertEqual(cf.font_size, 21.0)

    def test_set_phase_reference_colour(self):
        cf = nmrConfig.NmrConfig()
        cf.set_phase_reference_colour("Black")
        self.assertEqual(cf.phase_reference_colour, "Black")
        cf.set_phase_reference_colour("Red")
        self.assertEqual(cf.phase_reference_colour, "Red")

    def test_set_values(self):
        cf = nmrConfig.NmrConfig()
        cf.set_values("auto_plot", "no")
        self.assertEqual(cf.auto_plot, False)
        cf.set_values("auto_plot", "yes")
        self.assertEqual(cf.auto_plot, True)


if __name__ == "__main__":
    unittest.main()
