#!/usr/bin/env python


"""

test_nmrHsqc

author(s): C. Ludwig
origin: 26-07-2023


"""

import unittest
import metabolabpy.nmr.nmrHsqc as nmrHsqc
import os
import numpy as np


class nmrDataTestCase(unittest.TestCase):

    def test___str__(self):
        hsqc = nmrHsqc.NmrHsqc()
        self.assertEqual(hsqc.__str__(), "NMR HSQC multiplet analysis")
        # end test___str__

    def test_read_metabolite_information(self):
        hsqc = nmrHsqc.NmrHsqc()
        hsqc.read_metabolite_information('L-LacticAcid')
        idx1 = hsqc.metabolite_information.find(':')
        idx2 = hsqc.metabolite_information.find('\n')
        self.assertEqual(hsqc.metabolite_information[idx1 + 1:idx2].replace(' ', ''), 'L-LacticAcid')
        # end test_read_metabolite_information

    def test_set_metabolite_list(self):
        hsqc = nmrHsqc.NmrHsqc()
        hsqc.set_metabolite_list()
        self.assertTrue('L-LacticAcid' in hsqc.metabolite_list)
        # end test_set_metabolite_list

    def test_set_metabolite_information(self):
        metabolite = 'L-LacticAcid'
        hsqc = nmrHsqc.NmrHsqc()
        hsqc.read_metabolite_information(metabolite)
        hsqc.set_metabolite_information(metabolite, hsqc.metabolite_information)
        self.assertTrue(metabolite in hsqc.hsqc_data.keys())
        # end test_set_metabolite_information

    def test_set_peak_information(self):
        metabolite = 'L-LacticAcid'
        hsqc = nmrHsqc.NmrHsqc()
        hsqc.read_metabolite_information(metabolite)
        hsqc.set_metabolite_information(metabolite, hsqc.metabolite_information)
        hsqc.cur_metabolite = metabolite
        hsqc.set_peak_information()
        self.assertEqual(hsqc.hsqc_data[metabolite].spin_systems[0]['c13_shifts'][0][0], 71.2399)
        # end test_set_peak_information


if __name__ == "__main__":
    unittest.main()
