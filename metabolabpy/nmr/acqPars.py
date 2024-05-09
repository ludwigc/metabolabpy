"""
NMR spectrum acquisition parameters
"""

import numpy as np
import os
from metabolabpy.nmr import acqRegEx
from metabolabpy.nmr import acqProcparRegEx


class AcqPars:

    def __init__(self):
        self.acqus_text = str('')
        self.acqu2s_text = str('')
        self.acqu3s_text = str('')
        self.byte_order = -1
        self.sw = np.array([0, 0, 0], dtype='float64')
        self.sw_h = np.array([0, 0, 0], dtype='float64')
        self.sfo1 = 0.0
        self.sfo2 = 0.0
        self.sfo3 = 0.0
        self.sfo4 = 0.0
        self.sfo5 = 0.0
        self.sfo6 = 0.0
        self.sfo7 = 0.0
        self.sfo8 = 0.0
        self.bf1 = 0.0
        self.bf2 = 0.0
        self.bf3 = 0.0
        self.bf4 = 0.0
        self.bf5 = 0.0
        self.bf6 = 0.0
        self.bf7 = 0.0
        self.bf8 = 0.0
        self.o1 = 0.0
        self.o2 = 0.0
        self.o3 = 0.0
        self.o4 = 0.0
        self.o5 = 0.0
        self.o6 = 0.0
        self.o7 = 0.0
        self.o8 = 0.0
        self.holder = 0
        self.n_data_points = np.array([0, 0, 0], dtype='int')
        self.aq_mode = 0
        self.decim = 0
        self.dspfvs = 0
        self.group_delay = 0.0
        self.dig_mod = 0
        self.transients = 0
        self.steady_state_scans = 0
        self.relaxation_delay = 0.0
        self.spin_rate = 0.0
        self.nndp = 0.0
        self.pulse_program = str('')
        self.pul_prog_name = str('')
        self.instrument = str('')
        self.data_type = 0
        self.solvent = str('')
        self.probe = str('')
        self.title = str('')
        self.origin = str('')
        self.owner = str('')
        self.meta_info = str('')
        self.aunm = str('')
        self.temperature = -273.15
        self.cnst = np.array([], 'float64')
        self.delay = np.array([], 'float64')
        self.pulse = np.array([], 'float64')
        self.pcpd = np.array([], 'float64')
        self.power_level = np.array([], 'float64')
        self.power_level_watt = np.array([], 'float64')
        self.power_level_max = np.array([], 'float64')
        self.shaped_power = np.array([], 'float64')
        self.shaped_power_watt = np.array([], 'float64')
        self.spoal = np.array([], 'float64')
        self.spoffs = np.array([], 'float64')
        self.cpd_prog = np.array([], dtype='str')
        self.gp_name = np.array([], dtype='str')
        self.vc_list = str('')
        self.vd_list = str('')
        self.vp_list = str('')
        self.va_list = str('')
        self.vt_list = str('')
        self.nuc1 = str('')
        self.nuc2 = str('')
        self.nuc3 = str('')
        self.nuc4 = str('')
        self.nuc5 = str('')
        self.nuc6 = str('')
        self.nuc7 = str('')
        self.nuc8 = str('')
        self.gpx = np.array([], dtype='float64')
        self.gpy = np.array([], dtype='float64')
        self.gpz = np.array([], dtype='float64')
        self.increments = np.array([], dtype='float64')
        self.nus_list = str('')
        self.nus_amount = 0.0
        self.nus_seed = 0
        self.nus_jsp = 0
        self.nus_t2 = 0.0
        self.nus_td = 0
        self.over_flow = 0
        self.pynm = str('')
        self.spc_frequency = np.array([0, 0], dtype='float64')
        self.spc_s_freq = np.array([0, 0], dtype='float64')
        self.spc_nucleus = np.array(['     ', '     '], dtype='str')
        self.spc_offset = np.array([0, 0], dtype='float64')
        self.acq_t0 = np.array([0, 0], dtype='int')
        self.fn_mode = 0
        self.inf = np.array([], dtype='float64')
        self.reg_ex = acqRegEx.AcqRegEx()
        self.reg_ex_varian = acqProcparRegEx.AcqProcparRegEx()
        self.manufacturer = ''
        self.ni = 0
        self.ni2 = 0
        self.np = 1
        self.np2 = 1
        self.phase = np.array([], dtype='int')
        self.phase2 = np.array([], dtype='int')
        self.autopos = ''
        # end __init__

    def __str__(self):  # pragma: no cover
        return self.acqus_text
        # end __str__

    def parse_reg_ex(self):
        try:
            self.sfo1 = float(self.reg_ex.sfo1.findall(self.acqus_text)[0])
        except:
            self.sfo1 = float(self.reg_ex.sfo1i.findall(self.acqus_text)[0])

        try:
            self.sfo2 = float(self.reg_ex.sfo2.findall(self.acqus_text)[0])
        except:
            self.sfo2 = float(self.reg_ex.sfo2i.findall(self.acqus_text)[0])

        try:
            self.sfo3 = float(self.reg_ex.sfo3.findall(self.acqus_text)[0])
        except:
            self.sfo3 = float(self.reg_ex.sfo3i.findall(self.acqus_text)[0])

        try:
            self.sfo4 = float(self.reg_ex.sfo4.findall(self.acqus_text)[0])
        except:
            self.sfo4 = float(self.reg_ex.sfo4i.findall(self.acqus_text)[0])

        try:
            self.sfo5 = float(self.reg_ex.sfo5.findall(self.acqus_text)[0])
        except:
            self.sfo5 = float(self.reg_ex.sfo5i.findall(self.acqus_text)[0])

        try:
            self.sfo6 = float(self.reg_ex.sfo6.findall(self.acqus_text)[0])
        except:
            self.sfo6 = float(self.reg_ex.sfo6i.findall(self.acqus_text)[0])

        try:
            self.sfo7 = float(self.reg_ex.sfo7.findall(self.acqus_text)[0])
        except:
            self.sfo7 = float(self.reg_ex.sfo7i.findall(self.acqus_text)[0])

        try:
            self.sfo8 = float(self.reg_ex.sfo8.findall(self.acqus_text)[0])
        except:
            self.sfo8 = float(self.reg_ex.sfo8i.findall(self.acqus_text)[0])

        self.bf1 = float(self.reg_ex.bf1.findall(self.acqus_text)[0])
        self.bf2 = float(self.reg_ex.bf2.findall(self.acqus_text)[0])
        self.bf3 = float(self.reg_ex.bf3.findall(self.acqus_text)[0])
        self.bf4 = float(self.reg_ex.bf4.findall(self.acqus_text)[0])
        self.bf5 = float(self.reg_ex.bf5.findall(self.acqus_text)[0])
        self.bf6 = float(self.reg_ex.bf6.findall(self.acqus_text)[0])
        self.bf7 = float(self.reg_ex.bf7.findall(self.acqus_text)[0])
        self.bf8 = float(self.reg_ex.bf8.findall(self.acqus_text)[0])
        self.o1 = float(self.reg_ex.o1.findall(self.acqus_text)[0])
        self.o2 = float(self.reg_ex.o2.findall(self.acqus_text)[0])
        self.o3 = float(self.reg_ex.o3.findall(self.acqus_text)[0])
        self.o4 = float(self.reg_ex.o4.findall(self.acqus_text)[0])
        try:
            self.o5 = float(self.reg_ex.o5.findall(self.acqus_text)[0])
            self.o6 = float(self.reg_ex.o6.findall(self.acqus_text)[0])
            self.o7 = float(self.reg_ex.o7.findall(self.acqus_text)[0])
            self.o8 = float(self.reg_ex.o8.findall(self.acqus_text)[0])
        except:
            pass

        self.sw[0] = float(self.reg_ex.sw.findall(self.acqus_text)[0])
        self.sw_h[0] = float(self.reg_ex.sw_h.findall(self.acqus_text)[0])
        self.n_data_points[0] = int(self.reg_ex.td.findall(self.acqus_text)[0])
        self.decim = int(self.reg_ex.decim.findall(self.acqus_text)[0])
        self.dspfvs = int(self.reg_ex.dspfvs.findall(self.acqus_text)[0])
        if self.acqus_text.find("$GRPDLY=") > -1:
            self.group_delay = float(self.reg_ex.grpdly.findall(self.acqus_text)[0][0])
            if self.group_delay < 0:
                self.group_delay = 0

        self.byte_order = int(self.reg_ex.byte_order.findall(self.acqus_text)[0])
        self.aq_mode = int(self.reg_ex.aq_mode.findall(self.acqus_text)[0])
        self.dig_mod = int(self.reg_ex.dig_mod.findall(self.acqus_text)[0])
        self.transients = int(self.reg_ex.transients.findall(self.acqus_text)[0])
        self.steady_state_scans = int(self.reg_ex.steady_state_scans.findall(self.acqus_text)[0])
        self.relaxation_delay = float(self.reg_ex.relaxation_delay.findall(self.acqus_text)[0]) if (
                    len(self.reg_ex.relaxation_delay.findall(self.acqus_text)) > 0) else 0.0
        self.spin_rate = int(self.reg_ex.spin_rate.findall(self.acqus_text)[0])
        self.pul_prog_name = self.reg_ex.pul_prog.findall(self.acqus_text)[0]
        self.aunm = self.reg_ex.aunm.findall(self.acqus_text)[0]
        try:
            self.autopos = self.reg_ex.autopos.findall(self.acqus_text)[0]
        except:
            self.autopos = 'N/A'

        try:
            self.holder = self.reg_ex.holder.findall(self.acqus_text)[0]
        except:
            self.holder = -1

        self.instrument = self.reg_ex.instrument.findall(self.acqus_text)[0]
        self.data_type = int(self.reg_ex.data_type.findall(self.acqus_text)[0])
        self.solvent = self.reg_ex.solvent.findall(self.acqus_text)[0]
        self.probe = self.reg_ex.probe.findall(self.acqus_text)[0]
        self.probe = self.probe.replace('<', '')
        self.probe = self.probe.replace('>', '')
        self.title = self.reg_ex.title.findall(self.acqus_text)[0]
        self.origin = self.reg_ex.origin.findall(self.acqus_text)[0]
        self.owner = self.reg_ex.owner.findall(self.acqus_text)[0]
        try:
            self.meta_info = self.reg_ex.meta_info.findall(self.acqus_text)[0]
        except:
            pass

        try:
            self.temperature = float(self.reg_ex.temperature.findall(self.acqus_text)[0])
        except:
            self.temperature = 300

        dd = self.reg_ex.cnst.search(self.acqus_text)
        dd = self.acqus_text[dd.span()[0]:]
        dd = dd[dd.find('\n') + 1:]
        dd = dd[:dd.find('##$')]
        self.cnst = np.array(dd.split(), dtype='float64')
        dd = self.reg_ex.delay.search(self.acqus_text)
        dd = self.acqus_text[dd.span()[0]:]
        dd = dd[dd.find('\n') + 1:]
        dd = dd[:dd.find('##$')]
        self.delay = np.array(dd.split(), dtype='float64')
        try:
            dd = self.reg_ex.cpd_prog.search(self.acqus_text)
            dd = self.acqus_text[dd.span()[0]:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[:dd.find('##$')]
            self.cpd_prog = np.array(dd.split(), dtype='str')
        except:
            dd = np.resize(np.array([], dtype=str), 9)
            for k in range(9):
                dd[k] = '<>'

            self.cpd_prog = dd

        try:
            dd = self.reg_ex.gp_name.search(self.acqus_text)
            dd = self.acqus_text[dd.span()[0]:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[:dd.find('##$')]
            self.gp_name = np.array(dd.split(), dtype='str')
        except:
            dd = np.resize(np.array([], dtype=str), 32)
            for k in range(32):
                dd[k] = '<>'

            self.gp_name = dd

        dd = self.reg_ex.gpx.search(self.acqus_text)
        dd = self.acqus_text[dd.span()[0]:]
        dd = dd[dd.find('\n') + 1:]
        dd = dd[:dd.find('##$')]
        self.gpx = np.array(dd.split(), dtype='float64')
        dd = self.reg_ex.gpy.search(self.acqus_text)
        dd = self.acqus_text[dd.span()[0]:]
        dd = dd[dd.find('\n') + 1:]
        dd = dd[:dd.find('##$')]
        self.gpy = np.array(dd.split(), dtype='float64')
        dd = self.reg_ex.gpz.search(self.acqus_text)
        dd = self.acqus_text[dd.span()[0]:]
        dd = dd[dd.find('\n') + 1:]
        dd = dd[:dd.find('##$')]
        self.gpz = np.array(dd.split(), dtype='float64')
        dd = self.reg_ex.increments.search(self.acqus_text)
        dd = self.acqus_text[dd.span()[0]:]
        dd = dd[dd.find('\n') + 1:]
        dd = dd[:dd.find('##$')]
        self.increments = np.array(dd.split(), dtype='float64')
        dd = self.reg_ex.pulse.search(self.acqus_text)
        dd = self.acqus_text[dd.span()[0]:]
        dd = dd[dd.find('\n') + 1:]
        dd = dd[:dd.find('##$')]
        self.pulse = np.array(dd.split(), dtype='float64')
        dd = self.reg_ex.pcpd.search(self.acqus_text)
        dd = self.acqus_text[dd.span()[0]:]
        dd = dd[dd.find('\n') + 1:]
        dd = dd[:dd.find('##$')]
        self.pcpd = np.array(dd.split(), dtype='float64')
        dd = self.reg_ex.power_level.search(self.acqus_text)
        dd = self.acqus_text[dd.span()[0]:]
        dd = dd[dd.find('\n') + 1:]
        dd = dd[:dd.find('##$')]
        self.power_level = np.array(dd.split(), dtype='float64')
        if self.acqus_text.find("$PLW=") > -1:
            dd = self.reg_ex.power_level_watt.search(self.acqus_text)
            dd = self.acqus_text[dd.span()[0]:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[:dd.find('##$')]
            self.power_level_watt = np.array(dd.split(), dtype='float64')
            dd = self.reg_ex.power_level_max.search(self.acqus_text)
            dd = self.acqus_text[dd.span()[0]:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[:dd.find('##$')]
            self.power_level_max = np.array(dd.split(), dtype='float64')

        dd = self.reg_ex.shaped_power.search(self.acqus_text)
        dd = self.acqus_text[dd.span()[0]:]
        dd = dd[dd.find('\n') + 1:]
        dd = dd[:dd.find('##$')]
        self.shaped_power = np.array(dd.split(), dtype='float64')
        if self.acqus_text.find("$SPW=") > -1:
            dd = self.reg_ex.shaped_power_watt.search(self.acqus_text)
            dd = self.acqus_text[dd.span()[0]:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[:dd.find('##$')]
            self.shaped_power_watt = np.array(dd.split(), dtype='float64')

        dd = self.reg_ex.spoal.search(self.acqus_text)
        dd = self.acqus_text[dd.span()[0]:]
        dd = dd[dd.find('\n') + 1:]
        dd = dd[:dd.find('##$')]
        self.spoal = np.array(dd.split(), dtype='float64')
        dd = self.reg_ex.spoffs.search(self.acqus_text)
        dd = self.acqus_text[dd.span()[0]:]
        dd = dd[dd.find('\n') + 1:]
        dd = dd[:dd.find('##$')]
        self.spoffs = np.array(dd.split(), dtype='float64')
        if self.acqus_text.find("$INF=") > -1:
            dd = self.reg_ex.inf.search(self.acqus_text)
            dd = self.acqus_text[dd.span()[0]:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[:dd.find('##$')]
            self.inf = np.array(dd.split(), dtype='float64')

        self.vc_list = self.reg_ex.vc_list.findall(self.acqus_text)[0]
        self.vd_list = self.reg_ex.vd_list.findall(self.acqus_text)[0]
        self.vp_list = self.reg_ex.vp_list.findall(self.acqus_text)[0]
        self.va_list = self.reg_ex.va_list.findall(self.acqus_text)[0]
        self.vt_list = self.reg_ex.vt_list.findall(self.acqus_text)[0]
        self.nuc1 = self.reg_ex.nuc1.findall(self.acqus_text)[0]
        self.nuc2 = self.reg_ex.nuc2.findall(self.acqus_text)[0]
        self.nuc3 = self.reg_ex.nuc3.findall(self.acqus_text)[0]
        self.nuc4 = self.reg_ex.nuc4.findall(self.acqus_text)[0]
        self.nuc5 = self.reg_ex.nuc5.findall(self.acqus_text)[0]
        self.nuc6 = self.reg_ex.nuc6.findall(self.acqus_text)[0]
        self.nuc7 = self.reg_ex.nuc7.findall(self.acqus_text)[0]
        self.nuc8 = self.reg_ex.nuc8.findall(self.acqus_text)[0]
        if self.acqus_text.find("$nus_list=") > -1:
            self.nus_list = self.reg_ex.nus_list.findall(self.acqus_text)[0]
            try:
                self.nus_amount = float(self.reg_ex.nus_amount.findall(self.acqus_text)[0])
                self.nus_seed = int(self.reg_ex.nus_seed.findall(self.acqus_text)[0])
                self.nus_jsp = int(self.reg_ex.nus_jsp.findall(self.acqus_text)[0])
                self.nus_t2 = float(self.reg_ex.nus_t2.findall(self.acqus_text)[0][0])
            except:
                pass

        self.over_flow = int(self.reg_ex.over_flow.findall(self.acqus_text)[0])
        if self.acqus_text.find("$PYNM=") > -1:
            self.pynm = self.reg_ex.pynm.findall(self.acqus_text)[0]

        self.spc_frequency[0] = float(self.reg_ex.bf1.findall(self.acqus_text)[0])
        self.spc_s_freq[0] = float(self.reg_ex.sfo1.findall(self.acqus_text)[0])
        self.spc_nucleus[0] = self.reg_ex.nuc1.findall(self.acqus_text)[0]
        self.spc_offset[0] = float(self.reg_ex.o1.findall(self.acqus_text)[0])
        try:
            self.sw[1] = float(self.reg_ex.sw.findall(self.acqu2s_text)[0])
        except:
            pass

        try:
            self.sw_h[1] = float(self.reg_ex.sw_h.findall(self.acqu2s_text)[0])
        except:
            pass

        try:
            self.spc_frequency[1] = float(self.reg_ex.bf1.findall(self.acqu2s_text)[0])
        except:
            pass

        try:
            self.spc_s_freq[1] = float(self.reg_ex.sfo1.findall(self.acqu2s_text)[0])
        except:
            pass

        try:
            self.spc_nucleus[1] = self.reg_ex.nuc1.findall(self.acqu2s_text)[0]
        except:
            pass

        try:
            self.spc_offset[1] = float(self.reg_ex.o1.findall(self.acqu2s_text)[0])
        except:
            pass

        try:
            self.acq_t0[1] = int(self.reg_ex.acq_t0.findall(self.acqu2s_text)[0])
        except:
            pass

        try:
            self.n_data_points[1] = int(self.reg_ex.td.findall(self.acqu2s_text)[0])
        except:
            pass

        try:
            self.nus_td = int(self.reg_ex.nus_td.findall(self.acqu2s_text)[0])
        except:
            pass

        try:
            self.fn_mode = int(self.reg_ex.fn_mode.findall(self.acqu2s_text)[0])
        except:
            pass

        try:
            self.acq_t0[0] = int(self.reg_ex.acq_t0.findall(self.acqus_text)[0])
        except:
            pass

        try:
            self.sw[2] = float(self.reg_ex.sw.findall(self.acqu3s_text)[0]) if (
                        len(self.reg_ex.sw.findall(self.acqu3s_text)) > 0) else 0.0
            self.sw_h[2] = float(self.reg_ex.sw_h.findall(self.acqu3s_text)[0]) if (
                        len(self.reg_ex.sw_h.findall(self.acqu3s_text)) > 0) else 0.0
            self.n_data_points[2] = int(self.reg_ex.td.findall(self.acqu3s_text)[0]) if (
                        len(self.reg_ex.td.findall(self.acqu3s_text)) > 0) else 0

        except:
            pass

        # end parse_reg_ex

    def parse_reg_ex_varian(self):
        dd = self.reg_ex_varian.sfo1.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            self.sfo1 = float(dd)

        dd = self.reg_ex_varian.sfo2.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            self.sfo2 = float(dd)

        dd = self.reg_ex_varian.sfo3.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            self.sfo3 = float(dd)

        dd = self.reg_ex_varian.o1.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            self.o1 = float(dd)

        dd = self.reg_ex_varian.o2.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            self.o2 = float(dd)

        dd = self.reg_ex_varian.o3.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            self.o3 = float(dd)

        self.bf1 = self.sfo1 - self.o1 / 1000.0
        self.bf2 = self.sfo2 - self.o2 / 1000.0
        self.bf3 = self.sfo3 - self.o3 / 1000.0
        self.o1 = 0
        self.o2 = 0
        self.o3 = 0
        dd = self.reg_ex_varian.sw_h.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            self.sw_h[0] = float(dd)
            self.sw[0] = self.sw_h[0] / self.sfo1

        dd = self.reg_ex_varian.sw2_h.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            self.sw_h[1] = float(dd)
            self.sw[1] = self.sw_h[1] / self.sfo2

        dd = self.reg_ex_varian.sw3_h.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            dd2 = self.reg_ex_varian.sw2_h.search(self.acqus_text)
            if hasattr(dd2, 'span'):
                self.sw_h[2] = float(dd)
                self.sw[2] = self.sw_h[2] / self.sfo3
            else:
                self.sw_h[1] = float(dd)
                self.sw[1] = self.sw_h[1] / self.sfo3

        dd = self.reg_ex_varian.td.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            self.n_data_points[0] = int(dd)  # int(int(dd)/2)

        dd = self.reg_ex_varian.phase.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[:dd.find(' ')]
            self.np = int(dd)

        dd = self.reg_ex_varian.phase2.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[:dd.find(' ')]
            self.np2 = int(dd)

        dd = self.reg_ex_varian.ni.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            self.ni = int(dd)

        self.ni = max(self.ni, 1)
        dd = self.reg_ex_varian.ni2.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            self.ni2 = int(dd)

        self.ni2 = max(self.ni2, 1)
        dd = self.reg_ex_varian.phase.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            self.phase = np.array(dd.split(), dtype='int')

        dd = self.reg_ex_varian.phase2.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            self.phase2 = np.array(dd.split(), dtype='int')

        dd = self.reg_ex_varian.transients.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            self.transients = int(dd)

        dd = self.reg_ex_varian.steady_state_scans.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            self.steady_state_scans = int(dd)

        dd = self.reg_ex_varian.nucleus.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            dd = dd.replace(" ", "")
            dd = dd.replace("\"", "")
            self.spc_nucleus = dd

        dd = self.reg_ex_varian.pul_prog_name.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            dd = dd.replace(" ", "")
            dd = dd.replace("\"", "")
            self.pul_prog_name = dd

        dd = self.reg_ex_varian.temperature.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            self.temperature = float(dd) + 273.15

        dd = self.reg_ex_varian.nuc1.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            dd = dd.replace(" ", "")
            dd = dd.replace("\"", "")
            self.nuc1 = dd

        dd = self.reg_ex_varian.nuc2.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            dd = dd.replace(" ", "")
            dd = dd.replace("\"", "")
            self.nuc2 = dd

        dd = self.reg_ex_varian.nuc3.search(self.acqus_text)
        if hasattr(dd, 'span'):
            dd = self.acqus_text[dd.span()[0] + 1:]
            dd = dd[dd.find('\n') + 1:]
            dd = dd[dd.find(' '):dd.find('\n')]
            dd = dd.replace(" ", "")
            dd = dd.replace("\"", "")
            self.nuc3 = dd

        # end parse_reg_ex_varian

    def read(self, spc_dir, acqus_name='acqus'):
        acqus_name = spc_dir + os.sep + acqus_name
        acqu2s_name = spc_dir + os.sep + 'acqu2s'
        acqu3s_name = spc_dir + os.sep + 'acqu3s'
        procpar_name = spc_dir + os.sep + 'procpar'
        procpar_name2 = spc_dir + os.sep + 'PROCPAR'
        if os.path.isfile(acqus_name):
            try:
                f = open(acqus_name, "r")
                self.acqus_text = f.read()
                f.close()

            except:
                f = open(acqus_name, "r", encoding='latin-1')
                self.acqus_text = f.read()
                f.close()

            self.manufacturer = 'Bruker'

        if os.path.isfile(acqu2s_name):
            try:
                f = open(acqu2s_name, "r")
                self.acqu2s_text = f.read()
                f.close()

            except:
                f = open(acqu2s_name, "r", encoding='latin-1')
                self.acqu2s_text = f.read()
                f.close()

            self.manufacturer = 'Bruker'

        if os.path.isfile(acqu3s_name):
            try:
                f = open(acqu3s_name, "r")
                self.acqu3s_text = f.read()
                f.close()

            except:
                f = open(acqu3s_name, "r", encoding='latin-1')
                self.acqu3s_text = f.read()
                f.close()

            self.manufacturer = 'Bruker'

        if os.path.isfile(procpar_name):
            f = open(procpar_name, "r")
            self.acqus_text = f.read()
            f.close()
            self.byte_order = 0
            self.manufacturer = 'Varian'

        if os.path.isfile(procpar_name2):
            f = open(procpar_name2, "r")
            self.acqus_text = f.read()
            f.close()
            self.byte_order = 0
            self.manufacturer = 'Varian'

        if self.manufacturer == 'Bruker':
            self.parse_reg_ex()
            if self.group_delay == 0.0:
                self.set_group_delay()

        if self.manufacturer == 'Varian':
            self.parse_reg_ex_varian()
            self.n_data_points[1] = self.ni * self.ni2 * len(self.phase) * len(self.phase2) * 2
            self.group_delay = 0.0
            self.decim = 0
            self.dspfvs = 0

        # end read

    def set_group_delay(self):
        decims = np.array([2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048])
        dspfirm10 = np.array(
            [179, 201, 533, 709, 1097, 1449, 2225, 2929, 4481, 5889, 8993, 11809, 18017, 23649, 36065, 47329, 72161,
             94689, 144353, 189409, 288737])
        dspfirm11 = np.array(
            [184, 219, 384, 602, 852, 1668, 2312, 3368, 4656, 6768, 9344, 13568, 18560, 27392, 36992, 55040, 73856,
             110336, 147584, 220928, 295040])
        dspfirm12 = np.array(
            [184, 219, 384, 602, 852, 1668, 2292, 3368, 4616, 6768, 9264, 13568, 18560, 27392, 36992, 55040, 73856,
             110336, 147584, 220928, 295040])
        dspfirm13 = np.array([11, 17, 23, 35, 47, 71, 95, 143, 191, 287, 383, 575])
        dspfirm14 = np.array([60, 90, 118, 179, 244, 360, 492, 724, 980, 1444, 1958, 2886, 3912, 5768, 7820, 11532])
        dspfirm15 = np.array([0, 0, 58, 152, 202, 318, 418, 642, 842, 1290, 1690, 2586, 3386])
        dspfirm = [dspfirm10, dspfirm11, dspfirm12, dspfirm13, dspfirm14, dspfirm15]
        addr = np.where(decims == self.decim)[0][0]
        dly = dspfirm[self.dspfvs - 10][addr]
        if self.decim == 3 and self.dspfvs == 15 and self.sw_h > 104000.0:
            self.group_delay = 55.0
        else:
            self.group_delay = dly / self.decim / 2.0

        # end set_group_delay
