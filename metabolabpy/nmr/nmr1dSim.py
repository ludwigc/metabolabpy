"""
NMR 1d Simulation

"""

import numpy as np
from metabolabpy.nmr import nmr1dData  # pragma: no cover
import os


class Nmr1dSim:

    def __init__(self):
        self.range_h = 0.1            # [ppm]
        self.autopick_range_h = 0.07  # [ppm]
        self.threshold = 0.05         # [%]
        self.j_hh = 10                # [Hz]
        self.metabolite_list = []
        self.metabolite_information = ''
        self.nmr1d_data = {}
        self.cur_metabolite = ''
        self.echo_time = 1.95
        self.delta_h1_low = 0.15
        self.delta_h1_high = 0.3
        # end __init__

    def __str__(self):
        str_str = "NMR HSQC multiplet analysis"
        return str_str
        # end __str__

    def read_metabolite_information(self, metabolite_name=''):
        if len(metabolite_name) == 0:
            self.metabolite_information = ''
            return

        #print(os.path.dirname(__file__))
        file_name = os.path.join(os.path.dirname(__file__), 'metabolites', metabolite_name + '.mlInfo')
        fid = open(file_name, 'r')
        self.metabolite_information = fid.read()
        fid.close()
        # end read_metabolite_information

    def set_metabolite_list(self): # GUI item, not finished
        met_list = next(os.walk(os.path.join(os.path.dirname(__file__), 'metabolites')), (None, None, []))[2]
        met_list = sorted(met_list)
        self.metabolite_list = []
        for k in met_list:
            idx = k.find('.mlInfo')
            if idx > -1:
                self.metabolite_list.append(k[:idx])

        # end set_metabolite_list

    def set_metabolite_information(self, metabolite_name='', metabolite_information=''):
        if len(metabolite_name) == 0 or len(metabolite_information) == 0:
            return

        if metabolite_name not in self.nmr1d_data.keys():
            self.nmr1d_data[metabolite_name] = nmr1dData.Nmr1dData()
            self.nmr1d_data[metabolite_name].init_data(metabolite_information)

        # end set_metabolite_information

    def set_peak_information(self): # GUI item, not finished
        n_carbons = len(self.hsqc_data[self.cur_metabolite].c13_shifts)
        c13_idx = self.hsqc_data[self.cur_metabolite].h1_index[self.cur_peak - 1]
        self.hsqc_data[self.cur_metabolite].spin_systems[self.cur_peak - 1]['c13_shifts'] = []
        cs = []
        for l in range(self.hsqc_data[self.cur_metabolite].n_bonds + 1):
            cs.append(c13_idx - l)
            cs.append(c13_idx + l)

        used = []
        unique = [used.append(x) for x in cs if x not in used]
        used2 = np.array(used)
        used2 = used2[np.where(np.array(used) > 0)[0]]
        used2 = used2[np.where(used2 < n_carbons + 1)[0]].tolist()[1:]
        used2 = self.sub_lists(used2)
        for k in range(len(used2)):
            used2[k].insert(0, used[0])

        used2.sort(key=len)
        c13_idx = []
        j_cc = []
        c13_nc0 = []
        c13_shifts = []
        for k in range(len(used2)):
            if len(used2[k]) == 1:
                c13_idx.append(used2[k])
                j_cc.append([])
                c13_nc0.append([self.hsqc_data[self.cur_metabolite].c13_nc[k]])
                c13_shifts.append([self.hsqc_data[self.cur_metabolite].c13_shifts[self.hsqc_data[self.cur_metabolite].h1_index[self.cur_peak - 1] - 1]])
            else:
                c13_idx2 = []
                j_cc2 = []
                c13_nc2 = []
                c13_shifts2 = []
                for l in range(len(used2[k]) - 1):
                    if used2[k][0] > used2[k][l + 1]:
                        idx1 = np.where(self.hsqc_data[self.cur_metabolite].j_nuc1 == used2[k][l + 1])[0]
                        idx2 = np.where(self.hsqc_data[self.cur_metabolite].j_nuc2 == used2[k][0])[0]
                        ddd1 = self.hsqc_data[self.cur_metabolite].j_nuc1 == used2[k][l + 1]
                        ddd2 = self.hsqc_data[self.cur_metabolite].j_nuc2 == used2[k][0]
                    else:
                        idx1 = np.where(self.hsqc_data[self.cur_metabolite].j_nuc1 == used2[k][0])[0]
                        idx2 = np.where(self.hsqc_data[self.cur_metabolite].j_nuc2 == used2[k][l + 1])[0]
                        ddd1 = self.hsqc_data[self.cur_metabolite].j_nuc1 == used2[k][0]
                        ddd2 = self.hsqc_data[self.cur_metabolite].j_nuc2 == used2[k][l + 1]

                    if len(idx1) > 0 and len(idx2) > 0:
                        has_j = False
                        idx0 = 0
                        for m in range(len(ddd1)):
                            if ddd1[m] == True and ddd2[m] == True:
                                if abs(self.hsqc_data[self.cur_metabolite].j_nuc1[m] - self.hsqc_data[self.cur_metabolite].j_nuc2[m]) <= self.hsqc_data[self.cur_metabolite].n_bonds:
                                    has_j = True
                                    idx0 = min(m, len(self.hsqc_data[self.cur_metabolite].c13_nc) - 1)

                        if has_j:
                            if len(c13_idx2) == 0:
                                c13_idx2.append(used2[k][0])

                            if len(c13_nc2) == 0:
                                c13_nc2.append(self.hsqc_data[self.cur_metabolite].c13_nc[idx0])

                            if len(c13_shifts2) == 0:
                                c13_shifts2.append(self.hsqc_data[self.cur_metabolite].c13_shifts[self.hsqc_data[self.cur_metabolite].h1_index[self.cur_peak - 1] - 1])

                            c13_idx2.append(used2[k][l + 1])
                            c13_nc2.append(self.hsqc_data[self.cur_metabolite].c13_nc[l + 1])
                            j_cc2.append(self.hsqc_data[self.cur_metabolite].j_cc[idx0])
                            c13_shifts2.append(self.hsqc_data[self.cur_metabolite].c13_shifts[c13_idx2[l + 1] - 1])

                if len(c13_idx2) > 0 and len(j_cc2) > 0 and c13_idx2 not in c13_idx:
                    c13_idx.append(c13_idx2)
                    j_cc.append(j_cc2)
                    c13_nc0.append(c13_nc2)
                    c13_shifts.append(c13_shifts2)

        if len(self.hsqc_data[self.cur_metabolite].c13_picked[self.cur_peak - 1]) > 0:
            for k in range(len(c13_shifts)):
                c13_shifts[k][0] = np.mean(self.hsqc_data[self.cur_metabolite].c13_picked[self.cur_peak - 1])

        self.hsqc_data[self.cur_metabolite].spin_systems[self.cur_peak - 1]['c13_idx'] = c13_idx
        self.hsqc_data[self.cur_metabolite].spin_systems[self.cur_peak - 1]['j_cc'] = j_cc
        self.hsqc_data[self.cur_metabolite].spin_systems[self.cur_peak - 1]['c13_nc'] = c13_nc0
        self.hsqc_data[self.cur_metabolite].spin_systems[self.cur_peak - 1]['c13_shifts'] = c13_shifts
        c13_offset = []
        for k in range(len(c13_idx)):
            key = ' '.join(str(e) for e in c13_idx[k])
            if key in self.hsqc_data[self.cur_metabolite].c13_offset.keys():
                c13_offset.append(self.hsqc_data[self.cur_metabolite].c13_offset[key])
            else:
                c13_offset.append(0)

        self.hsqc_data[self.cur_metabolite].spin_systems[self.cur_peak - 1]['c13_offset'] = c13_offset
        if 'contribution' not in self.hsqc_data[self.cur_metabolite].spin_systems[self.cur_peak - 1].keys() or len(self.hsqc_data[self.cur_metabolite].spin_systems[self.cur_peak - 1]['c13_idx']) != len(self.hsqc_data[self.cur_metabolite].spin_systems[self.cur_peak - 1]['contribution']):
            contribution = []
            contribution.append(100)
            for k in range(len(c13_shifts) - 1):
                contribution.append(0)

            self.hsqc_data[self.cur_metabolite].spin_systems[self.cur_peak - 1]['contribution'] = contribution

        # end set_peak_information

    def sub_lists(self, l): # GUI item, not finished
        lists = [[]]
        for i in range(len(l) + 1):
            for j in range(i):
                lists.append(l[j: i])

        return lists
