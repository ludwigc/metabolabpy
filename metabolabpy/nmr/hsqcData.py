"""
NMR HSQC data

"""

import numpy as np
from metabolabpy.nmr import nmrConfig  # pragma: no cover
import os


class HsqcData:

    def __init__(self):
        self.name = ''
        self.alt_name = ''
        self.source = ''
        self.hmdb = []
        self.smpdb = []
        self.kegg = []
        self.c13_shifts = np.array([])
        self.h1_shifts = np.array([])
        self.c13_picked = np.array([[]])
        self.h1_picked = np.array([[]])
        self.c13_picked_lib = np.array([[]])
        self.h1_picked_lib = np.array([[]])
        self.hsqc = np.array([])
        self.co_hsqc = np.array([])
        self.j_cc = np.array([])
        self.j_nuc1 = np.array([])
        self.j_nuc2 = np.array([])
        self.c13_intensities = np.array([])
        self.h1_index = np.array([])
        self.h1_number = np.array([])
        self.h1_suffix = np.array([])
        self.c13_index = np.array([])
        # end __init__

    def __str__(self):  # pragma: no cover
        str_str = "NMR HSQC data\n"
        str_str += "_________________________________________________________\n"
        str_str += "Name               : {}\n".format(self.name)
        str_str += "AltName            : {}\n".format(self.alt_name)
        str_str += "Source             : {}\n".format(self.source)
        str_str += "HMDB               : {}\n".format(self.hmdb)
        str_str += "SMPDB              : {}\n".format(self.smpdb)
        str_str += "KEGG               : {}\n".format(self.kegg)
        str_str += "13C chemical shifts: {}\n".format(self.c13_shifts)
        str_str += "13C picked         : {}\n".format(self.c13_picked)
        str_str += "13C index          : {}\n".format(self.c13_index)
        str_str += " 1H chemical shifts: {}\n".format(self.h1_shifts)
        str_str += " 1H picked         : {}\n".format(self.h1_picked)
        str_str += " 1H index          : {}\n".format(self.h1_index)
        str_str += "Number of 1H nuclei: {}\n".format(self.h1_number)
        str_str += " 1H suffix         : {}\n".format(self.h1_suffix)
        str_str += "HSQC active        : {}\n".format(self.hsqc)
        str_str += "CO-HSQC active     : {}\n".format(self.co_hsqc)
        str_str += "Jcc                : {}\n".format(self.j_cc)
        str_str += "Jnuc1              : {}\n".format(self.j_nuc1)
        str_str += "Jnuc2              : {}\n".format(self.j_nuc2)
        str_str += "13C intensities    : {}\n".format(self.c13_intensities)
        return str_str
        # end __str__

    def init_data(self, metabolite_information=''):
        if metabolite_information == None:
            return

        if len(metabolite_information) == 0:
            return

        #print('Init data')
        idx1 = metabolite_information.find('Name')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        self.name = mi[idx1 + 1:idx2].strip()
        idx1 = metabolite_information.find('AltName')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        self.alt_name = mi[idx1 + 1:idx2].strip()
        idx1 = metabolite_information.find('Source')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        self.source = mi[idx1 + 1:idx2].strip()
        idx1 = metabolite_information.find('HMDB ')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        self.hmdb = mi[idx1 + 1:idx2].strip().split()
        idx1 = metabolite_information.find('SMPDB ')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        self.smpdb = mi[idx1 + 1:idx2].strip().split()
        idx1 = metabolite_information.find('KEGG ')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        self.kegg = mi[idx1 + 1:idx2].strip().split()
        idx1 = metabolite_information.find('C13chemicalShift')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        self.c13_shifts = np.matrix(mi[idx1 + 1:idx2].strip()).A[0]
        for k in range(len(self.c13_shifts)):
            self.c13_picked = np.copy(np.append(self.c13_picked, [np.array(np.array([]))]))

        idx1 = metabolite_information.find('H1chemicalShift')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        self.h1_shifts = np.matrix(mi[idx1 + 1:idx2].strip()).A[0]
        for k in range(len(self.h1_shifts)):
            self.h1_picked = np.copy(np.append(self.h1_picked, [np.array(np.array([]))]))

        idx1 = metabolite_information.find('C13Intensities')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        self.c13_intensities = np.matrix(mi[idx1 + 1:idx2].strip()).A[0]
        idx1 = metabolite_information.find('C13HSQC')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        self.h1_picked = [[] for n in range(len(self.h1_shifts))]
        self.c13_picked = [[] for n in range(len(self.h1_shifts))]
        self.hsqc = np.matrix(mi[idx1 + 1:idx2].strip(), dtype=int).A[0]
        idx1 = metabolite_information.find('CO_HSQC')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        self.co_hsqc = np.matrix(mi[idx1 + 1:idx2].strip(), dtype=int).A[0]
        idx1 = metabolite_information.find('jCC')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        j_cc_info = mi[idx1 + 1:idx2].strip().split()
        cc = []
        nuc1 = []
        nuc2 = []
        for k in j_cc_info:
            cc.append(float(k.split(';')[2]))
            nuc1.append(int(k.split(';')[0]))
            nuc2.append(int(k.split(';')[1]))

        self.j_cc = np.array(cc)
        self.j_nuc1 = np.array(nuc1)
        self.j_nuc2 = np.array(nuc2)
        idx1 = metabolite_information.find('C13Index')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        c13idx = mi[idx1 + 1:idx2].strip().split()
        self.c13_index = []
        for k in c13idx:
            self.c13_index.append(int(k.split(';')[2]))

        idx1 = metabolite_information.find('H1Index')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        h1idx = mi[idx1 + 1:idx2].strip().split()
        self.h1_index = []
        self.h1_suffix = []
        self.h1_number = []
        for k in h1idx:
            self.h1_index.append(int(k.split(';')[0]))
            self.h1_number.append(int(k.split(';')[1]))
            suffix = chr(int(k.split(';')[2]) + 96)
            suffix = suffix.replace('`', '')
            suffix = suffix.replace('l', 'ab')
            suffix = suffix.replace('Û', 'abc')
            self.h1_suffix.append(suffix)

        # end init_data

