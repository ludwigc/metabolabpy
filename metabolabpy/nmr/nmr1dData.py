"""
NMR 1D data

"""

import numpy as np
from metabolabpy.nmr import nmrConfig  # pragma: no cover
import os

class Nmr1dData:

    def __init__(self):
        self.name = ''
        self.alt_name = ''
        self.source = ''
        self.iupac = ''
        self.inchi_identifier = ''
        self.inchi_key = ''
        self.smiles = ''
        self.formula = ''
        self.mass = 0.0
        self.pka = ''
        self.hmdb = []
        self.smpdb = []
        self.kegg = []
        self.chebi = ''
        self.cid = ''
        self.h1_shifts = np.array([])
        self.h1_offsets = np.array([])
        self.j_hh = np.array([])
        self.j_nuc1 = np.array([])
        self.j_nuc2 = np.array([])
        self.h1_index = np.array([])
        self.h1_number = np.array([])
        self.sim_spc = np.array([[]])
        self.spin_system = []
        self.intensity = 0
        self.r2 = []
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
        str_str += " 1H chemical shifts: {}\n".format(self.h1_shifts)
        str_str += " 1H index          : {}\n".format(self.h1_index)
        str_str += "Number of 1H nuclei: {}\n".format(self.h1_number)
        str_str += "Jhh                : {}\n".format(self.j_hh)
        str_str += "1H intensity       : {}\n".format(self.intensity)
        str_str += "_________________________________________________________\n"
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
        idx1 = metabolite_information.find('IUPAC')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        self.iupac = mi[idx1 + 1:idx2].strip()
        idx1 = metabolite_information.find('InChi_Identifier')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        self.inchi_identifier = mi[idx1 + 1:idx2].strip()
        idx1 = metabolite_information.find('InChi_key')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        self.inchi_key = mi[idx1 + 1:idx2].strip()
        idx1 = metabolite_information.find('SMILES')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        self.smiles = mi[idx1 + 1:idx2].strip()
        idx1 = metabolite_information.find('Formula')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        self.formula = mi[idx1 + 1:idx2].strip()
        idx1 = metabolite_information.find('Mass')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        self.mass = float(mi[idx1 + 1:idx2].strip())
        idx1 = metabolite_information.find('pKa')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        self.pKa = mi[idx1 + 1:idx2].strip()
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
        idx1 = metabolite_information.find('CHEBI')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        self.chebi = mi[idx1 + 1:idx2].strip()
        idx1 = metabolite_information.find('CID')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        self.cid = mi[idx1 + 1:idx2].strip()
        idx1 = metabolite_information.find('H1chemicalShift')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        self.h1_shifts = np.array(mi[idx1 + 1:idx2].strip().split(), dtype='float64')
        self.r2 = []
        for k in range(len(self.h1_shifts)):
            self.r2.append(4.0)

        idx1 = metabolite_information.find('jHH')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        j_hh_info = mi[idx1 + 1:idx2].strip().split()
        hh = []
        nuc1 = []
        nuc2 = []
        for k in j_hh_info:
            hh.append(float(k.split(';')[2]))
            nuc1.append(int(k.split(';')[0]))
            nuc2.append(int(k.split(';')[1]))

        self.j_hh = np.array(hh)
        self.j_nuc1 = np.array(nuc1)
        self.j_nuc2 = np.array(nuc2)
        idx1 = metabolite_information.find('H1Index')
        mi = metabolite_information[idx1:]
        idx1 = mi.find(':')
        idx2 = mi.find('\n')
        h1idx = mi[idx1 + 1:idx2].strip().split()
        self.h1_index = []
        self.h1_number = []
        for k in h1idx:
            self.h1_index.append(int(k.split(';')[0]))
            self.h1_number.append(int(k.split(';')[1]))

        # end init_data

    def set_r2(self, new_r2=4.0):
        for k in range(len(self.r2)):
            self.r2[k] = new_r2

        # end set_r2

