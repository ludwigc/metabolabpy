"""
NMRPipe processed spectrum class
"""

import os
import numpy as np
import re


class NmrPipeData:

    def __init__(self):
        self.spc = np.array([[]], dtype='float64')
        self.header = np.array([], dtype='float64')
        self.fdf2label = 0.0
        self.fdf2apod = 0.0
        self.fdf2sw = 0.0
        self.fdf2obs = 0.0
        self.fdf2orig = 0.0
        self.fdf2units = 0.0
        self.fdf2quadflag = 0.0
        self.fdf2ftflag = 0.0
        self.fdf2aqsign = 0.0
        self.fdf2lb = 0.0
        self.fdf2car = 0.0
        self.fdf2center = 0.0
        self.fdf2offppm = 0.0
        self.fdf2p0 = 0.0
        self.fdf2p1 = 0.0
        self.fdf2apodcode = 0.0
        self.fdf2apodq1 = 0.0
        self.fdf2apodq2 = 0.0
        self.fdf2apodq3 = 0.0
        self.fdf2c1 = 0.0
        self.fdf2zf = 0.0
        self.fdf2x1 = 0.0
        self.fdf2xn = 0.0
        self.fdf2ftsize = 0.0
        self.fdf2tdsize = 0.0
        self.fdf1label = 0.0
        self.fdf1apod = 0.0
        self.fdf1sw = 0.0
        self.fdf1obs = 0.0
        self.fdf1orig = 0.0
        self.fdf1units = 0.0
        self.fdf1ftflag = 0.0
        self.fdf1aqsign = 0.0
        self.fdf1lb = 0.0
        self.fdf1quadflag = 0.0
        self.fdf1car = 0.0
        self.fdf1center = 0.0
        self.fdf1offppm = 0.0
        self.fdf1p0 = 0.0
        self.fdf1p1 = 0.0
        self.fdf1apodcode = 0.0
        self.fdf1apodq1 = 0.0
        self.fdf1apodq2 = 0.0
        self.fdf1apodq3 = 0.0
        self.fdf1c1 = 0.0
        self.fdf1zf = 0.0
        self.fdf1x1 = 0.0
        self.fdf1xn = 0.0
        self.fdf1ftsize = 0.0
        self.fdf1tdsize = 0.0
        self.fdf3label = 0.0
        self.fdf3apod = 0.0
        self.fdf3obs = 0.0
        self.fdf3sw = 0.0
        self.fdf3orig = 0.0
        self.fdf3ftflag = 0.0
        self.fdf3aqsign = 0.0
        self.fdf3size = 0.0
        self.fdf3quadflag = 0.0
        self.fdf3units = 0.0
        self.fdf3p0 = 0.0
        self.fdf3p1 = 0.0
        self.fdf3car = 0.0
        self.fdf3center = 0.0
        self.fdf3offppm = 0.0
        self.fdf3apodcode = 0.0
        self.fdf3apodq1 = 0.0
        self.fdf3apodq2 = 0.0
        self.fdf3apodq3 = 0.0
        self.fdf3c1 = 0.0
        self.fdf3zf = 0.0
        self.fdf3x1 = 0.0
        self.fdf3xn = 0.0
        self.fdf3ftsize = 0.0
        self.fdf3tdsize = 0.0
        self.fdmagic = 0.0
        self.fdfltformat = 0.0
        self.fdfltorder = 0.0
        self.fdsize = 0.0
        self.fdrealsize = 0.0
        self.fdspecnum = 0.0
        self.fdquadflag = 0.0
        self.fd2dphase = 0.0
        self.fdtemperature = 0.0
        self.strip = np.array([0.0, 0.0])
        self.full_size = 0.0
        self.strip_size = 0.0
        self.x_size = 1.0
        self.y_size = 1.0
        self.z_size = 1.0
        # end __init__

    def read_header(self):
        self.fdf2label = self.header[16]
        self.fdf2apod = self.header[95]
        self.fdf2sw = self.header[100]
        self.fdf2obs = self.header[119]
        self.fdf2orig = self.header[101]
        self.fdf2units = self.header[152]
        self.fdf2quadflag = self.header[56]
        self.fdf2ftflag = self.header[220]
        self.fdf2aqsign = self.header[64]
        self.fdf2lb = self.header[111]
        self.fdf2car = self.header[66]
        self.fdf2center = self.header[79]
        self.fdf2offppm = self.header[480]
        self.fdf2p0 = self.header[109]
        self.fdf2p1 = self.header[110]
        self.fdf2apodcode = self.header[413]
        self.fdf2apodq1 = self.header[415]
        self.fdf2apodq2 = self.header[416]
        self.fdf2apodq3 = self.header[417]
        self.fdf2c1 = self.header[418]
        self.fdf2zf = self.header[108]
        self.fdf2x1 = self.header[257]
        self.fdf2xn = self.header[258]
        self.fdf2ftsize = self.header[96]
        self.fdf2tdsize = self.header[386]
        self.fdf1label = self.header[18]
        self.fdf1apod = self.header[428]
        self.fdf1sw = self.header[229]
        self.fdf1obs = self.header[218]
        self.fdf1orig = self.header[249]
        self.fdf1units = self.header[234]
        self.fdf1ftflag = self.header[222]
        self.fdf1aqsign = self.header[475]
        self.fdf1lb = self.header[243]
        self.fdf1quadflag = self.header[55]
        self.fdf1car = self.header[67]
        self.fdf1center = self.header[80]
        self.fdf1offppm = self.header[481]
        self.fdf1p0 = self.header[245]
        self.fdf1p1 = self.header[246]
        self.fdf1apodcode = self.header[414]
        self.fdf1apodq1 = self.header[420]
        self.fdf1apodq2 = self.header[421]
        self.fdf1apodq3 = self.header[422]
        self.fdf1c1 = self.header[423]
        self.fdf1zf = self.header[437]
        self.fdf1x1 = self.header[259]
        self.fdf1xn = self.header[260]
        self.fdf1ftsize = self.header[98]
        self.fdf1tdsize = self.header[387]
        self.fdf3label = self.header[20]
        self.fdf3apod = self.header[50]
        self.fdf3obs = self.header[10]
        self.fdf3sw = self.header[11]
        self.fdf3orig = self.header[12]
        self.fdf3ftflag = self.header[13]
        self.fdf3aqsign = self.header[476]
        self.fdf3size = self.header[15]
        self.fdf3quadflag = self.header[51]
        self.fdf3units = self.header[58]
        self.fdf3p0 = self.header[60]
        self.fdf3p1 = self.header[61]
        self.fdf3car = self.header[68]
        self.fdf3center = self.header[81]
        self.fdf3offppm = self.header[482]
        self.fdf3apodcode = self.header[400]
        self.fdf3apodq1 = self.header[401]
        self.fdf3apodq2 = self.header[402]
        self.fdf3apodq3 = self.header[403]
        self.fdf3c1 = self.header[404]
        self.fdf3zf = self.header[438]
        self.fdf3x1 = self.header[261]
        self.fdf3xn = self.header[262]
        self.fdf3ftsize = self.header[200]
        self.fdf3tdsize = self.header[388]
        self.fdmagic = self.header[0]
        self.fdfltformat = self.header[1]
        self.fdfltorder = self.header[2]
        self.fdsize = self.header[99]
        self.fdrealsize = self.header[97]
        self.fdspecnum = self.header[219]
        self.fdquadflag = self.header[106]
        self.fd2dphase = self.header[256]
        self.fdtemperature = self.header[157]
        if (self.fdf2x1 == 0.0) & (self.fdf2xn == 0.0):
            self.strip_size = self.fdf2ftsize
        else:
            self.strip = np.array([self.fdf2x1, self.fdf2xn])
            self.full_size = self.fdf2ftsize
            self.strip_size = abs(self.fdf2x1 - self.fdf2xn) + 1.0

        self.x_size = self.fdf1ftsize
        self.y_size = self.strip_size
        # end read_header

    def read_pipe(self, p_name, f_name):
        spc3d = re.compile('.+\d+.+')
        f_name2 = p_name + os.sep + f_name
        f = open(f_name2, 'rb')
        self.header = np.resize(self.header, 512)
        self.header = np.fromfile(f, dtype=np.float32, count=512)
        self.read_header()
        if len(spc3d.findall(f_name)) > 0:
            self.x_size = self.fdspecnum
            
        self.spc = np.resize(self.spc, (int(self.x_size), int(self.y_size)))
        spc = np.array([[]], dtype='float64')
        spc = np.fromfile(f, dtype=np.float32, count=int(self.x_size * self.y_size))
        spc = spc.reshape(int(self.x_size), int(self.y_size))
        self.spc = spc
        f.close()
        # end read_pipe
