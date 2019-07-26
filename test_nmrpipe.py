import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as pl
import numpy as np
from metabolabpy.nmr import nmrData
import os
import time

spcName = 'test.dat'
dsNum   = '50'
dsName  = '/Users/ludwigc/dataDir/pyData'


#npd =  nmrpipeData.NmrPipeData()
#npd.readPipe(pName, fName)
nd = nmrData.NmrData()
nd.dataSetName   = dsName
nd.dataSetNumber = dsNum
nd.readSpc()
nd.readPipe2D(dsName + os.sep + dsNum + '.proc', spcName)

minLev                   = 0.01
maxLev                   = 0.08
nLev                     = 17
mm                       = nd.spc.real.max()
posLev                   = np.linspace( minLev*mm, maxLev*mm,nLev)
negLev                   = np.linspace(-maxLev*mm,-minLev*mm,nLev)
#pl.contour(nd.spc.real, posLev, colors = 'b', linestyles = 'solid')
#pl.contour(nd.spc.real, negLev, colors = 'r', linestyles = 'solid')
#pl.show()
#pl.plot(spc)
#spc2 = nd.spc
#pl.subplot(2,1,1)
#pl.contour(nd.ppm1, nd.ppm2, nd.spc.real, posLev, colors = 'b', linestyles = 'solid')
#pl.contour(nd.ppm1, nd.ppm2, nd.spc.real, negLev, colors = 'r', linestyles = 'solid')
#pl.gca().invert_xaxis()
#pl.gca().invert_yaxis()
#spc0 = nd.spc
print("Start!")
t = time.time()
nd.phase2a(-70, 0, 0)
elapsed1 = time.time() - t
print("time elapsed: {:4.2f}".format(elapsed1))
#print("Start!")
#t = time.time()
#spc1 = nd.phase2a(spc0,-70, 0, 0)
#elapsed2 = time.time() - t
#print("time elapsed: {:4.2f}".format(elapsed2))
#nd.spc = nd.phase2a(nd.spc, -70, 0, 0)
#pl.subplot(2,1,1)
pl.contour(nd.ppm1, nd.ppm2, nd.spc.real, posLev, colors = 'b', linestyles = 'solid')
pl.contour(nd.ppm1, nd.ppm2, nd.spc.real, negLev, colors = 'r', linestyles = 'solid')
pl.gca().invert_xaxis()
pl.gca().invert_yaxis()
#pl.xlim(1.5344, 1.2765)
#pl.ylim(25.525, 15.998)
#pl.subplot(2,1,2)
#pl.contour(nd.ppm1, nd.ppm2, spc1.real, posLev, colors = 'b', linestyles = 'solid')
#pl.contour(nd.ppm1, nd.ppm2, spc1.real, negLev, colors = 'r', linestyles = 'solid')
#pl.gca().invert_xaxis()
#pl.gca().invert_yaxis()
#pl.xlim(1.5344, 1.2765)
#pl.ylim(25.525, 15.998)
pl.show()
