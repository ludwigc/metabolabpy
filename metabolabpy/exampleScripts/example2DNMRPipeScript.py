msg = self.clear()
dataPath = 'Replace with directory containing Bruker and NMRPipe processed data'		# select directory containing bruker data, interactive for file dialog
dataSets = []						# add comma separated list of experiment numbers (e.g. [1, 2, 3])
procDataName = 'test.dat'				# name of NMRPipe processed data file (resides in e.g. 1.proc)
self.readNMRPipeSpcs(dataPath,dataSets, procDataName)				# reading Bruker spectra
lds = len(self.nd.nmrdat[0])
for k in range(lds):						#
    self.nd.nmrdat[0][k].disp.nLevels = 34			# setting number of contour levels
    self.nd.nmrdat[0][k].disp.minLevel = 0.01			# setting minimum contour level
    self.nd.nmrdat[0][k].disp.maxLevel = 0.08			# setting maximum contour level

self.nd.s = 0						# make first data set current data set
self.nd.e = 0						# make first experiment current experiment
self.plotSpc()
