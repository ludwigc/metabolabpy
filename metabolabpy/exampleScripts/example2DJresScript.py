msg = self.clear()
dataPath = 'Replace with directory containing Bruker data'		# select directory containing bruker data, interactive for file dialog
dataSets = []						# add comma separated list of experiment numbers (e.g. [1, 2, 3])
self.readSpcs(dataPath,dataSets)				# reading Bruker spectra
lds = len(self.nd.nmrdat[0])
for k in range(lds):						#
    self.nd.nmrdat[0][k].disp.yLabel = '1H'			# setting axis label
    self.nd.nmrdat[0][k].disp.axisType2 = 'Hz'			# setting axis units for indirect dimension to [Hz]
    self.nd.nmrdat[0][k].proc.windowType = np.array([5, 3, 0])	# setting window functions (SEM for direct dimension, SINE for indirect dimension)
    self.nd.nmrdat[0][k].proc.lb[0] = 0.5				# setting line broadening for SEM
    self.nd.nmrdat[0][k].disp.nLevels = 34			# setting number of contour levels
    self.nd.nmrdat[0][k].disp.minLevel = 0.001			# setting minimum contour level
    self.nd.nmrdat[0][k].disp.maxLevel = 0.01			# setting maximum contour level

self.nd.s = 0						# make first data set current data set
self.nd.e = 0						# make first experiment current experiment
r = np.linspace(0,0.6,lds)					# define colour gradient for spectrum plotting
g = np.linspace(0,0,lds)					# define colour gradient for spectrum plotting
b = np.linspace(0.6,0,lds)					# define colour gradient for spectrum plotting
for k in range(lds):						# set plot colour for each spectrum
    self.nd.nmrdat[0][k].disp.posColRGB = (r[k], g[k], b[k])

msg = self.nd.ftAll()					# Fourier Transform all NMR spectra
self.nd.pjres(2,'skyline')					# calculate skyline projection in data set 2
self.nd.s = 1						# change to data set 2
self.nd.e = 0						# change to experiment 1 in data set 2
self.plotSpc()
