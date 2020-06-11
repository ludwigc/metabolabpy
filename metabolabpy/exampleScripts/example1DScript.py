msg = self.clear()
dataPath = 'Replace with directory containing Bruker data'		# select directory containing bruker data, interactive for file dialog
dataSets = []						# add comma separated list of experiment numbers (e.g. [1, 2, 3])
self.readSpcs([dataPath],dataSets)				# reading Bruker spectra
msg = self.nd.setZeroFill([131072])				# zero fill to 131072 data points
msg = self.nd.setLb([0.3])					# 0.3 Hz line broadening
msg = self.nd.setWindowType([1])				# exponential window function
msg = self.nd.ftAll()					# Fourier Transform all NMR spectra
msg = self.nd.autorefAll()					# automatically reference to TMSP
self.nd.e = 0						# make first spectrum the current experiment
lds = len(self.nd.nmrdat[0])
r = np.linspace(0,0.6,lds)					# define colour gradient for spectrum plotting
g = np.linspace(0,0,lds)					# define colour gradient for spectrum plotting
b = np.linspace(0.6,0,lds)					# define colour gradient for spectrum plotting
for k in range(lds):						# set plot colour for each spectrum
    self.nd.nmrdat[0][k].disp.posColRGB = (r[k], g[k], b[k])

self.setPhRefExp(1)						# select reference spectrum for manual phase correction
self.plotSpc()
