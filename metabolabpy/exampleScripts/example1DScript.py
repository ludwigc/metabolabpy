msg = self.clear()
dataPath = 'Replace with directory containing Bruker data'		# select directory containing bruker data, interactive for file dialog
dataSets = []						# add comma separated list of experiment numbers (e.g. [1, 2, 3])
self.readSpcs(dataPath,dataSets)				# reading Bruker spectra
msg = self.nd.setZeroFill([131072])				# zero fill to 131072 data points
msg = self.nd.setLb([0.3])					# 0.3 Hz line broadening
msg = self.nd.setWindowType([1])				# exponential window function
msg = self.nd.ftAll()					# Fourier Transform all NMR spectra
msg = self.nd.autorefAll()					# automatically reference to TMSP
msg = self.nd.autophase1dAll()				# automatically phase correct NMR spectra (comment out if not needed with # in front)
msg = self.enableBaseline()					# this option is needed if the above automatic phase correction is carried out, comment out if automatic phase correction is not performed
self.nd.e = 0						# make first spectrum the current experiment
r = np.linspace(0,0.6,len(dataSets))				# define colour gradient for spectrum plotting
g = np.linspace(0,0,len(dataSets))				# define colour gradient for spectrum plotting
b = np.linspace(0.6,0,len(dataSets))				# define colour gradient for spectrum plotting
for k in range(len(dataSets)):					# set plot colour for each spectrum
    self.nd.nmrdat[0][k].disp.posColRGB = (r[k], g[k], b[k])

self.setPhRefExp(1)						# select reference spectrum for manual phase correction
msg = self.nd.pp.init(len(dataSets))				# initialise data pre-processing options
self.nd.pp.plotSelect = np.array([])				# select spectra to plot (e.g. [0, 1, 2])
self.nd.pp.classSelect = np.array([])				# associate classes with each spectrum (e.g. ["WT", "Mutation", "Mutation"])
self.nd.pp.excludeStart = np.array([])				# define regions to exlcude from statistical data analysis (start), e.g. [10.0, 4.5, -2.5]
self.nd.pp.excludeEnd   = np.array([])				# define regions to exlcude from statistical data analysis (end), e.g. [12.0, 5.5, 0.5]
