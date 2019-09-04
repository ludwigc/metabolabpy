lds = len(self.nd.nmrdat[0])
msg = self.nd.pp.init(lds)					# initialise data pre-processing options
self.nd.pp.plotSelect = np.array([])				# select spectra to plot (e.g. [0, 1, 2])
self.nd.pp.classSelect = np.array([])				# associate classes with each spectrum (e.g. ["WT", "Mutation", "Mutation"])
self.nd.pp.excludeStart = np.array([])				# define regions to exlcude from statistical data analysis (start), e.g. [10.0, 4.5, -2.5]
self.nd.pp.excludeEnd   = np.array([])				# define regions to exlcude from statistical data analysis (end), e.g. [12.0, 5.5, 0.5]
