kz = self.clear(True)
data_path = 'Replace with directory containing Bruker data'	# select directory containing bruker data, interactive for file dialog
data_sets = []										# add comma separated list of experiment numbers (e.g. [1, 2, 3])
self.read_spcs([data_path],data_sets)			# reading Bruker spectra
lds = len(self.nd.nmrdat[0])
for k in range(lds):									#
    self.nd.nmrdat[0][k].display.n_levels = 34					# setting number of contour levels
    self.nd.nmrdat[0][k].display.min_level = 0.01					# setting minimum contour level
    self.nd.nmrdat[0][k].display.max_level = 0.08					# setting maximum contour level
    self.nd.nmrdat[0][k].proc.invert_matrix[1] = True					# adjust as necessary

self.ft_all()
self.nd.s = 0										# make first data set current data set
self.nd.e = 0										# make first experiment current experiment
self.plot_spc(keep_zoom=kz)

