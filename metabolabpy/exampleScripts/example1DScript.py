msg = self.clear()
data_path = 'Replace with directory containing Bruker data'		# select directory containing bruker data, interactive for file dialog
data_sets = []								# add comma separated list of experiment numbers (e.g. [1, 2, 3])
self.read_spcs([data_path],data_sets)				# reading Bruker spectra
msg = self.nd.set_zero_fill([131072])				# zero fill to 131072 data points
msg = self.nd.set_lb([0.3])						# 0.3 Hz line broadening
msg = self.nd.set_window_type(['Exponential'])			# window function ('None', 'Exponential', 'Gaussian', 'Sine', 'QSine', 'SEM')
msg = self.nd.set_water_suppression('None')				# set water suppression ('None', 'Conv', 'Poly', 'WaveWat')
msg = self.nd.set_autobaseline(False)				# perform automatic baseline correction (False, True) 
msg = self.nd.ft_all()							# Fourier Transform all NMR spectra
msg = self.nd.auto_ref_all()						# automatically reference to TMSP
self.nd.e = 0						 		# make first spectrum the current experiment
self.set_ph_ref_exp(1)							# select reference spectrum for manual phase correction
zoom_checked = self.w.keepZoom.isChecked()
self.w.keepZoom.setChecked(False)
self.plot_spc()
self.w.keepZoom.setChecked(zoom_checked)
