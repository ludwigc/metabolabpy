spectrum = self.nd.nmrdat[0][0] # Select set and experiment we would like to adjust 
spectrum.auto_add_baseline_points() # Automatically select baseline points for spline correction

new_baseline = spectrum.calc_spline_baseline() # Calculate the new baseline
spectrum.spc[0] = spectrum.spc[0] - new_baseline # Re-assign the initial spectrum
