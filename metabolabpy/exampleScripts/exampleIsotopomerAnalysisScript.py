#!/usr/bin/env python3
import numpy as np

from metabolabpytools import isotopomerAnalysis
ia = isotopomerAnalysis.IsotopomerAnalysis()
ia.read_hsqc_multiplets('/Users/ludwigc/Documents/hsqcData.xlsx')
ia.read_gcms_data('/Users/ludwigc/Documents/gcmsData.xlsx')
ia.read_nmr1d_data('/Users/ludwigc/Documents/nmr1dData.xlsx')
isotopomers = {}
isotopomers['L-LacticAcid'] = [[0, 0, 1], [0, 1, 1]]
isotopomers['L-Alanine'] = [[0, 0, 1], [0, 1, 1]]
for k in range(len(ia.metabolites)):
	if ia.metabolites[k] in isotopomers.keys():
		ia.fit_all_exps(metabolite=ia.metabolites[k], fit_isotopomers=isotopomers[ia.metabolites[k]])
	else:
		print(f'{ia.metabolites[k]} not in isotopomers list')


# implement as example script in metabolabpy
# write simulated distribution back to metabolbapy
# and simulate for printing
# add checkbox to print any simulated data