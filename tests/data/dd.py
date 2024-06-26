from metabolabpy.nmr import nmrDataSet
p_name = "nmrData"
e_name = "10"
nd = nmrDataSet.NmrDataSet()  # create nmrDataSet object
nd.read_spc(p_name, e_name)  # check if Bruker data can be read
nd.nmrdat[0][0].proc.ph0[0] = -85.94556
nd.nmrdat[0][0].proc.ph1[0] = 6.64201
nd.nmrdat[0][0].proc.lb[0] = 0.3
nd.ft()
nd.nmrdat[0][0].autobaseline1d()
nd.nmrdat[0][0].auto_ref()
e_name = "11"  # 1D NMR data in exp 1
nd.read_spc(p_name, e_name)  # check if Bruker data can be read
nd.nmrdat[0][1].proc.ph0[0] = -79.29513383337031
nd.nmrdat[0][1].proc.ph1[0] = 7.790498486770286
nd.nmrdat[0][1].proc.lb[0] = 0.3
nd.e = 1
nd.ft()
nd.pp.flag_scale_spectra = True
nd.pp.scale_spectra_ref_spc = 1
print(f'nd.nmrdat[0][1].spc[0].real.max(): {nd.nmrdat[0][1].spc[0].real.max()}')

nd.pp.scale_pqn = True
nd.data_pre_processing()
print(f'nd.nmrdat[0][1].spc[0].real.max(): {nd.nmrdat[0][1].spc[0].real.max()}')
nd.reset_data_pre_processing()
print(f'nd.nmrdat[0][1].spc[0].real.max(): {nd.nmrdat[0][1].spc[0].real.max()}')
nd.pp.scale_pqn = False
nd.pp.flag_scale_spectra = True
nd.pp.preserve_overall_scale = True
nd.data_pre_processing()
print(f'nd.nmrdat[0][1].spc[0].real.max(): {nd.nmrdat[0][1].spc[0].real.max()}')
nd.reset_data_pre_processing()
nd.pp.preserve_overall_scale = False
nd.data_pre_processing()
print(f'nd.nmrdat[0][1].spc[0].real.max(): {nd.nmrdat[0][1].spc[0].real.max()}')



