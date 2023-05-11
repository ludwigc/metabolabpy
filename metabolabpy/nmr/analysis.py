import numpy as np

def sat_rec(time_vector=[], intensity=0.0, rate=1.0):
    if len(time_vector) == 0:
        return []

    time_vector = np.array(time_vector)
    data_vector = np.copy(intensity*(1.0 - np.exp(-time_vector*rate)))
    return data_vector
    # end sat_rec

def sat_rec2(time_vector=[], parameters=[1.0, 1.0]):
    if len(time_vector) == 0:
        return []

    intensity = parameters[0]
    rate = parameters[1]
    time_vector = np.array(time_vector)
    data_vector = np.copy(intensity*(1.0 - np.exp(-time_vector*rate)))
    return data_vector
    # end sat_rec

def fct_sat_rec(fit_parameters=[], time_data=[], intensity_data=[]):
    if len(fit_parameters) != 2 or len(time_data) == 0 or len(intensity_data) == 0 or len(time_data) != len(intensity_data):
        return []

    time_data = np.array(time_data)
    intensity_data = np.array(intensity_data)
    intensity = fit_parameters[0]
    rate = fit_parameters[1]
    sim_data = sat_rec(time_data, intensity, rate)
    err = (sim_data - intensity_data)**2
    return err.sum()
    # end fct_sat_rec

def exp_decay(time_vector=[], intensity=0.0, rate=1.0, offset=0.0):
    if len(time_vector) == 0:
        return []

    time_vector = np.array(time_vector)
    data_vector = np.copy(intensity*np.exp(-time_vector*rate) + offset)
    return data_vector
    # end exp_decay

def fct_exp_decay(fit_parameters=[], time_data=[], intensity_data=[]):
    if len(fit_parameters) != 3 or len(time_data) == 0 or len(intensity_data) == 0 or len(time_data) != len(intensity_data):
        return []

    time_data = np.array(time_data)
    intensity_data = np.array(intensity_data)
    intensity = fit_parameters[0]
    rate = fit_parameters[1]
    offset = fit_parameters[2]
    sim_data = exp_decay(time_data, intensity, rate, offset)
    err = (sim_data - intensity_data)**2
    return err.sum()
    # end fct_exp_decay


