import numpy as np

def sat_rec(time_vector=[], intensity=0.0, rate=1.0):
    if len(time_vector) == 0:
        return []

    time_vector = np.array(time_vector)
    data_vector = np.copy(intensity*(1.0 - np.exp(-time_vector*rate)))
    return data_vector
    # end sat_rec
