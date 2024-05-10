import numpy as np
import math
from numba import njit

@njit
def phase3a(mat, ph0, ph1, npts):
    ph0_1 = -ph0 * math.pi / 180.0
    ph1_1 = -ph1 * math.pi / 180.0
    for k in range(int(npts)):
        frac = float(k) / float(npts)
        ph = ph0_1 + frac * ph1_1
        t = complex(math.cos(ph) * mat[k].real + math.sin(ph) * mat[k].imag,
                    -math.sin(ph) * mat[k].real + math.cos(ph) * mat[k].imag)
        mat[k] = t

    #mat -= np.mean(np.array([np.mean(mat[512:1024].real), np.mean(mat[-1024:-512].real)]))
    return mat
    # end phase3

@njit
def phase3(mat, ph0, ph1, npts):
    ph0_1 = -ph0 * math.pi / 180.0
    ph1_1 = -ph1 * math.pi / 180.0
    for k in range(int(npts)):
        frac = float(k) / float(npts)
        ph = ph0_1 + frac * ph1_1
        t = complex(math.cos(ph) * mat[k].real + math.sin(ph) * mat[k].imag,
                    -math.sin(ph) * mat[k].real + math.cos(ph) * mat[k].imag)
        mat[k] = t

    #mat -= np.mean(np.array([np.mean(mat[512:1024].real), np.mean(mat[-1024:-512].real)]))
    return mat
    # end phase3

@njit
def objective_function(phase, spc, start_peak, end_peak):
    ph0 = phase[0] * math.pi / 180.0
    ph1 = phase[1] * math.pi / 180.0
    npts = len(spc)
    phase_row = ph0 + np.linspace(0, npts - 1, npts) * ph1 / (npts - 1)
    spc2 = np.zeros(len(spc), dtype=complex)
    spc2.real = np.copy(spc.real * np.cos(phase_row) - spc.imag * np.sin(phase_row))
    spc2.imag = np.copy(spc.real * np.sin(phase_row) + spc.imag * np.cos(phase_row))
    of = np.zeros(len(start_peak))
    max_components = 7
    ddiff = int((max_components - 1) / 2)
    for k in range(len(of)):
        if end_peak[k] - start_peak[k] > max_components:
            of[k] = np.abs(np.sum(spc2[start_peak[k]:start_peak[k] + ddiff].real) - np.sum(
                spc2[end_peak[k] - ddiff:end_peak[k]].real))

    return np.sum(of)
    # end objective_function

@njit
def penalty_function(phase, start_peak, end_peak, spc):
    ph0 = phase[0] * math.pi / 180.0
    ph1 = phase[1] * math.pi / 180.0
    npts = len(spc)
    phase_row = ph0 + np.linspace(0, npts - 1, npts) * ph1 / (npts - 1)
    spc2 = np.zeros(npts, dtype=complex)
    spc2.real = np.copy(spc.real * np.cos(phase_row) - spc.imag * np.sin(phase_row))
    spc2.imag = np.copy(spc.real * np.sin(phase_row) + spc.imag * np.cos(phase_row))
    pf = 0.0
    for k in range(len(start_peak)):
        npts2 = end_peak[k] - start_peak[k]
        baseline = np.linspace(spc2[start_peak[k]].real, spc2[end_peak[k]].real, npts2)
        spc2[start_peak[k]:end_peak[k]].real -= baseline
        pf += np.sum((spc2[start_peak[k]:end_peak[k]].real - np.abs(spc2[start_peak[k]:end_peak[k]].real)) ** 2)

    return pf
    # end penalty_function

