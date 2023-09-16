import numpy as np
import math
from numba import njit

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

    mat -= np.mean(np.array([np.mean(mat[512:1024].real), np.mean(mat[-1024:-512].real)]))
    return mat
    # end phase3
