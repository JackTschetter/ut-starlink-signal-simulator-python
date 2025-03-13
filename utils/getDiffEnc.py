''' Ported from Matlab to Python by Jack R. Tschetter on 03/13/2025.
    np.angle() is used to compute atan2(imag(di), real(di)) in Matlab.
    The rounding logic is preserved.
    The NaN row is prepended using np.vstack.
'''

import numpy as np

def get_diff_enc(x):
    """
    getDiffEnc returns the differentially decoded sequence x.

    Parameters:
    x : ndarray
        N x K matrix representing a differential 4QAM encoded sequence.
        If only Ns < N elements of the sequence are known, x should still 
        be N long, with np.nan in place of unknown elements.

    Returns:
    s : ndarray
        N x 1 vector of differentially decoded elements. If an element 
        could not be differentially decoded, it takes a value of np.nan.
    """
    N, K = x.shape

    aim1 = x[:-1, :]
    ai = x[1:, :]
    di = ai * np.conj(aim1)
    thetai = np.angle(di)
    s = np.round(thetai * 2 / np.pi).astype(int)

    # Handling the equivalent phase values
    s[s == -2] = 2
    s[s == -1] = 3

    # Prepending NaN for the first row
    s = np.vstack((np.full((1, K), np.nan), s))

    return s
