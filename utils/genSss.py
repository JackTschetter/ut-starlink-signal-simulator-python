''' Ported from Matlab to Python by Jack R. Tschetter on 03/13/2025.
    The following key considerations were made.

    1. ) Handling the Frequency Domain Sequence (sssVecFull):
        MATLAB loads it from sssVec.mat, but in Python, we need to:
            Either load it from a .mat file using scipy.io.loadmat
            Or generate it if the expected pattern is known.
    2. ) IFFT for Time Domain Representation:
        MATLAB uses ifft(sssFreqDomainVec), which is directly mapped to numpy.fft.ifft.
    3. ) Cyclic Prefix Handling:
        MATLAB prepends the last Ng samples to the start of sssTimeDomainVec.
'''

import numpy as np
from scipy.io import loadmat

def gen_sss(sssVecFull=None):
    """
    Generates the Starlink downlink secondary synchronization sequence (SSS).

    Parameters:
    sssVecFull (ndarray, optional): N-by-1 SSS vector in the frequency domain.
                                    If None, it should be loaded externally.

    Returns:
    sssTimeDomainVec (ndarray): (N + Ng)-by-1 SSS in the time domain with cyclic prefix.
    sssFreqDomainVec (ndarray): N-by-1 SSS in the frequency domain.
    """
    # Number of subcarriers
    N = 1024
    # Duration of cyclic prefix (CP)
    Ng = 32

    # Load sssVecFull if not provided
    if sssVecFull is None:
        try:
            mat_data = loadmat("sssVec.mat")  # Adjust the path if necessary
            sssVecFull = mat_data["sssVecFull"].flatten()
        except FileNotFoundError:
            raise ValueError("sssVecFull must be provided or loaded from a file.")

    # Validation checks
    if len(sssVecFull) != N:
        raise ValueError("Frequency domain SSS expression must be of length N.")
    
    if not np.all(np.abs(sssVecFull[2:-2]) == 1):  # MATLAB indexing (3:end-2) -> Python (2:-2)
        raise ValueError("sssVecFull expected to have 4QAM symbols with unit coordinates.")

    # Frequency domain SSS
    sssFreqDomainVec = sssVecFull

    # Convert to Time Domain (IFFT)
    sssTimeDomainVec = (N / np.sqrt(N)) * np.fft.ifft(sssFreqDomainVec)

    # Add Cyclic Prefix (CP)
    sssTimeDomainVec = np.concatenate((sssTimeDomainVec[-Ng:], sssTimeDomainVec))

    return sssTimeDomainVec, sssFreqDomainVec
