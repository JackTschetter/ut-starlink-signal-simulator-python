''' Ported from Matlab to Python by Jack R. Tschetter on 03/13/2025.
    The following key considerations were made.
    MATLAB uses 1-based indexing, while Python uses 0-based indexing, so the ciVec adjustments were made accordingly.
    MATLAB's flipud(lfsrSeq) is np.flipud(lfsr_seq) in Python.
    The cumsum function in NumPy is used to compute the phase shift for pk_seg_vec.
    The sequence is expanded using np.tile() to replicate it NpssSeg - 1 times.
'''

import numpy as np

def gen_pss():
    """
    Generates the Starlink downlink primary synchronization sequence (PSS).
    
    Returns:
        pss (ndarray): (N + Ng) x 1 PSS sequence, sampled at 240 MHz.
    """
    # Parameters
    Ng = 32  # Cyclic prefix length
    NpssSeg = 8  # Number of sub-segments in the PSS
    n = 7
    m = 2**n - 1

    # Fibonacci LFSR initialization
    ci_vec = np.array([3, 7]) - 1  # Convert to zero-based indexing
    a0_vec = np.array([0, 0, 1, 1, 0, 1, 0], dtype=int)
    
    # Generate LFSR sequence
    lfsr_seq = np.zeros(m, dtype=int)
    for idx in range(m):
        buffer = a0_vec[ci_vec]
        val = np.mod(np.sum(buffer), 2)  # XOR operation
        a0_vec = np.insert(a0_vec[:-1], 0, val)
        lfsr_seq[idx] = val

    # Process sequence
    lfsr_seq_mod = np.flipud(lfsr_seq)
    lfsr_seq_mod = np.insert(lfsr_seq_mod, 0, 0)  # Append leading zero
    seq = 2 * lfsr_seq_mod - 1  # Map {0,1} -> {-1,1}

    # Generate synchronization sequence
    pk_seg_vec = np.exp(-1j * np.pi / 4 - 1j * 0.5 * np.pi * np.cumsum(seq))
    pk_vec = np.concatenate([-pk_seg_vec, np.tile(pk_seg_vec, NpssSeg - 1)])

    # Add cyclic prefix
    pk_cp = pk_vec[-Ng:]
    pss = np.concatenate([-pk_cp, pk_vec])

    return pss
