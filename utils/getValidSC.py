''' Ported from Matlab to Python by Jack R. Tschetter on 03/13/2025.
    The following key considerations were made.

    1. ) Computes valid Starlink OFDM subcarrier indices
    2. ) Uses Starlink-specific subcarrier mapping (513-1024, 1-512)
    3. ) Handles edge cases when the frequency range extends beyond the channel
'''

import numpy as np

def get_valid_sc(Fcr, Fsr):
    """
    Returns the valid subcarrier indices based on the provided 
    center and sampling frequencies.

    Parameters:
    Fcr : float
        Receiver center frequency in Hz.
    Fsr : float
        Receiver sampling frequency in Hz.

    Returns:
    idxs : ndarray
        Array of valid subcarrier indices (1 to 1024).
    """
    
    Fs = 240e6  # Starlink signal bandwidth
    F = Fs / 1024  # Starlink OFDM subcarrier bandwidth

    # Compute the closest Starlink channel center
    ch_idx = round(((Fcr / 1e9) - 10.7 - (F / 2 / 1e9)) / 0.25 + 0.5)
    Fcii = (10.7e9 + F / 2 + 250e6 * (ch_idx - 0.5))

    # Compute subcarrier channel centers
    Fsc = np.arange(Fcii - round(Fs / 2 - F / 2), Fcii + round(Fs / 2 - F / 2) + F, F)

    # Compute capture frequency range
    Fcr_start = Fcr - round(Fsr / 2)
    Fcr_end = Fcr + round(Fsr / 2)

    # Find valid subcarrier index range
    Fsc_start_idx = np.where((Fsc - F / 2 - Fcr_start) > 0)[0]
    if Fsc_start_idx.size == 0:
        raise ValueError("Frequency range does not match: beginning of band is outside a Starlink channel.")
    Fsc_start_idx = Fsc_start_idx[0]

    Fsc_end_idx = np.where((Fsc + F / 2 - Fcr_end) > -1)[0]
    Fsc_end_idx = min(Fsc_end_idx[0], 1023) if Fsc_end_idx.size > 0 else 1023

    # Map indices from 1-1024 using Starlink indexing convention
    idxs = np.arange(Fsc_start_idx, Fsc_end_idx + 1)
    idxs = (idxs + 511) % 1024 + 1  # Wrap around Starlink OFDM indexing

    return idxs
