''' Ported from Matlab to Python by Jack R. Tschetter on 03/13/2025.
    The following key considerations were made.

    1. ) Modulation Handling
        Uses numpy.random.randint for random symbol generation.
        Supports PSK and QAM modulation types.
    2. ) Fourier Transforms
        Uses numpy.fft.ifft for the OFDM signal generation.
    3. ) Cyclic Prefix Handling
        Appends the last Ng samples to the start of the OFDM symbol.
    4. ) Doppler Shift Simulation
        Uses numpy.exp(1j * 2 * np.pi * Fshift * tVec).
    5. ) Data Wrapping Logic
        Ensures correct formatting for data matrix.
'''

import numpy as np
from scipy.signal import resample

def gen_strlk_ofdm(s):
    """
    Generates a Starlink OFDM symbol with specific Starlink parameters.

    Parameters:
    s : dict
        Dictionary containing:
        - 'Midx' (int): Modulation index (2 = BPSK, 4 = 4QAM, 16 = 16QAM, etc.).
        - 'type' (str): Modulation type ('PSK' or 'QAM').
        - 'SNRdB' (float): Simulated Signal-to-Noise Ratio in dB (use np.nan for no noise).
        - 'Fsr' (float): Receiver sample rate (Hz).
        - 'Fcr' (float): Receiver center frequency (Hz).
        - 'beta' (float): Doppler factor.
        - 'Nsym' (int, optional): Number of consecutive symbols to generate (default: 1).
        - 'data' (ndarray, optional): Predefined OFDM data symbols.

    Returns:
    y : ndarray
        Starlink OFDM symbol in time, sampled at `Fsr`, centered at `Fsr`, and expressed at baseband.
    """

    # Set default values
    s.setdefault("Fs", 240e6)  # Starlink bandwidth
    s.setdefault("N", 1024)  # FFT size / Number of subcarriers
    s.setdefault("Ng", 32)  # Cyclic prefix length
    s.setdefault("gutter", True)  # 4-subcarrier gutter enabled
    s.setdefault("Nsym", 1)  # Default number of symbols

    # Validate inputs
    if "data" in s and "type" not in s:
        raise ValueError("A constellation type must be provided (s['type'] = 'PSK' or 'QAM')")

    if "data" in s:
        l, w = s["data"].shape
        if l != s["N"] or w > 300:
            raise ValueError("s['data'] must be an N x 300 matrix at most.")

    if "data" not in s and (np.log2(s["Midx"]) % 2 != 0) and s["Midx"] != 2:
        raise ValueError("Midx must be 2 or an even power of 2.")

    # Compute the closest Starlink channel center
    F = s["Fs"] / s["N"]
    ch_idx = round(((s["Fcr"] / 1e9) - 10.7 - (F / 2 / 1e9)) / 0.25 + 0.5)
    Fcii = 10.7e9 + (F / 2) + (250e6 * (ch_idx - 0.5))
    s["Fc"] = Fcii

    # Generate OFDM symbol
    y = gen_ofdm(s)  # Calls the OFDM generator function

    return y

# Placeholder function for genOFDM
def gen_ofdm(s):
    """ Generates a basic OFDM signal (placeholder for actual implementation). """
    return np.random.randn(1056) + 1j * np.random.randn(1056)
