''' Ported from Matlab to Python by Jack R. Tschetter on 03/13/2025.
    The following key considerations were made.
    Struct Handling → Use a Python dictionary (dict) for input.
    Random Symbol Generation → Use numpy.random.randint.
    Modulation → Use scipy.signal for PSK and QAM modulation.
    Fourier Transforms → Use numpy.fft.ifft and numpy.fft.fftshift.
    Cyclic Prefix Handling → Use NumPy slicing.
    Doppler Shift Simulation → Use numpy.exp with complex exponentials.
    AWGN Channel Simulation → Use numpy.random.randn.
    Resampling → Use scipy.signal.resample.
'''

import numpy as np
from scipy.signal import resample
from scipy.special import erfc

def gen_ofdm(s):
    """
    Generates an OFDM signal with optional Doppler shift and noise.

    Parameters:
    s : dict
        Dictionary containing:
        - 'Fs' (float): Channel bandwidth (Hz).
        - 'N' (int): Number of subcarriers.
        - 'Ng' (int): Number of guard subcarriers (cyclic prefix).
        - 'Midx' (int): Modulation order (2 = BPSK, 4 = QPSK, 16 = 16QAM, etc.).
        - 'type' (str): Modulation type ('PSK' or 'QAM').
        - 'SNRdB' (float): Signal-to-noise ratio in dB. Set to np.nan for no noise.
        - 'Fsr' (float): Receiver sample rate (Hz).
        - 'Fcr' (float): Receiver center frequency (Hz).
        - 'Fc' (float): OFDM signal center frequency (Hz).
        - 'beta' (float): Doppler factor. Doppler shift is FD = -beta * Fc.
        - 'gutter' (bool, optional): Enables a 4-subcarrier gap at the center.
        - 'Nsym' (int, optional): Number of consecutive symbols to generate (default=1).
        - 'data' (ndarray, optional): Predefined OFDM data symbols.

    Returns:
    yVec : ndarray
        Generated OFDM signal in time domain, sampled at `Fsr` Hz.
    """
    
    # Default parameters
    if "gutter" not in s:
        s["gutter"] = False
    if "Nsym" not in s:
        s["Nsym"] = 1

    # Validate Midx (Modulation Order)
    if "data" not in s and (np.log2(s["Midx"]) % 2 != 0) and s["Midx"] != 2:
        raise ValueError("Midx must be 2 or an even power of 2.")
    
    # Generate or validate data symbols
    if "data" not in s:
        x = np.random.randint(0, s["Midx"], (s["N"], s["Nsym"]))
    else:
        if s["data"].shape[0] != s["N"] or s["data"].shape[1] > s["Nsym"]:
            raise ValueError("s['data'] must be an N x Nsym matrix.")
        x = s["data"]

    # OFDM Symbol Duration
    T = s["N"] / s["Fs"]  # Symbol duration (without cyclic prefix)
    Tg = s["Ng"] / s["Fs"]  # Guard interval duration
    Tsym = T + Tg  # Total OFDM symbol duration
    F = s["Fs"] / s["N"]  # Subcarrier spacing

    # Generate Modulated Symbols
    XVec = np.zeros_like(x, dtype=complex)
    
    if s["type"].upper() == "PSK":
        XVec = np.exp(1j * (np.pi / s["Midx"]) * (2 * x + 1))  # PSK Modulation
    elif s["type"].upper() == "QAM":
        XVec = (2 * (x % np.sqrt(s["Midx"])) - np.sqrt(s["Midx"]) + 1) + \
               1j * (2 * (x // np.sqrt(s["Midx"])) - np.sqrt(s["Midx"]) + 1)
        XVec /= np.sqrt(np.mean(np.abs(XVec) ** 2))  # Normalize for unit average power
    else:
        raise ValueError('Type must be "QAM" or "PSK".')

    # Handle symbol repetition if needed
    w = XVec.shape[1]
    if s["Nsym"] - w > 0:
        Nreps = s["Nsym"] // w
        Nremain = s["Nsym"] % w
        XVec = np.tile(XVec, (1, Nreps))
        XVec = np.hstack((XVec, XVec[:, :Nremain]))

    # Apply "Gutter" (Central Subcarrier Gap)
    if s["gutter"]:
        Ngut = 4  # Starlink-like 4-subcarrier gap
        XVec = np.fft.fftshift(XVec, axes=0)
        XVec[(s["N"] - Ngut) // 2 : (s["N"] - Ngut) // 2 + Ngut, :] = 0
        XVec = np.fft.fftshift(XVec, axes=0)

    # Convert to Time Domain (IFFT)
    Mx = np.sqrt(s["N"]) * np.fft.ifft(XVec, axis=0)

    # Add Cyclic Prefix
    MxCP = np.vstack((Mx[-s["Ng"]:, :], Mx))

    # Convert to Serial Stream
    xVec = MxCP.flatten()

    # Apply Doppler Shift and Frequency Offset
    if s["beta"] != 0 or s["Fc"] != s["Fcr"]:
        tVec = np.arange(len(xVec)) / s["Fs"]
        FD = -s["beta"] * s["Fc"]
        Fshift = FD + s["Fc"] - s["Fcr"]
        xVec *= np.exp(1j * 2 * np.pi * Fshift * tVec)

    # Add AWGN Noise
    yVec = xVec.copy()
    if not np.isnan(s["SNRdB"]):
        SNR_linear = 10 ** (s["SNRdB"] / 10)
        sigmaIQ = np.sqrt(1 / (2 * SNR_linear))
        noise = sigmaIQ * (np.random.randn(len(xVec)) + 1j * np.random.randn(len(xVec)))
        yVec += noise

    # Resample to Simulate Receiver Capture
    if s["Fsr"] != s["Fs"] and len(yVec) > 0:
        tVec = np.arange(len(yVec)) / s["Fs"]
        yVec = resample(yVec, int(len(yVec) * s["Fsr"] / s["Fs"]))

    return yVec
