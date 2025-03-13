''' Ported from Matlab to Python by Jack R. Tschetter on 03/13/2025.
    The following key considerations were made.

    1. ) Noise Handling
        Computes noise variance (sigma2_w) based on SNR.
        Adds AWGN noise using numpy.random.randn().
    2. ) Frame Generation
        Calls genStrlkFrame() for each frame.
        Uses numpy.random.binomial() to determine frame presence when prob is set.
    3. ) Doppler and Frequency Shift
        Uses numpy.exp(1j * 2 * np.pi * fShift * tVec) to apply Doppler shifts.
    4. ) Resampling
        Uses scipy.signal.resample() for handling sample rate conversion.       
'''

import numpy as np
from scipy.signal import resample

def gen_strlk_stream(s):
    """
    Generates a stream of Starlink frames with OFDM symbols.

    Parameters:
    s : dict
        Dictionary containing:
        - 'Midx' (int): Modulation index (2 = BPSK, 4 = 4QAM, etc.).
        - 'type' (str): Modulation type ('PSK' or 'QAM').
        - 'SNRdB' (float): Signal-to-noise ratio in dB.
        - 'Fsr' (float): Receiver sample rate (Hz).
        - 'Fcr' (float): Receiver center frequency (Hz).
        - 'beta' (float): Doppler factor.
        - 'Tdur' (float): Duration of the stream.
        - 'data' (ndarray, optional): Predefined OFDM data.
        - 'prob' (float, optional): Probability of frame presence.
        - 'present' (ndarray, optional): Vector of 1s and 0s indicating frame presence.
        - 'tau' (float, optional): Time offset for first frame slot.
        - 'doppVec' (ndarray, optional): Phase history for Doppler shift.
        - 'sigma2_w' (float, optional): Noise variance.

    Returns:
    y : ndarray
        Starlink OFDM symbol stream in time, sampled at `Fsr`, centered at `Fsr`,
        and expressed at baseband.
    """

    # Compute Noise Parameters
    SNR = 10 ** (s["SNRdB"] / 10)
    sigma2_w = 1 / SNR
    A = 1
    if "sigma2_w" in s:
        sigma2_w = s["sigma2_w"]
        A = np.sqrt(SNR * sigma2_w)

    # Ensure proper `prob` and `tau` defaults
    prob = s.get("prob", 1)
    tau = s.get("tau", 0)

    # Constants
    Fs = 240e6  # Starlink signal bandwidth
    Tframe = 1 / 750  # Starlink Frame duration
    Nfr = int(np.ceil(s["Tdur"] / Tframe))  # Number of whole Starlink Frames in stream
    Ns = int(np.floor(s["Tdur"] * Fs))  # Total samples in stream
    Nframe = int(np.floor(Tframe * Fs))  # Samples per Starlink frame

    # Handle Doppler Phase History Conflicts
    if "beta" in s and "doppVec" in s:
        raise ValueError("You can specify either 'beta' or 'doppVec', not both.")

    if "doppVec" in s and len(s["doppVec"]) != Ns:
        raise ValueError(f"'doppVec' must be {Ns} samples long.")

    # Determine Frame Presence
    if "present" in s:
        if "prob" in s:
            raise ValueError("Specify either 'prob' or 'present', not both.")
        present = np.array(s["present"])
        if len(present) < Nfr:
            present = np.pad(present, (0, Nfr - len(present)), mode='constant', constant_values=0)
        else:
            present = present[:Nfr]
    else:
        present = np.random.binomial(1, prob, Nfr)

    # Generate Starlink Stream
    y = np.zeros(Ns, dtype=complex)

    for ii in range(Nfr):
        if present[ii]:  # Place frame if present
            fr = {
                "SNRdB": np.nan,
                "Fsr": Fs,
                "Fcr": get_closest_fch(s["Fsr"]),
                "beta": 0
            }
            if "type" in s:
                fr["type"] = s["type"]
            if "data" in s:
                fr["data"] = s["data"]

            frame = gen_strlk_frame(fr)

            start_idx = ii * Nframe
            end_idx = start_idx + len(frame)
            y[start_idx:end_idx] = frame

    # Apply Initial Time Offset
    Ntau = int(np.round(tau * Fs))
    y = np.concatenate((np.zeros(Ntau, dtype=complex), y))

    # Apply Doppler and Frequency Shift
    if "beta" in s:
        if not (s["beta"] == 0 and get_closest_fch(s["Fcr"]) == s["Fcr"]):
            FD = -s["beta"] * get_closest_fch(s["Fcr"])  # Doppler shift
            fShift = FD + get_closest_fch(s["Fcr"]) - s["Fcr"]  # Total frequency shift

            tVec = np.arange(len(y)) / Fs
            y *= np.exp(1j * 2 * np.pi * fShift * tVec)

    elif "doppVec" in s:
        fShift = get_closest_fch(s["Fcr"]) - s["Fcr"]  # Total frequency shift
        tVec = np.arange(len(y)) / Fs
        offsetVec = 2 * np.pi * s["doppVec"] / Fs
        Phihist = np.cumsum(offsetVec)
        y *= np.exp(1j * (Phihist + 2 * np.pi * fShift * tVec))

    # Add AWGN Noise
    if not np.isnan(s["SNRdB"]):
        sigmaIQ = np.sqrt(sigma2_w / 2)
        noise = sigmaIQ * (np.random.randn(len(y)) + 1j * np.random.randn(len(y)))
        y += noise

    # Resample if necessary
    if s["Fsr"] < Fs:
        y = resample(y, int(len(y) * s["Fsr"] / Fs))

    return y

# Helper Functions
def get_closest_fch(Fc):
    """ Returns the closest Starlink channel frequency. """
    F = 240e6 / 1024  # Subcarrier spacing
    ch_idx = round(((Fc / 1e9) - 10.7 - (F / 2 / 1e9)) / 0.25 + 0.5)
    Fcii = 10.7e9 + (F / 2) + (250e6 * (ch_idx - 0.5))
    return Fcii

def gen_strlk_frame(s):
    """ Placeholder for Starlink frame generation. """
    return np.random.randn(1056) + 1j * np.random.randn(1056)
