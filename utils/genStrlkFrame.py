''' Ported from Matlab to Python by Jack R. Tschetter on 03/13/2025.
    The following key considerations were made.
    
    1. ) Handling Input Parameters → Uses a Python dictionary (dict) for flexibility.
    2. ) OFDM Frame Composition:
        Generates PSS (genPss()), SSS (genSss()), and data frames (genStrlkOFDM()).
        Ensures cyclic prefix and noise handling.
    3. ) Doppler Shift and Frequency Adjustments:
        Uses numpy.exp(1j * 2 * np.pi * fShift * tVec) to apply shifts.
    4. ) Noise Addition (AWGN):
        Uses numpy.random.randn() for Gaussian noise generation.
    5. )Resampling for Receiver Capture:
        Uses scipy.signal.resample() if Fsr < Fs.
'''

import numpy as np
from scipy.signal import resample

def gen_strlk_frame(s):
    """
    Generates a Starlink frame with PSS, SSS, and data, incorporating Doppler shifts and noise.

    Parameters:
    s : dict
        Dictionary containing:
        - 'SNRdB' (float): Signal-to-noise ratio in dB (use np.nan for no noise).
        - 'Fsr' (float): Receiver sample rate (Hz).
        - 'Fcr' (float): Receiver center frequency (Hz).
        - 'beta' (float): Doppler factor.
        - 'data' (ndarray, optional): Predefined OFDM data symbols.
        - 'type' (str, optional): Modulation type ('PSK' or 'QAM').
        - 'doppVec' (ndarray, optional): Phase history for Doppler shift.
        - 'sigma2_w' (float, optional): Noise variance.
        - 'channel' (object, optional): Communications channel object.

    Returns:
    frame : ndarray
        Starlink OFDM frame in time, sampled at `Fsr`, centered at `Fsr`, expressed at baseband.
    """
    
    # Calculate Noise Parameters
    SNR = 10 ** (s["SNRdB"] / 10)
    sigma2_w = 1 / SNR
    A = 1
    if "sigma2_w" in s:
        sigma2_w = s["sigma2_w"]
        A = np.sqrt(SNR * sigma2_w)

    # Frame Guard Interval Parameters
    Tfg = (68 / 15) * 1e-6  # Frame guard interval
    s["Fs"] = 240e6  # Starlink signal bandwidth
    Nfg = round(Tfg * s["Fs"])
    s["N"] = 1024  # Number of subcarriers
    Fdelta = 250e6  # Starlink channel spacing
    Fcii = get_closest_fch(s["Fcr"])  # Get closest Starlink channel center

    # Backup initial values
    SNRdB, beta, Fsr, Fcr = s["SNRdB"], s["beta"], s["Fsr"], s["Fcr"]
    
    # Modify settings for initial signal generation
    s.update({"SNRdB": np.nan, "beta": 0, "Fcr": Fcii, "Fsr": s["Fs"]})

    # Generate Starlink Frame Components
    PSS = gen_pss()  # Primary Synchronization Sequence
    SSS = gen_sss()  # Secondary Synchronization Sequence

    if "data" not in s:
        # Generate predefined data blocks based on observed Starlink structure
        s.update({"Nsym": 4, "type": "PSK", "Midx": 4})
        Data = gen_strlk_ofdm(s)  # First block (4PSK)

        s.update({"Nsym": 4, "type": "QAM", "Midx": 16})
        Data = np.vstack((Data, gen_strlk_ofdm(s)))  # Second block (16QAM)

        s.update({"Nsym": 292, "type": "QAM", "Midx": 4})
        Data = np.vstack((Data, gen_strlk_ofdm(s)))  # Remaining blocks (4QAM)
    else:
        s["Nsym"] = 300
        Data = gen_strlk_ofdm(s)

    # Frame Guard
    Fg = np.zeros(Nfg)

    # Restore previous settings
    s.update({"Fcr": Fcr, "Fsr": Fsr, "beta": beta, "SNRdB": SNRdB})

    # Concatenate full frame
    frame = np.concatenate((PSS, SSS, Data, Fg))

    # Apply Doppler Shift and Frequency Bias
    if "beta" in s:
        if "channel" in s:
            FD = -s["beta"] * get_closest_fch(s["Fcr"])  # Doppler shift
            s["channel"].DirectPathDopplerShift = FD
            frame = s["channel"](frame)
        else:
            FD = -s["beta"] * get_closest_fch(s["Fcr"])
            fShift = FD + get_closest_fch(s["Fcr"]) - s["Fcr"]
            tVec = np.arange(len(frame)) / s["Fs"]
            frame *= np.exp(1j * 2 * np.pi * fShift * tVec)
    
    elif "doppVec" in s:
        if "channel" in s:
            fShift = get_closest_fch(s["Fcr"]) - s["Fcr"]
            offsetVec = 2 * np.pi * (s["doppVec"] - s["doppVec"][0]) / s["Fs"]
            s["channel"].DirectPathDopplerShift = s["doppVec"][0]
            frame = s["channel"](frame)
            Phihist = np.cumsum(offsetVec)
            frame *= np.exp(1j * (Phihist + 2 * np.pi * fShift * tVec))
        else:
            fShift = get_closest_fch(s["Fcr"]) - s["Fcr"]
            offsetVec = 2 * np.pi * s["doppVec"] / s["Fs"]
            Phihist = np.cumsum(offsetVec)
            frame *= np.exp(1j * (Phihist + 2 * np.pi * fShift * tVec))

    # Scale the frame signal
    frame *= A

    # Add AWGN Noise to PSS, SSS, and Frame Guard
    if not np.isnan(s["SNRdB"]):
        sigmaIQ = np.sqrt(sigma2_w / 2)
        noise = sigmaIQ * (np.random.randn(len(frame)) + 1j * np.random.randn(len(frame)))
        frame += noise

    # Resample if necessary
    if s["Fsr"] < s["Fs"]:
        tVec = np.arange(len(frame)) / s["Fs"]
        frame = resample(frame, int(len(frame) * s["Fsr"] / s["Fs"]))

    return frame

# Helper Functions
def get_closest_fch(Fcr):
    """ Returns the closest Starlink channel frequency. """
    return round(Fcr / 1e6) * 1e6  # Example approximation

def gen_pss():
    """ Placeholder for PSS generation. """
    return np.random.randn(1056) + 1j * np.random.randn(1056)

def gen_sss():
    """ Placeholder for SSS generation. """
    return np.random.randn(1056) + 1j * np.random.randn(1056)

def gen_strlk_ofdm(s):
    """ Placeholder for Starlink OFDM generation. """
    return np.random.randn(1056 * s["Nsym"]) + 1j * np.random.randn(1056 * s["Nsym"])
