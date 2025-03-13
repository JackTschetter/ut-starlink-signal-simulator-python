import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import resample
from scipy.fftpack import fft

import os
import sys
sys.path.append(os.path.abspath("../utils"))  # Adjust the path as needed

from genStrlkOFDM import gen_strlk_ofdm
from getClosestFch import get_closest_fch
from plotSpec import plot_spec
from plotPwr import plot_pwr

# ------------------ Options ------------------ #
specPlot_en = True
pwrPlot_en = True
constPlot_en = True

# ------------------ Parameters ------------------ #
NFFT = 2 ** 10  # FFT points for spectrogram

# ------------------ OFDM Input Parameters ------------------ #
s = {
    "Nsym": 100,       # Number of OFDM symbols
    "SNRdB": np.nan,   # Signal to noise ratio (not used in this script)
    "Fsr": 200e6,      # Receiver sample rate
    "Fcr": 12.075e9,   # Receiver center frequency
    "beta": 0,         # Carrier frequency offset
    "Midx": 4,         # Modulation index
    "type": "QAM",     # Modulation type
    "Fc": 11.805e9     # Channel center frequency
}

# ------------------ Generate Signal ------------------ #
y = gen_strlk_ofdm(s)  # Generate OFDM signal

# ------------------ Spectrogram ------------------ #
if specPlot_en:
    Stitle = f"Spectrogram of {s['Nsym']} Starlink symbols, centered at closest channel center"
    F = 240e6 / 1024
    chIdx = round((s["Fcr"] / 1e9 - 10.7 - F / 2 / 1e9) / 0.25 + 0.5)
    Fcii = 10.7e9 + F / 2 + 250e6 * (chIdx - 0.5)
    tstart = 0
    tdur = len(y) / 240e6
    plot_spec(y, tstart, tdur, s["Fcr"], s["Fsr"], Fcii, 240e6, NFFT, Stitle)

# ------------------ Power Plot ------------------ #
if pwrPlot_en:
    plot_pwr(y, 20, 0, len(y) / s["Fsr"], s["Fsr"], NFFT)

# ------------------ Constellation Plot ------------------ #
sc_mag_threshold = 0.95  # Threshold for subcarrier selection

if constPlot_en:
    # Resample and center the signal to the channel
    if s["Fsr"] != 240e6:
        tVec = np.arange(len(y)) / s["Fsr"]
        y = resample(y, int(len(y) * 240e6 / s["Fsr"]))

    # Frequency Shift Correction
    f_shift = get_closest_fch(s["Fcr"]) - s["Fcr"]
    if f_shift != 0:
        tVec = np.arange(len(y)) / 240e6
        y *= np.exp(-1j * 2 * np.pi * f_shift * tVec)

    # Zero padding
    yVec = np.concatenate([y, np.zeros(1056 * np.ceil(len(y) / 1056) - len(y))])

    # Remove Cyclic Prefix
    y = yVec.reshape((-1, 1056))
    y = y[:, 32:]

    # Take FFT
    Y = (1 / np.sqrt(1024)) * fft(y, axis=1)

    # Find good subcarriers
    idx = np.where(np.abs(Y) / np.max(np.abs(Y)) > sc_mag_threshold)

    # Set gutter values to NaN
    Y[:2, :] = np.nan
    Y[-2:, :] = np.nan

    # Plot Constellation Diagram
    print("Press Enter to start constellation plotting loop...")
    input()  # Wait for user input before starting the loop
    plt.figure()
    for ii in range(Y.shape[1]):
        plt.scatter(Y[idx][:, ii].real, Y[idx][:, ii].imag, marker="*")
        plt.title(f"Constellation Plot - Symbol {ii+1}")
        plt.xlabel("In-phase")
        plt.ylabel("Quadrature")
        plt.grid(True)
        plt.pause(1)  # Pause for 1 second before plotting next symbol
        plt.clf()
