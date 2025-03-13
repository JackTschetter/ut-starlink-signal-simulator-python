''' Ported from Matlab to Python by Jack R. Tschetter on 03/13/2025.
    The following key considerations were made.

    1. ) Uses SciPy's Welch Method for PSD
        Uses scipy.signal.welch() to estimate power spectral density (PSD).
        Applies Kaiser windowing to match MATLAB behavior.
    2. ) Handles FFT Processing
        Uses numpy.fft.fft() for computing the FFT-based power estimation.
    3. ) Power vs. Time Calculation
        Segments data into blocks, computes power per block, and plots power vs. time.
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch, windows

def plot_pwr(data, Nblocks, tstart, tdur, Fs, NFFT):
    """
    Plots the Power Spectral Density (PSD) and Power vs. Time of the provided signal.

    Parameters:
    data : ndarray
        Vector of complex samples.
    Nblocks : int
        Number of blocks for power vs time plot.
    tstart : float
        Start time (relative to first sample).
    tdur : float
        Duration of data to process.
    Fs : float
        Sampling rate of the signal.
    NFFT : int
        Number of FFT points for the PSD plot.

    Returns:
    fig : matplotlib.figure.Figure
        The generated figure.
    """

    if Nblocks > 10 * tdur * Fs:
        raise ValueError("Provided number of blocks is too high for the specified duration.")

    # Extract the required segment
    seek_offset = int(np.floor(tstart * Fs))
    Nblock = int(np.floor(tdur * Fs))
    data = data[seek_offset : seek_offset + Nblock]

    # Time Vector
    tVec = np.arange(len(data)) / Fs + tstart

    # Compute PSD using Welch's Method
    fVec, Syy = welch(
        data, 
        window=windows.kaiser(NFFT, 3), 
        nperseg=NFFT, 
        noverlap=NFFT//2, 
        nfft=NFFT, 
        fs=Fs, 
        return_onesided=False, 
        scaling="density"
    )

    # Convert to dB scale
    Syy_dB = 10 * np.log10(Syy)

    # Plot PSD
    fig, axs = plt.subplots(2, 1, figsize=(10, 6))

    axs[0].plot(fVec / 1e6, Syy_dB)
    axs[0].grid(True)
    axs[0].set_xlabel("Frequency (MHz)")
    axs[0].set_ylabel("Power Density (dB/Hz)")
    axs[0].set_title("Power Spectral Density Estimate")

    # Power vs. Time Calculation
    L_seg = len(data) // Nblocks
    Pwrdb = np.zeros(Nblocks)

    for i in range(Nblocks):
        segment = data[L_seg * i : L_seg * (i + 1)]
        if len(segment) == 0:
            continue
        _, Syy_seg = welch(segment, window=windows.kaiser(round(L_seg / 2), 3), nfft=NFFT, fs=Fs, return_onesided=False)
        Pwrdb[i] = max(10 * np.log10(np.sum(Syy_seg * Fs)), -180)  # dB/Hz * Hz = dB

    # Time axis for power vs time plot
    tVec_blocks = np.linspace(0, tdur, Nblocks) + tstart

    axs[1].plot(tVec_blocks, Pwrdb)
    axs[1].grid(True)
    axs[1].set_xlabel("Time (s)")
    axs[1].set_ylabel("dB")
    axs[1].set_title("Power vs Time")

    plt.tight_layout()
    plt.show()

    return fig
