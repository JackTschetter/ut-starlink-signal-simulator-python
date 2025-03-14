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

from scipy.signal import welch, windows
import numpy as np
import matplotlib.pyplot as plt

def plot_pwr(data, Nblocks, tstart, tdur, Fs, NFFT):
    """Plots the PSD and Power vs Time of the provided signal."""
    
    if Nblocks > 10 * tdur * Fs:
        raise ValueError("Provided number of blocks is too high for duration specified.")

    # Extract correct portion of data
    seek_offset = int(np.floor(tstart * Fs))
    Nblock = int(np.floor(tdur * Fs))
    data = data[seek_offset:seek_offset + Nblock]

    # Define frequency axis
    tVec = np.arange(len(data)) / Fs + tstart

    # **Fix 1: Ensure nperseg is valid**
    nperseg = min(NFFT, len(data) // 8)

    # Compute Welch PSD
    fVec, Syy = welch(
        data,
        window=windows.kaiser(nperseg, 3),
        nperseg=nperseg,  # **Fix 2: Explicitly set nperseg**
        nfft=NFFT,
        fs=Fs,
        return_onesided=False
    )

    # Plot Power Spectral Density
    plt.subplot(2, 1, 1)
    plt.plot(fVec / 1e6, 10 * np.log10(Syy))
    plt.grid(True)
    plt.xlabel("Frequency (MHz)")
    plt.ylabel("Power density (dB/Hz)")
    plt.title("Power Spectral Density Estimate")

    # **Fix 3: Ensure L_seg is valid**
    L_seg = max(1, int(len(data) / Nblocks))  # Prevent zero division

    # Compute Power vs Time
    Pwrdb = np.zeros(Nblocks)
    for i in range(Nblocks):
        segment = data[L_seg * i:L_seg * (i + 1)]
        if len(segment) == 0:
            continue
        
        # **Fix 4: Set nperseg for block-wise Welch computation**
        block_nperseg = min(round(L_seg / 2), NFFT)
        
        _, Syy_seg = welch(
            segment,
            window=windows.kaiser(block_nperseg, 3),
            nperseg=block_nperseg,
            nfft=NFFT,
            fs=Fs,
            return_onesided=False
        )
        Pwrdb[i] = max(10 * np.log10(np.sum(Syy_seg * Fs)), -180)  # dB/Hz * Hz = dB

    # Plot Power vs Time
    plt.subplot(2, 1, 2)
    tVec = np.linspace(0, tdur, Nblocks) + tstart
    plt.plot(tVec, Pwrdb)
    plt.grid(True)
    plt.xlabel("Time (s)")
    plt.ylabel("dB")
    plt.title("Power vs Time")
    plt.show()
