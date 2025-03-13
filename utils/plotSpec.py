''' Ported from Matlab to Python by Jack R. Tschetter on 03/13/2025.
    The following key considerations were made.

    1. ) Handles Data Resampling
        Uses scipy.signal.resample() if the provided sampling rate (Fsr) differs from the desired spectrogram rate (Fs).
    2. ) Applies Frequency Shift if Needed
        If Fc ≠ Fcr, applies a complex exponential to shift the center frequency.
    3. ) Computes and Plots the Spectrogram
        Uses scipy.signal.spectrogram() with a Kaiser window.
        Normalizes data to match MATLAB’s behavior.
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import spectrogram, resample, windows

def plot_spec(data, tstart, tdur, Fcr, Fsr, Fc, Fs, NFFT, Stitle=None):
    """
    Plots a spectrogram of the provided signal.

    Parameters:
    data : ndarray
        Vector of complex samples.
    tstart : float
        Start time (relative to first sample).
    tdur : float
        Duration of data to process. If 0, plots entire signal.
    Fcr : float
        Center frequency of received signal.
    Fsr : float
        Sampling rate of provided signal.
    Fc : float
        Desired center frequency for the plot.
    Fs : float
        Desired sampling rate for spectrogram.
    NFFT : int
        Number of FFT points for the spectrogram.
    Stitle : str, optional
        Title for the spectrogram plot.

    Returns:
    fig : matplotlib.figure.Figure
        The generated spectrogram figure.
    """

    # Default title
    if Stitle is None:
        Stitle = f"Spectrogram centered at {Fc / 1e6} MHz"

    # If duration is 0, use the entire signal
    if tdur == 0:
        tdur = len(data) / Fsr

    # Extract the required segment
    seek_offset = int(np.floor(tstart * Fsr))
    Nblock = int(np.floor(tdur * Fsr))
    data = data[seek_offset : seek_offset + Nblock]

    # Time vector
    tyVec = np.arange(len(data)) / Fsr

    # Resample if necessary
    if Fs != Fsr:
        data = resample(data, int(len(data) * (Fs / Fsr)))
        tyVec = np.linspace(0, tdur, len(data))

    # Shift frequency if needed
    if Fc != Fcr:
        Fshift = Fcr - Fc
        data *= np.exp(1j * 2 * np.pi * Fshift * tyVec)

    # Compute Spectrogram
    fVec, tVec, Sxx = spectrogram(
        data / np.sqrt(np.mean(np.abs(data) ** 2)),  # Normalize power
        fs=Fs,
        window=windows.kaiser(NFFT, 0.5),
        nperseg=NFFT,
        noverlap=NFFT // 2,
        nfft=NFFT,
        scaling="density",
        mode="psd"
    )

    # Convert to dB scale
    Sxx_dB = 10 * np.log10(Sxx)
    Sxx_dB[Sxx_dB < -90] = -90  # Thresholding at -90 dB

    # Plot Spectrogram
    fig, ax = plt.subplots(figsize=(10, 5))
    c = ax.pcolormesh(tVec + tstart, fVec / 1e6, Sxx_dB, shading="auto", cmap="jet", vmin=-90)
    plt.colorbar(c, label="Power Density (dB/Hz)")

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Frequency (MHz)")
    ax.set_title(Stitle)

    plt.show()
    return fig
