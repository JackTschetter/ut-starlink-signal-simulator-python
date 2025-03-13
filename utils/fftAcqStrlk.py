''' Ported from Matlab to Python by Jack R. Tschetter on 03/13/2025.
    The following key considerations were made.
    Struct Handling → Used Python dictionaries or a dataclass for input/output.
    Matrix & Vector Ops → Used NumPy for efficient calculations.
    Resampling → Used scipy.signal.resample for equivalent MATLAB resample.
    Fourier Transforms → Used numpy.fft.fft and numpy.fft.ifft.
    Looping and Indexing → Maintained efficiency by vectorizing where possible.
'''

import numpy as np
import scipy.signal as signal

def fft_acq_strlk(s):
    """
    Performs FFT acquisition based on PSS and SSS for Starlink signals.

    Parameters:
    s : dict
        Dictionary containing:
        - Fsr (float): Receiver sampling frequency in Hz
        - Fcr (float): Receiver center frequency in Hz
        - fmax (float): Max frequency for search in Hz
        - fstep (float): Frequency step in Hz
        - fstepFine (float): Fine frequency step in Hz
        - fmin (float): Min frequency for search in Hz
        - y (ndarray): Ns x 1 vector of baseband data
        - Pfa (float): False alarm probability

    Returns:
    data : dict
        Dictionary containing:
        - grid (ndarray): Acquisition grid
        - tauvec (ndarray): Vector of considered delay times
        - fdvec (ndarray): Vector of considered Doppler frequencies
        - fdfine (float): Fine Doppler frequency estimate
        - tau (float): Time offset estimate
        - Nc (int): Length of local replica
        - sigma2_c (float): Variance of local replica
        - pkVec (ndarray): Local replica used
    """

    # Constants
    N = 1024  # Number of subcarriers
    Ng = 32   # CP length
    Ns = N + Ng  # Samples per symbol
    Ndsym = 300  # Data-carrying symbols
    Fs = 240e6  # Bandwidth
    Nk = int((1 / 750) * Fs)  # Samples per frame

    # Generate known sequences (PSS, SSS)
    PSS = gen_pss()  # Implement gen_pss()
    SSS = gen_sss()  # Implement gen_sss()
    c = np.concatenate((PSS, SSS))  # Known frame
    pkVec = c
    c = np.concatenate((c, np.zeros(Nk - len(c))))  # Zero-padding

    # Resample to full bandwidth if needed
    y = s["y"]
    if s["Fsr"] != Fs:
        tVec = np.arange(len(y)) / s["Fsr"]
        y = signal.resample(y, int(len(y) * Fs / s["Fsr"]))

    buffer = 100
    if len(y) >= Nk - buffer:
        y = np.concatenate((y, np.zeros(Nk - len(y))))
    else:
        raise ValueError("Error! Data should be long enough to contain at least 1 frame.")

    # Remove receiver bias
    Fshift = get_closest_fch(s["Fcr"]) - s["Fcr"]
    if Fshift != 0:
        tVec = np.arange(len(y)) / Fs
        y *= np.exp(-1j * 2 * np.pi * Fshift * tVec)

    # Derived parameters
    fdvec = np.arange(s["fmin"], s["fmax"] + s["fstep"], s["fstep"])
    tauvec = np.arange(Nk) / Fs
    known = c[c != 0]
    Nc = len(known)
    sigma2_c = np.var(known)
    Ngrid = len(fdvec) * len(tauvec)

    C = np.fft.fft(c)
    Nac = len(y) // Nk
    y = y[:Nac * Nk]  # Truncate to a multiple of frame length

    # Initialize acquisition grid
    grids = np.zeros((len(fdvec), Nk))
    Sk = np.zeros((len(fdvec), len(c), Nac), dtype=complex)

    # Generate acquisition grid
    for jj in range(Nac):
        for ii, fd in enumerate(fdvec):
            pre_idcs = slice(jj * Nk, (jj + 1) * Nk)
            tau_pre = tauvec[:Nk]
            x = y[pre_idcs]

            beta = -fd / get_closest_fch(s["Fcr"])
            FsD = (1 + beta) * Fs
            x, tau_post = signal.resample(x, len(x), t=tau_pre)

            th_hat = 2 * np.pi * fd * tau_post
            x_tilde = x * np.exp(-1j * th_hat)
            x_tilde = np.concatenate((x_tilde, np.zeros(Nk - len(x_tilde))))[:Nk]

            X_tilde = np.fft.fft(x_tilde)
            Zr = X_tilde * np.conj(C)
            sk = np.fft.ifft(Zr)
            Sk[ii, :, jj] = sk[:len(c)]

    grids[:, :] = np.sum(np.abs(Sk), axis=2)

    # Process acquisition
    tau_prob = np.max(grids, axis=0)
    mx_val = np.max(tau_prob)
    tau_idx = np.argmax(tau_prob)
    tau = tauvec[tau_idx]
    fdcoarse = fdvec[np.argmax(grids[:, tau_idx])]

    # Fine Doppler search
    fdfine = fdcoarse
    if s["fstep"] != 1:
        fdvec_fine = np.arange(fdcoarse - s["fstep"] // 2, fdcoarse + s["fstep"] // 2 + s["fstepFine"], s["fstepFine"])
        thmat = 2 * np.pi * (tauvec - tauvec[Nk - 1])[:, None] * fdvec_fine
        z_fine = np.zeros(len(fdvec_fine))

        ctau = np.roll(c, tau_idx - 1)
        for jj in range(Nac):
            iidum = slice(jj * Nk, (jj + 1) * Nk)
            x1_shift = y[iidum, None] * np.exp(-1j * thmat[:len(iidum), :])
            z_fine += np.abs(x1_shift.T @ ctau)

        fdfine = fdvec_fine[np.argmax(z_fine)]

    # Output results
    data = {
        "grid": grids,
        "tauvec": tauvec,
        "fdvec": fdvec,
        "fdfine": fdfine,
        "tau": tau,
        "Nc": Nc,
        "sigma2_c": sigma2_c,
        "pkVec": pkVec,
        "Sk": Sk
    }

    return data

# Placeholder functions for missing dependencies
def gen_pss():
    """ Placeholder for PSS generation. """
    return np.random.randn(512) + 1j * np.random.randn(512)

def gen_sss():
    """ Placeholder for SSS generation. """
    return np.random.randn(512) + 1j * np.random.randn(512)

def get_closest_fch(fcr):
    """ Placeholder for function that determines the closest frequency channel. """
    return round(fcr / 1e6) * 1e6  # Example approximation
