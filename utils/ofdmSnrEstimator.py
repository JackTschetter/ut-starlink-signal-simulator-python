''' Ported from Matlab to Python by Jack R. Tschetter on 03/13/2025.
    The following key considerations were made.
    1. ) Handles OFDM Reshaping
        Uses numpy.reshape() to segment y into OFDM symbols.
    2. ) Iterative Maximum Likelihood Estimation
        Uses NumPy operations to compute noise variance (sigma_2) and SNR (snr_dB).
    3. ) Threshold-based Iteration Control
        Iterates until the threshold condition is met.
'''

import numpy as np
from scipy.stats import norm

def ofdm_snr_estimator(s):
    """
    Estimates the SNR from a stream of OFDM symbols using the 
    two-step non-data-aided estimation method.

    Parameters:
    s : dict
        Dictionary containing:
        - 'N' (int): Number of subcarriers.
        - 'Ng' (int): Length of cyclic prefix (CP).
        - 'Nsym' (int): Number of consecutive symbols in `y`.
        - 'y' (ndarray): (Nsym * (N + Ng), 1) vector of OFDM samples in time.
        - 'threshold' (float, optional): Threshold for ML iteration (default: 1e-6).

    Returns:
    snr_dB : float
        Estimated SNR in dB.
    sigma_2 : float
        Estimated noise variance.
    nitter : int
        Number of iterations in ML estimation.
    """

    # Set default threshold if not provided
    threshold = s.get("threshold", 1e-6)

    # Reshape `y` into OFDM symbols
    y = s["y"].reshape((s["N"] + s["Ng"], -1))

    u = np.arange(1, s["Ng"] + 1)  # CP indices
    Ju = np.zeros(len(u))
    ksi = np.zeros(len(u))
    Lhat = np.zeros(len(u))

    for uu in range(len(u) - 1, -1, -1):
        Ju[uu] = (1 / (2 * s["Nsym"] * (s["Ng"] - uu + 1))) * np.sum(
            np.abs(y[uu:s["Ng"], :] - y[s["N"] + uu:s["N"] + s["Ng"], :]) ** 2
        )

        if uu == len(u) - 1:
            ksi[uu] = Ju[uu]
        else:
            ksi[uu] = Ju[uu] - (1 - 1 / (s["Ng"] - uu + 1)) * Ju[uu]

        # Given L = uu
        sigma_2 = Ju[uu]
        Lhat[uu] = np.prod(
            norm.pdf(ksi[uu:], sigma_2 / (s["Ng"] - uu + 1), sigma_2**2 / (s["Nsym"] * (s["Ng"] - uu + 1) ** 2))
        ) ** (1 / (s["Ng"] - uu + 1))

    # Find minimum Lhat
    idx = np.argmin(Lhat)
    sigma_2_hat = Ju[idx]

    # Compute signal power estimate
    s_hat = (1 / (s["Nsym"] * s["N"])) * np.sum(
        np.real(y[0 : (s["Ng"] + Lhat[idx]), :] * np.conj(y[s["N"] : (s["Ng"] + Lhat[idx] + s["N"]), :]))
    )

    # Compute energy per symbol
    Ey_hat = (1 / (s["Nsym"] * (s["N"] - s["Ng"]))) * np.sum(np.abs(y) ** 2)

    # Initial estimates
    theta_c_hat = np.array([s_hat, sigma_2_hat, Ey_hat])

    # Iterative Estimation Process
    nitter = 0
    while True:
        varShat = ((s["Ng"] + Lhat[idx]) * (Ey_hat**2 + s_hat**2)) / (2 * s["Nsym"] * s["Ng"]**2)
        varsigma2hat = sigma_2_hat**2 / (s["Nsym"] * (s["Ng"] - Lhat[idx]))
        varEyhat = Ey_hat**2 / (s["Nsym"] * (s["N"] - Lhat[idx]))
        covShatsigma2hat = -sigma_2_hat**2 / (2 * s["Nsym"] * s["Ng"])

        # Construct covariance matrix
        C = np.array([
            [varShat, covShatsigma2hat, 0],
            [covShatsigma2hat, varsigma2hat, 0],
            [0, 0, varEyhat]
        ])

        # Define transformation matrix
        H = np.array([
            [1, 0],
            [0, 1],
            [1, 1]
        ])

        # Compute new estimates
        thetar = np.linalg.inv(H.T @ np.linalg.inv(C) @ H) @ H.T @ np.linalg.inv(C) @ theta_c_hat
        e = np.abs(thetar - theta_c_hat[:2]) ** 2
        theta_c_hat[:2] = thetar
        nitter += 1

        if np.sum(e < threshold) == 2:
            break

    # Compute final SNR estimate
    sigma_2 = thetar[1]
    SNR_hat = thetar[0] / thetar[1]
    snr_dB = 10 * np.log10(SNR_hat)

    return snr_dB, sigma_2, nitter
