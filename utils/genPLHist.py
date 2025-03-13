''' Ported from Matlab to Python by Jack R. Tschetter on 03/13/2025.
    The following key considerations were made.
    Struct Handling → Use Python dictionary (dict) for input.
    Physical Constants → Use scipy.constants.speed_of_light for accuracy.
    Interpolation → Use scipy.interpolate.interp1d for interp1.
    Vector Norms → Use numpy.linalg.norm for distance calculation.
    Path Loss Calculation → Implement path loss model using the exponent parameter.
'''

import numpy as np
from scipy.interpolate import interp1d
from scipy.constants import speed_of_light
from datetime import timedelta

def gen_pl_hist(s):
    """
    Generates a path loss history interpolated from the provided parameters.

    Parameters:
    s : dict
        Dictionary containing:
        - 'epoch_datetime' (ndarray): Nx1 array of datetime objects.
        - 'epoch_state' (ndarray): Nx6 array of ECEF state [x,y,z,vx,vy,vz] in km.
        - 'start_datetime' (datetime): Start of the desired path loss history.
        - 'tdur' (float): Duration of the simulated path loss history (seconds).
        - 'Fc' (float): Center frequency of the signal (Hz).
        - 'Fs' (float): Sampling frequency for path loss history (Hz).
        - 'rx_state' (ndarray): 1x6 array of receiver ECEF state [x,y,z,vx,vy,vz] in km.
        - 'exp' (float): Path loss exponent.

    Returns:
    epoch_datetimes : ndarray
        Corresponding timestamps for the computed path loss.
    epoch_state_interp : ndarray
        Interpolated satellite ECEF states.
    path_loss : ndarray
        Path loss over time in dB.
    """

    c_km_s = speed_of_light / 1e3  # Convert to km/s
    lamda = c_km_s / s["Fc"]  # Wavelength in km
    Nsamples = int(s["Fs"] * s["tdur"])  # Total samples

    # Generate time vector
    epoch_datetimes = np.array([s["start_datetime"] + timedelta(seconds=i / s["Fs"]) for i in range(Nsamples)])

    # Interpolate the satellite ECEF state
    epoch_state_interp = np.zeros((6, Nsamples))
    for i in range(6):
        interp_func = interp1d(s["epoch_datetime"], s["epoch_state"][:, i], kind="linear", fill_value="extrapolate")
        epoch_state_interp[i, :] = interp_func(epoch_datetimes)

    # Compute relative position & distance between satellite and receiver
    d = np.linalg.norm(s["rx_state"][:3] - epoch_state_interp[:3, :], axis=0)  # Distance in km

    # Compute Path Loss using Free-Space Path Loss Model
    path_loss = 20 * np.log10(4 * np.pi * d / lamda) + 10 * s["exp"] * np.log10(d)

    return epoch_datetimes, epoch_state_interp, path_loss
