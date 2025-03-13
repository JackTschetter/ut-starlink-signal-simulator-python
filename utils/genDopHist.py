''' Ported from Matlab to Python by Jack R. Tschetter on 03/13/2025.
    The following key considerations were made.
    Struct Handling → Python dictionary (dict) for input.
    Physical Constants → scipy.constants.speed_of_light for accuracy.
    Interpolation → scipy.interpolate.interp1d for interp1.
    Vector Norms and Dot Products → numpy.linalg.norm and numpy.dot.
'''

import numpy as np
from scipy.interpolate import interp1d
from scipy.constants import speed_of_light
from datetime import timedelta

def gen_dop_hist(s):
    """
    Generates a Doppler history interpolated from the provided parameters.

    Parameters:
    s : dict
        Dictionary containing:
        - 'epoch_datetime' (ndarray): Nx1 array of datetime objects.
        - 'epoch_state' (ndarray): Nx6 array of ECEF state [x,y,z,vx,vy,vz] in km.
        - 'start_datetime' (datetime): Start of the desired Doppler history.
        - 'tdur' (float): Duration of the simulated Doppler history (seconds).
        - 'Fc' (float): Center frequency of the signal (Hz).
        - 'Fs' (float): Sampling frequency for Doppler history (Hz).
        - 'rx_state' (ndarray): 1x6 array of receiver ECEF state [x,y,z,vx,vy,vz] in km.

    Returns:
    fdop : ndarray
        Doppler frequency shift over time.
    epoch_datetimes : ndarray
        Corresponding timestamps for Doppler shift values.
    epoch_state_interp : ndarray
        Interpolated satellite ECEF states.
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

    # Compute relative position & velocity between satellite and receiver
    d = np.linalg.norm(s["rx_state"][:3] - epoch_state_interp[:3, :], axis=0)  # Distance
    rG = (s["rx_state"][:3] - epoch_state_interp[:3, :]) / d  # Unit vector

    v_rel = s["rx_state"][3:6] - epoch_state_interp[3:6, :]  # Relative velocity

    # Compute Doppler shift using df = f_app - f_c
    fdop = -1 / lamda * np.einsum('ij,ij->j', rG, v_rel)  # Dot product across columns

    return fdop, epoch_datetimes, epoch_state_interp
