''' Ported from Matlab to Python by Jack R. Tschetter on 03/13/2025.
    The following key considerations were made.

    1. ) Uses numpy.fromfile() for Binary Reading
        Reads 16-bit integer samples (int16).
        Handles complex numbers by interleaving real and imaginary parts.
    2. ) Implements fseek Equivalent
        Uses file.seek() to jump to the correct start time (tSeek).
    3. ) Checks for Data Availability
        Ensures there is enough data to read the requested duration.
'''

import numpy as np

def read_from_bin(filepath, Fsr, tDur, tSeek):
    """
    Reads complex IQ samples from a binary file.

    Parameters:
    filepath : str
        Path to the binary file.
    Fsr : float
        Sampling rate of the data (Hz).
    tDur : float
        Duration to read from the file (seconds).
    tSeek : float
        Starting time in the file (seconds).

    Returns:
    y : ndarray
        (tDur * Fsr) x 1 vector of complex samples.
    """

    # Compute number of samples to read
    Ns = int(np.floor(tDur * Fsr))

    # Compute byte offset (each sample is 4 bytes: 2 bytes for I and 2 bytes for Q)
    seek_offset = int(np.floor(tSeek * Fsr)) * 4

    # Open file and seek to position
    with open(filepath, "rb") as f:
        f.seek(seek_offset, 0)  # Seek from start of file
        data = np.fromfile(f, dtype=np.int16, count=2 * Ns)  # Read interleaved I/Q

    # Ensure we read the full requested data
    if len(data) < 2 * Ns:
        raise ValueError("Insufficient data in file for requested duration.")

    # Convert interleaved I/Q data to complex format
    y = data[0::2] + 1j * data[1::2]

    return y
