''' Ported from Matlab to Python by Jack R. Tschetter on 03/13/2025.
    The following key considerations were made.
    Uses the same formula as MATLAB to compute the closest Starlink channel.
    Ensures floating-point precision when calculating ch_idx.
    Returns the nearest valid Starlink channel center frequency.
'''

def get_closest_fch(Fc):
    """
    Returns the center frequency of the closest Starlink OFDM channel.

    Parameters:
    Fc : float
        Center frequency of receiver in Hz.

    Returns:
    Fcii : float
        Center frequency of the closest Starlink OFDM channel in Hz.
    """
    F = 240e6 / 1024  # Subcarrier spacing
    ch_idx = round(((Fc / 1e9) - 10.7 - (F / 2 / 1e9)) / 0.25 + 0.5)
    Fcii = 10.7e9 + (F / 2) + (250e6 * (ch_idx - 0.5))
    
    return Fcii
