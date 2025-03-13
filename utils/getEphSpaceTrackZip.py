''' Ported from Matlab to Python by Jack R. Tschetter on 03/13/2025.
    The following key considerations were made.

    1. ) Handles ZIP Extraction
        Uses Python’s built-in zipfile module to extract and read files.
    2. ) Validates File Names
        Ensures a max of 3 files, all with .zip extensions.
        Checks naming consistency.
    3. ) Parses Ephemeris Data
        Extracts date_created, ephemeris_start, ephemeris_stop, step_size.
        Reads epoch state vectors and covariance matrices.
    4. ) Manages Time Conversions
        Converts YYYYDOYhhmmss to standard datetime format.
'''

import zipfile
import os
import numpy as np
import datetime

def get_eph_space_track_zip(files, StarlinkID, verbosity=False):
    """
    Extracts and parses Starlink ephemeris data from ZIP files provided by space-track.org.

    Parameters:
    files : list
        List of filenames (max 3) corresponding to space-track.org downloads.
    StarlinkID : int
        ID of the Starlink satellite to extract data for.
    verbosity : bool, optional
        If True, prints detailed logs.

    Returns:
    data : dict
        Dictionary containing:
        - 'date_created': Datetime the ephemeris data was created (UTC)
        - 'ephemeris_start': Datetime of first ephemeris point (UTC)
        - 'ephemeris_stop': Datetime of last ephemeris point (UTC)
        - 'step_size': Time interval between ephemeris points (seconds)
        - 'epoch_datetime': List of timestamps for each ephemeris point (UTC)
        - 'epoch_state': NumPy array (6 x N) of [x, y, z, vx, vy, vz] states (ECI)
        - 'epoch_covariance': NumPy array (6 x 6 x N) covariance matrices
    """
    
    if len(files) > 3:
        raise ValueError("Will only accept up to 3 files.")

    # Validate file extensions
    for file in files:
        if not file.endswith(".zip"):
            raise ValueError(f"Invalid file format: {file}. Must be a .zip file.")

    # Ensure filenames only differ by the last digit (if multiple)
    filenames = [os.path.basename(f) for f in files]
    if len(set(fname[:-5] for fname in filenames)) != 1:
        raise ValueError("Filenames should differ only by the last digit.")

    # Pattern to search for Starlink ID
    pattern = f"_STARLINK-{StarlinkID}_"

    if verbosity:
        print(f"Looking for STARLINK-{StarlinkID}...")

    extracted_file = None

    # Search for the correct ephemeris file in the ZIP archives
    for file in files:
        with zipfile.ZipFile(file, 'r') as zip_ref:
            for file_name in zip_ref.namelist():
                if pattern in file_name:
                    extracted_file = file_name
                    if verbosity:
                        print(f"Found: {file_name}, extracting...")
                    zip_ref.extract(file_name, os.getcwd())  # Extract to current directory
                    break
        if extracted_file:
            break

    if not extracted_file:
        if verbosity:
            print(f"Warning: STARLINK-{StarlinkID} not found in provided files.")
        return {}

    # Parse the extracted ephemeris file
    data = {}
    epoch_datetime = []
    epoch_state = []
    epoch_covariance = []

    with open(extracted_file, 'r') as f:
        lines = f.readlines()

    header_flag = False
    iiEphemeris = 0

    for line in lines:
        if "created:" in line:
            data["date_created"] = datetime.datetime.strptime(line.replace("created:", "").strip(), "%Y-%m-%d %H:%M:%S.%f UTC")

        elif "ephemeris_start:" in line:
            parts = line.split("ephemeris_stop:")
            data["ephemeris_stop"] = datetime.datetime.strptime(parts[1].split("step_size:")[0].strip(), "%Y-%m-%d %H:%M:%S.%f UTC")
            parts = parts[0].split("ephemeris_start:")
            data["ephemeris_start"] = datetime.datetime.strptime(parts[1].strip(), "%Y-%m-%d %H:%M:%S.%f UTC")
            data["step_size"] = float(parts[0].split("step_size:")[1].strip())
            header_flag = True

        elif header_flag:
            values = line.strip().split()
            if len(values) == 7:
                # Convert ephemeris time
                yyyyDOYhhmmss = values[0]
                Y, DOY, H, M, S = int(yyyyDOYhhmmss[:4]), int(yyyyDOYhhmmss[4:7]), int(yyyyDOYhhmmss[7:9]), int(yyyyDOYhhmmss[9:11]), int(yyyyDOYhhmmss[11:13])
                fractional_sec = yyyyDOYhhmmss.split('.')[1] if '.' in yyyyDOYhhmmss else '0'
                Mo, D = datetime.datetime(Y, 1, 1) + datetime.timedelta(DOY - 1), datetime.datetime(Y, 1, 1).day + (DOY - 1)
                date = datetime.datetime(Y, Mo.month, D, H, M, S, int(fractional_sec))
                epoch_datetime.append(date)

                # Save ephemeris state vector
                state_vector = np.array(list(map(float, values[1:])))
                epoch_state.append(state_vector)

                # Parse covariance matrix (3 lines per entry)
                lower_diag = []
                for _ in range(3):
                    iiEphemeris += 1
                    cov_values = lines[iiEphemeris].strip().split()
                    if len(cov_values) != 6:
                        raise ValueError("Unexpected file length or corruption in covariance data.")
                    lower_diag.extend(list(map(float, cov_values)))

                # Convert lower diagonal into full symmetric matrix
                a = np.triu(np.ones((6, 6)))
                a[a > 0] = lower_diag
                epoch_covariance.append((a + a.T) / (np.eye(6) + 1))

    # Convert lists to NumPy arrays
    data["epoch_datetime"] = np.array(epoch_datetime)
    data["epoch_state"] = np.array(epoch_state).T
    data["epoch_covariance"] = np.array(epoch_covariance)

    if verbosity:
        print("Ephemeris file successfully processed.")

    return data
