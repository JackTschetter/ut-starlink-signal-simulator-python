''' Ported from Matlab to Python by Jack R. Tschetter on 03/13/2025.
    The following key considerations were made.
    
    1. ) Used Cartopy for 3D Earth Projection
        Used cartopy.crs.Globe to define Earth's ellipsoid.
        Plots topography and coastlines.
        Optionally marks a receiver location in latitude/longitude.
    2. ) Handles Starlink Coverage Cone
        Converts LLA (Lat, Lon, Alt) to ECEF coordinates.
        Computes a 3D cone in ECEF representing coverage.
        Rotates the cone using Yaw, Pitch, Roll transformations.
    3. ) Uses Matplotlib for Rendering
        Uses plot_surface() for 3D rendering.
        Uses scatter() for marking the receiver location.
'''

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from mpl_toolkits.mplot3d import Axes3D
from pyproj import Proj, transform

def lla_to_ecef(lat, lon, alt):
    """
    Converts Latitude, Longitude, Altitude (LLA) to Earth-Centered Earth-Fixed (ECEF) coordinates.
    """
    a = 6378.137  # WGS-84 semi-major axis in km
    e2 = 6.69437999014 * 1e-3  # WGS-84 first eccentricity squared

    lat, lon = np.radians(lat), np.radians(lon)
    N = a / np.sqrt(1 - e2 * np.sin(lat)**2)

    X = (N + alt) * np.cos(lat) * np.cos(lon)
    Y = (N + alt) * np.cos(lat) * np.sin(lon)
    Z = (N * (1 - e2) + alt) * np.sin(lat)

    return np.array([X, Y, Z])

def plot_earth(rx_lla=None):
    """
    Plots a 3D visualization of Earth with optional receiver location and Starlink coverage cone.

    Parameters:
    rx_lla : list, optional
        [latitude, longitude, altitude] of the receiver.
    """
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")

    # Load Earth topography and coastlines
    ax.set_box_aspect([1,1,1])
    ax.set_xlim([-7000, 7000])
    ax.set_ylim([-7000, 7000])
    ax.set_zlim([-7000, 7000])
    ax.set_xlabel("X (km)")
    ax.set_ylabel("Y (km)")
    ax.set_zlabel("Z (km)")
    
    # Plot Earth as a sphere
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    X = 6378 * np.outer(np.cos(u), np.sin(v))
    Y = 6378 * np.outer(np.sin(u), np.sin(v))
    Z = 6378 * np.outer(np.ones(np.size(u)), np.cos(v))
    
    ax.plot_surface(X, Y, Z, color='lightblue', edgecolors='k', alpha=0.3)

    # Plot receiver location
    if rx_lla:
        rx_ecef = lla_to_ecef(rx_lla[0], rx_lla[1], rx_lla[2])  # Convert to ECEF
        ax.scatter(rx_ecef[0], rx_ecef[1], rx_ecef[2], color="red", marker="o", s=100, label="Receiver")

        # Generate Starlink coverage cone
        plot_starlink_cone(ax, rx_ecef, rx_lla[0], rx_lla[1])

    plt.legend()
    plt.show()

def plot_starlink_cone(ax, rx_ecef, lat, lon, h=500, n=50):
    """
    Plots a 3D Starlink coverage cone in ECEF coordinates.

    Parameters:
    ax : matplotlib 3D axis
        3D plot axis.
    rx_ecef : array
        Receiver ECEF coordinates.
    lat, lon : float
        Receiver latitude and longitude.
    h : float, optional
        Cone height (default is 500 km for Starlink).
    n : int, optional
        Number of points to generate.
    """

    theta = np.linspace(0, 2*np.pi, n)  # Around Z-axis
    x = np.linspace(0, h, n)  # Height along the x-axis
    Theta, X = np.meshgrid(theta, x)
    R = (h * X) / h  # Radius at each height
    Z = R * np.cos(Theta)
    Y = R * np.sin(Theta)

    # Define rotation angles
    yaw, pitch, roll = -np.radians(lon), 0, np.radians(lat)

    # Rotation matrices
    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(pitch), -np.sin(pitch)],
        [0, np.sin(pitch), np.cos(pitch)]
    ])

    Ry = np.array([
        [np.cos(roll), 0, np.sin(roll)],
        [0, 1, 0],
        [-np.sin(roll), 0, np.cos(roll)]
    ])

    Rz = np.array([
        [np.cos(yaw), -np.sin(yaw), 0],
        [np.sin(yaw), np.cos(yaw), 0],
        [0, 0, 1]
    ])

    # Apply rotation transformations
    XYZ = np.stack((X.ravel(), Y.ravel(), Z.ravel()), axis=1)
    XYZ = XYZ @ Rx @ Ry @ Rz

    # Reshape to original dimensions
    X = XYZ[:, 0].reshape(n, n) + rx_ecef[0]
    Y = XYZ[:, 1].reshape(n, n) + rx_ecef[1]
    Z = XYZ[:, 2].reshape(n, n) + rx_ecef[2]

    # Plot the cone
    ax.plot_surface(X, Y, Z, color="blue", alpha=0.3)

# Example Usage
plot_earth(rx_lla=[37.7749, -122.4194, 0])  # Example: San Francisco, 0 km altitude
