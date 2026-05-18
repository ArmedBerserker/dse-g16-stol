import numpy as np
from scipy.interpolate import interp1d


def get_airfoil_heights(dat_file_path, xfs_pct, xrs_pct):
    coords = np.loadtxt(dat_file_path, skiprows=1)
    x, y = coords[:, 0], coords[:, 1]

    le_idx = np.argmin(x)
    # Reverse upper surface to ensure strictly increasing x for interpolation
    x_upper, y_upper = x[:le_idx + 1][::-1], y[:le_idx + 1][::-1]
    x_lower, y_lower = x[le_idx:], y[le_idx:]

    # Enforce bounds to prevent crazy behavior
    xfs_pct = np.clip(xfs_pct, x.min(), x.max())
    xrs_pct = np.clip(xrs_pct, x.min(), x.max())

    f_upper = interp1d(x_upper, y_upper, kind='linear', bounds_error=False, fill_value="extrapolate")
    f_lower = interp1d(x_lower, y_lower, kind='linear', bounds_error=False, fill_value="extrapolate")

    hfs_norm = float(f_upper(xfs_pct) - f_lower(xfs_pct))
    hrs_norm = float(f_upper(xrs_pct) - f_lower(xrs_pct))

    return hfs_norm, hrs_norm

#print(get_airfoil_heights("onze_airfoil.dat", 0.15, 0.65))

def calculate_section_inertia(tskin, tspar, A_str, n_str, xfs, xrs, hfs, hrs):
    if n_str % 2 != 0:
        raise ValueError("n_str must be even to maintain chordwise symmetry.")

    # Geometric Constants
    w = xrs - xfs
    theta = np.arctan(0.5 * (hfs - hrs) / w)
    slope_len = w / np.cos(theta)
    y_mid = 0.0  # Vertical datum

    n_per_side = n_str // 2

    # Stringer Placement
    # Placed along the sloped skins between spars
    x_str = np.linspace(xfs, xrs, n_per_side + 2)[1:-1]
    y_top_str = y_mid + 0.5 * (hfs + (x_str - xfs) * (hrs - hfs) / w)
    y_bot_str = y_mid - 0.5 * (hfs + (x_str - xfs) * (hrs - hfs) / w)

    # Areas
    Astr_total = A_str * n_str
    A_webs = (hfs + hrs) * tspar
    A_skins = 2 * (slope_len * tskin)
    Total_Area = A_webs + A_skins + Astr_total

    # Centroid (x_bar)
    mx_webs = (xfs * hfs * tspar) + (xrs * hrs * tspar)
    mx_skins = 2 * (tskin * slope_len * (xfs + xrs) / 2)
    mx_str = 2 * np.sum(x_str * A_str)
    x_bar = (mx_webs + mx_skins + mx_str) / Total_Area

    # Ixx (Vertical Bending Stiffness)
    Ixx_webs = (1 / 12) * tspar * (hfs ** 3 + hrs ** 3)
    # Integral for skins at +/- h/2
    Ixx_skins = (tskin * slope_len / 12) * (hfs ** 2 + hfs * hrs + hrs ** 2)
    Ixx_str = np.sum(A_str * y_top_str ** 2) + np.sum(A_str * y_bot_str ** 2)
    Ixx = Ixx_webs + Ixx_skins + Ixx_str

    # Izz (Chordwise Bending Stiffness)
    Izz_webs = (hfs * tspar * (xfs - x_bar) ** 2) + (hrs * tspar * (xrs - x_bar) ** 2)

    Izz_skins = 2 * tskin * slope_len * ((xfs ** 2 + xfs * xrs + xrs ** 2) / 3 - (xfs + xrs) * x_bar + x_bar ** 2)

    Izz_str = 2 * np.sum(A_str * (x_str - x_bar) ** 2)
    Izz = Izz_webs + Izz_skins + Izz_str

    # Torsional Constant (J) - Bredt-Batho
    Ae = 0.5 * (hfs + hrs) * w
    perimeter_integral = (hfs / tspar) + (hrs / tspar) + (2 * slope_len / tskin)
    J = (4 * Ae ** 2) / perimeter_integral

    # Shear Center Approx
    x_sc_approx = 0.5 * (xfs + xrs)

    return {
        "Ixx": Ixx,
        "Izz": Izz,
        "J": J,
        "x_bar": x_bar,
        "x_sc_midpoint": x_sc_approx
    }