from scipy.interpolate import interp1d
import numpy as np
from scipy.integrate import simpson

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


def calculate_section_inertia(tskin, tspar, A_str, n_str, xfs, xrs, hfs, hrs, materials, E_ref, G_ref):
    if n_str % 2 != 0:
        raise ValueError("n_str must be even to maintain chordwise symmetry.")

    # --- 1. Material Modular Ratios ---
    # Bending
    n_skin = materials['skin']['E'] / E_ref
    n_spar = materials['spar']['E'] / E_ref
    n_str_mod = materials['str']['E'] / E_ref

    # Torsion
    nG_skin = materials['skin']['G'] / G_ref
    nG_spar = materials['spar']['G'] / G_ref

    # Equivalent thicknesses and areas for bending
    tskin_eq = tskin * n_skin
    tspar_eq = tspar * n_spar
    Astr_eq = A_str * n_str_mod

    # Equivalent thicknesses for torsion
    tskin_G_eq = tskin * nG_skin
    tspar_G_eq = tspar * nG_spar

    # --- 2. Geometric Constants ---
    w = xrs - xfs
    theta = np.arctan(0.5 * (hfs - hrs) / w)
    slope_len = w / np.cos(theta)
    y_mid = 0.0  # Vertical datum

    n_per_side = n_str // 2

    # Stringer Placement (Placed along the sloped skins between spars)
    x_str = np.linspace(xfs, xrs, n_per_side + 2)[1:-1]
    y_top_str = y_mid + 0.5 * (hfs + (x_str - xfs) * (hrs - hfs) / w)
    y_bot_str = y_mid - 0.5 * (hfs + (x_str - xfs) * (hrs - hfs) / w)

    # --- 3. Equivalent Areas ---
    Astr_total = Astr_eq * n_str
    A_webs = (hfs + hrs) * tspar_eq
    A_skins = 2 * (slope_len * tskin_eq)
    Total_Area = A_webs + A_skins + Astr_total

    # --- 4. Centroid (x_bar) using equivalent areas ---
    mx_webs = (xfs * hfs * tspar_eq) + (xrs * hrs * tspar_eq)
    mx_skins = 2 * (tskin_eq * slope_len * (xfs + xrs) / 2)
    mx_str = 2 * np.sum(x_str * Astr_eq)
    x_bar = (mx_webs + mx_skins + mx_str) / Total_Area

    # --- 5. Ixx (Vertical Bending Stiffness) ---
    Ixx_webs = (1 / 12) * tspar_eq * (hfs ** 3 + hrs ** 3)
    Ixx_skins = 2 * (tskin_eq * slope_len / 12) * (hfs ** 2 + hfs * hrs + hrs ** 2)
    Ixx_str = np.sum(Astr_eq * y_top_str ** 2) + np.sum(Astr_eq * y_bot_str ** 2)
    Ixx_eq_total = Ixx_webs + Ixx_skins + Ixx_str

    # --- 6. Izz (Chordwise Bending Stiffness) ---
    Izz_webs = (hfs * tspar_eq * (xfs - x_bar) ** 2) + (hrs * tspar_eq * (xrs - x_bar) ** 2)
    Izz_skins = 2 * tskin_eq * slope_len * ((xfs ** 2 + xfs * xrs + xrs ** 2) / 3 - (xfs + xrs) * x_bar + x_bar ** 2)
    Izz_str = 2 * np.sum(Astr_eq * (x_str - x_bar) ** 2)
    Izz_eq_total = Izz_webs + Izz_skins + Izz_str

    # --- 7. Torsional Constant (J) - Bredt-Batho with Shear Modulus weighting ---
    Ae = 0.5 * (hfs + hrs) * w
    perimeter_integral = (hfs / tspar_G_eq) + (hrs / tspar_G_eq) + (2 * slope_len / tskin_G_eq)
    J_eq_total = (4 * Ae ** 2) / perimeter_integral

    # Shear Center Approx
    x_sc_approx = 0.5 * (xfs + xrs)

    return {
        "Ixx": Ixx_eq_total,
        "Izz": Izz_eq_total,
        "J": J_eq_total,
        "x_bar": x_bar,
        "x_sc_midpoint": x_sc_approx
    }


def compute_wingbox_properties(chords, airfoil_file, xfs_pct, xrs_pct, tskin_m, tspar_m, Astr_one_m, n_str, materials, E_ref, G_ref):
    # Computes equivalent inertia properties across all spanwise stations
    hfs_norm, hrs_norm = get_airfoil_heights(airfoil_file, xfs_pct, xrs_pct)

    Ixx_list = []
    Izz_list = []
    J_list = []
    x_sc_list = []

    for c in chords:
        local_hfs = hfs_norm * c
        local_hrs = hrs_norm * c
        local_xfs = xfs_pct * c
        local_xrs = xrs_pct * c

        section_props = calculate_section_inertia(
            tskin_m, tspar_m, Astr_one_m, n_str,
            local_xfs, local_xrs, local_hfs, local_hrs,
            materials, E_ref, G_ref
        )

        Ixx_list.append(section_props["Ixx"])
        Izz_list.append(section_props["Izz"])
        J_list.append(section_props["J"])
        x_sc_list.append(section_props["x_sc_midpoint"])

    return (
        np.array(Ixx_list),
        np.array(Izz_list),
        np.array(J_list),
        np.array(x_sc_list)
    )

def compute_wingbox_mass(y_stations, chords, airfoil_file, xfs_pct, xrs_pct, tskin, tspar, A_str, n_str, materials):
    hfs_norm, hrs_norm = get_airfoil_heights(airfoil_file, xfs_pct, xrs_pct)

    # Calculate physical dimensions at every station
    hfs = hfs_norm * chords
    hrs = hrs_norm * chords
    w = (xrs_pct - xfs_pct) * chords
    theta = np.arctan(0.5 * (hfs - hrs) / w)
    slope_len = w / np.cos(theta)

    # Calculate cross-sectional areas (m^2)
    A_webs = (hfs + hrs) * tspar
    A_skins = 2 * (slope_len * tskin)
    A_str_total = np.full_like(chords, A_str * n_str)

    # Calculate distributed mass profile (kg/m)
    m_prime = (A_webs * materials['spar']['rho'] +
               A_skins * materials['skin']['rho'] +
               A_str_total * materials['str']['rho'])

    # Integrate area over the span to get Volume (m^3)
    vol_webs = simpson(A_webs, x=y_stations)
    vol_skins = simpson(A_skins, x=y_stations)
    vol_str = simpson(A_str_total, x=y_stations)

    # Multiply by respective material densities (kg)
    mass_webs = vol_webs * materials['spar']['rho']
    mass_skins = vol_skins * materials['skin']['rho']
    mass_str = vol_str * materials['str']['rho']

    total_half_mass = mass_webs + mass_skins + mass_str

    return total_half_mass, mass_webs, mass_skins, mass_str, m_prime