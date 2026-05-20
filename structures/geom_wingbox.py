from scipy.interpolate import interp1d
import numpy as np
from scipy.integrate import simpson
import matplotlib.pyplot as plt

def get_airfoil_heights(dat_file_path, xfs_pct, xrs_pct):
    # Load coordinates, ignoring the header line
    coords = np.loadtxt(dat_file_path, skiprows=1)
    x, y = coords[:, 0], coords[:, 1]

    # Find the leading edge
    le_idx = np.argmin(x)

    x_upper_raw, y_upper_raw = x[:le_idx + 1], y[:le_idx + 1]
    x_lower_raw, y_lower_raw = x[le_idx:], y[le_idx:]

    # Sort both surfaces by x to guarantee strictly increasing arrays
    idx_u = np.argsort(x_upper_raw)
    x_upper, y_upper = x_upper_raw[idx_u], y_upper_raw[idx_u]

    idx_l = np.argsort(x_lower_raw)
    x_lower, y_lower = x_lower_raw[idx_l], y_lower_raw[idx_l]

    # Find overlapping bounds to prevent extrapolation errors at the trailing edge
    min_x = max(x_upper.min(), x_lower.min())
    max_x = min(x_upper.max(), x_lower.max())

    xfs_pct = np.clip(xfs_pct, min_x, max_x)
    xrs_pct = np.clip(xrs_pct, min_x, max_x)

    f_upper = interp1d(x_upper, y_upper, kind='linear', bounds_error=False, fill_value="extrapolate")
    f_lower = interp1d(x_lower, y_lower, kind='linear', bounds_error=False, fill_value="extrapolate")

    hfs_norm = float(f_upper(xfs_pct) - f_lower(xfs_pct))
    hrs_norm = float(f_upper(xrs_pct) - f_lower(xrs_pct))
    return hfs_norm, hrs_norm


def calculate_section_inertia(tskin, tspar_web, wspar_cap, tspar_cap, A_str, n_str, xfs, xrs, hfs, hrs, materials, E_ref, G_ref):
    if n_str % 2 != 0:
        raise ValueError("n_str must be even to maintain chordwise symmetry.")
    # This is a SYMMETRIC trapezoid hence y_bar = 0

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
    Astr_eq = A_str * n_str_mod

    # Equivalent thicknesses for torsion
    tskin_G_eq = tskin * nG_skin

    tspar_web_eq = tspar_web * n_spar
    A_spar_cap_eq = wspar_cap * tspar_cap * n_spar
    tspar_web_G_eq = tspar_web * nG_spar

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
    A_fs_eq = hfs * tspar_web_eq + 2 * A_spar_cap_eq
    A_rs_eq = hrs * tspar_web_eq + 2 * A_spar_cap_eq
    A_webs = A_fs_eq + A_rs_eq
    A_skins = 2 * (slope_len * tskin_eq)
    Total_Area = A_webs + A_skins + Astr_total

    # --- 4. Centroid (x_bar) using equivalent areas ---
    mx_webs = (xfs * A_fs_eq) + (xrs * A_rs_eq)
    mx_skins = 2 * (tskin_eq * slope_len * (xfs + xrs) / 2)
    mx_str = 2 * np.sum(x_str * Astr_eq)
    x_bar = (mx_webs + mx_skins + mx_str) / Total_Area

    # --- 5. Ixx (Vertical Bending Stiffness) ---
    Ixx_webs_only = (1 / 12) * tspar_web_eq * (hfs ** 3 + hrs ** 3)
    Ixx_caps = 2 * A_spar_cap_eq * ((hfs / 2) ** 2 + (hrs / 2) ** 2)
    Ixx_webs = Ixx_webs_only + Ixx_caps
    Ixx_skins = 2 * (tskin_eq * slope_len / 12) * (hfs ** 2 + hfs * hrs + hrs ** 2)
    Ixx_str = np.sum(Astr_eq * y_top_str ** 2) + np.sum(Astr_eq * y_bot_str ** 2)
    Ixx_eq_total = Ixx_webs + Ixx_skins + Ixx_str

    # --- 6. Izz (Chordwise Bending Stiffness) ---
    Izz_webs_steiner = (A_fs_eq * (xfs - x_bar) ** 2) + (A_rs_eq * (xrs - x_bar) ** 2)
    Izz_caps_intrinsic = 4 * ((1 / 12) * (wspar_cap ** 3) * tspar_cap * n_spar)
    Izz_webs = Izz_webs_steiner + Izz_caps_intrinsic
    Izz_skins = 2 * tskin_eq * slope_len * ((xfs ** 2 + xfs * xrs + xrs ** 2) / 3 - (xfs + xrs) * x_bar + x_bar ** 2)
    Izz_str = 2 * np.sum(Astr_eq * (x_str - x_bar) ** 2)
    Izz_eq_total = Izz_webs + Izz_skins + Izz_str

    # --- 7. Torsional Constant (J) - Bredt-Batho with Shear Modulus weighting ---
    Ae = 0.5 * (hfs + hrs) * w
    perimeter_integral = (hfs / tspar_web_G_eq) + (hrs / tspar_web_G_eq) + (2 * slope_len / tskin_G_eq)
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

def compute_wingbox_mass(y_stations, chords, airfoil_file, xfs_pct, xrs_pct, tskin,
    tspar_web, wspar_cap, tspar_cap, A_str, n_str, materials):
    hfs_norm, hrs_norm = get_airfoil_heights(airfoil_file, xfs_pct, xrs_pct)

    # Calculate physical dimensions at every station
    hfs = hfs_norm * chords
    hrs = hrs_norm * chords
    w = (xrs_pct - xfs_pct) * chords
    theta = np.arctan(0.5 * (hfs - hrs) / w)
    slope_len = w / np.cos(theta)

    # Calculate cross-sectional areas (m^2)
    A_webs = (hfs + hrs) * tspar_web + 4 * (wspar_cap * tspar_cap)
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


def compute_wingbox_properties(chords, airfoil_file, xfs_pct, xrs_pct, tskin_m, tspar_web, wspar_cap, tspar_cap, Astr_one_m, n_str, materials, E_ref, G_ref):
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
            tskin_m, tspar_web, wspar_cap, tspar_cap, Astr_one_m, n_str,
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


if "__main__" == __name__:
    # --- PLOTTING & VERIFICATION ---
    file_path = "onze_airfoil.dat"
    front_spar_x = 0.15
    rear_spar_x = 0.55

    hfs, hrs = get_airfoil_heights(file_path, front_spar_x, rear_spar_x)

    coords = np.loadtxt(file_path, skiprows=1)
    x, y = coords[:, 0], coords[:, 1]
    le_idx = np.argmin(x)

    x_upper_raw, y_upper_raw = x[:le_idx + 1], y[:le_idx + 1]
    x_lower_raw, y_lower_raw = x[le_idx:], y[le_idx:]

    idx_u = np.argsort(x_upper_raw)
    x_up, y_up = x_upper_raw[idx_u], y_upper_raw[idx_u]

    idx_l = np.argsort(x_lower_raw)
    x_lo, y_lo = x_lower_raw[idx_l], y_lower_raw[idx_l]

    f_upper = interp1d(x_up, y_up, kind='linear', bounds_error=False, fill_value="extrapolate")
    f_lower = interp1d(x_lo, y_lo, kind='linear', bounds_error=False, fill_value="extrapolate")

    # Sample the interpolated functions to draw smooth surfaces
    x_dense = np.linspace(0, 1, 200)
    y_up_dense = f_upper(x_dense)
    y_lo_dense = f_lower(x_dense)

    # Get discrete spar coordinates for drawing lines
    y_fs_up, y_fs_lo = float(f_upper(front_spar_x)), float(f_lower(front_spar_x))
    y_rs_up, y_rs_lo = float(f_upper(rear_spar_x)), float(f_lower(rear_spar_x))

    # Setup matplotlib figure
    plt.figure(figsize=(12, 5))

    # Plot raw upper and lower coordinates
    plt.scatter(x_up, y_up, color='gray', alpha=0.5, label='Upper Raw Data', zorder=2)
    plt.scatter(x_lo, y_lo, color='lightgray', alpha=0.5, label='Lower Raw Data', zorder=2)

    # Plot continuous interpolated profiles
    plt.plot(x_dense, y_up_dense, 'b-', linewidth=1.5, label='Interpolated Upper')
    plt.plot(x_dense, y_lo_dense, 'r-', linewidth=1.5, label='Interpolated Lower')

    # Draw the vertical spar lines using heights generated by your function
    plt.vlines(x=front_spar_x, ymin=y_fs_lo, ymax=y_fs_up, colors='g', linestyles='--', linewidth=2,
               label=f'Front Spar (Height: {hfs:.4f})', zorder=3)
    plt.vlines(x=rear_spar_x, ymin=y_rs_lo, ymax=y_rs_up, colors='m', linestyles='--', linewidth=2,
               label=f'Rear Spar (Height: {hrs:.4f})', zorder=3)

    # Mark intersection points
    plt.scatter([front_spar_x, front_spar_x, rear_spar_x, rear_spar_x],
                [y_fs_up, y_fs_lo, y_rs_up, y_rs_lo], color='black', edgecolor='white', s=50, zorder=4)

    plt.xlabel("Chord (x/c)", fontsize=11)
    plt.ylabel("Thickness (y/c)", fontsize=11)
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(loc='upper right')

    # Keep aspect ratio 1:1 so profile isn't distorted
    plt.axis('equal')

    plt.show()
