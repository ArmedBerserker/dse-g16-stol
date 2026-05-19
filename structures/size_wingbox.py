from geom_wingbox import *
from wingbox_helpers import *
import numpy as np

# --- Constants & Geometry Settings ---
load_factor = 3.8
g_accel = 9.81
xfs_pct = 0.15  # 15-20
xrs_pct = 0.65  # 55-65
tskin_m = 0.001
tspar_m = 0.003
Astr_one_m = 0.0001  # Cross-sectional area of ONE stringer (m^2)
n_str = 16  # Total number of stringers (MUST be even)

# --- Material Properties ---
# E = Young's Modulus (Pa), G = Shear Modulus (Pa), rho = Density (kg/m^3)
materials = {
    'skin': {'E': 73.1e9, 'G': 28.0e9, 'rho': 2780},
    'spar': {'E': 71.7e9, 'G': 26.9e9, 'rho': 2810},
    'str':  {'E': 71.7e9, 'G': 26.9e9, 'rho': 2810}
}

# Define the "Reference Material" to normalize the cross-section
E_ref = materials['spar']['E']
G_ref = materials['spar']['G']

# Wing Geometry
S = 25.675  # Surface Area (m^2)
AR = 9.0  # Aspect Ratio
taper = 0.4  # Taper Ratio

# Engine parameters
y_eng = 2.19 # Spanwise position of engine from root (m)
m_eng = 120  # Mass of the engine (kg)
thrust_eng = 1890  # Forward engine thrust force (N)
dx_eng = 1.0  # Forward coordinate relative to shear center/centroid (m, positive forward)
dz_eng = -0.375  # Vertical coordinate relative to shear center/centroid (m, positive up)
Fn_eng = m_eng * g_accel  # Weight force acting DOWN (Positive in local frame)
Ft_eng = -thrust_eng  # Thrust force acting FORWARD (Negative in local frame)


def compute_internal_loads(y_stations, chords, T_c4, Fn_aero, Ft_aero, m_prime, x_sc_stations, n):
    # 1. Calculate span (dy) for mass to point load conversion
    dy = np.zeros_like(y_stations)
    if len(y_stations) > 1:
        dy[1:-1] = (y_stations[2:] - y_stations[:-2]) / 2
        dy[0] = (y_stations[1] - y_stations[0]) / 2
        dy[-1] = (y_stations[-1] - y_stations[-2]) / 2
    Fn_weight = m_prime * g_accel * dy

    # 2. Aerodynamic Torsion Shift
    dx_sc = x_sc_stations - (0.25 * chords)
    T_sc_panel = T_c4 - (dx_sc * Fn_aero)

    # 3. Combine Forces for Vertical Shear and Bending
    Fn_total = Fn_aero + Fn_weight

    # 4. Shear and Torsion
    V_stations = np.cumsum(Fn_total[::-1])[::-1]
    T_stations = np.cumsum(T_sc_panel[::-1])[::-1]

    # 5. Bending Moments
    dy_matrix = y_stations[None, :] - y_stations[:, None]
    dy_matrix[dy_matrix < 0] = 0  # Zero out inboard contributions

    Mx_stations = np.sum(Fn_total * dy_matrix, axis=1)
    Mz_stations = np.sum(Ft_aero * dy_matrix, axis=1)

    # 6. Add Engine Point Load to All Inboard Stations
    inboard_mask = y_stations <= y_eng
    lever_arm_eng = y_eng - y_stations[inboard_mask]

    V_stations[inboard_mask] += Fn_eng
    Mx_stations[inboard_mask] += Fn_eng * lever_arm_eng
    Mz_stations[inboard_mask] += Ft_eng * lever_arm_eng
    T_stations[inboard_mask] += -(Fn_eng * dx_eng) + (Ft_eng * dz_eng)

    V_stations, Mx_stations, Mz_stations, T_stations = (V_stations * n,
    Mx_stations * n, Mz_stations *n , T_stations * n)

    return V_stations, Mx_stations, Mz_stations, T_stations


def main():
    # 1. Inputs
    b, y_tip, c_root, c_tip, C_mac = calculate_wing_geometry(S, AR, taper)

    # Rename variables to clarify they are aerodynamic-only
    y_stations, T_c4, Fn_aero, Ft_aero, xflr5_chords = load_aerodynamic_data("MainWing_a=3.50_v=72.00ms.txt")

    analytical_chords = c_root - (c_root - c_tip) * (y_stations / y_tip)

    # Chord length check
    max_error = np.max(np.abs(xflr5_chords - analytical_chords))
    if max_error > 0.01:
        print(f"WARNING: XFLR5 geometry deviates from analytical by max {max_error:.4f} m!")

    chords = analytical_chords

    # 2. Structural Properties
    Ixx_eq, Izz_eq, J_eq, x_sc_stations = compute_wingbox_properties(
        chords, "onze_airfoil.dat", xfs_pct, xrs_pct, tskin_m, tspar_m, Astr_one_m, n_str,
        materials, E_ref, G_ref
    )

    # 3. Mass Calculation (Must happen before internal loads to get m_prime)
    half_wing_mass, mass_webs, mass_skins, mass_str, m_prime = compute_wingbox_mass(
        y_stations, chords, "onze_airfoil.dat", xfs_pct, xrs_pct,
        tskin_m, tspar_m, Astr_one_m, n_str, materials
    )

    # 4. Internal Load Distribution
    V_stations, Mx_stations, Mz_stations, T_stations = compute_internal_loads(
        y_stations, chords, T_c4, Fn_aero, Ft_aero, m_prime, x_sc_stations, load_factor
    )

    # 5. Deflection Integration Using Reference Modulus
    twist_val = calculate_torsional_deflection(y_stations, T_stations, G_ref, J_eq)[-1]
    bending_n_val = calculate_bending_deflection(y_stations, Mx_stations, E_ref, Ixx_eq)[-1]
    bending_t_val = calculate_bending_deflection(y_stations, Mz_stations, E_ref, Izz_eq)[-1]

    # 6. Output Processing
    display_bending_normal = -1 * bending_n_val  # Up = Positive
    display_bending_tangential = -1 * bending_t_val  # Forward = Positive
    display_twist = twist_val  # Nose-Up = Positive
    total_mass = half_wing_mass * 2

    print("--- STRUCTURAL MASS ---")
    print(f"Skin Mass:      {mass_skins:.2f} kg")
    print(f"Spars Mass:     {mass_webs:.2f} kg")
    print(f"Stringers Mass: {mass_str:.2f} kg")
    print(f"TOTAL AIRCRAFT WINGBOX MASS: {total_mass:.2f} kg\n")

    print("--- DEFLECTIONS (Positive = Up / Forward / Nose-Up) ---")
    print(f"Total Tip Twist: {display_twist:.4f} deg")
    print(f"Total Tip Vertical Deflection: {display_bending_normal:.4f} mm")
    print(f"Total Tip Chordwise Deflection: {display_bending_tangential:.4f} mm")

    return (display_twist, display_bending_normal, display_bending_tangential, total_mass)


if __name__ == "__main__":
    main()