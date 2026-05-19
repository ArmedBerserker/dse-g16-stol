from geom_wingbox import *
from wingbox_helpers import *
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

# --- Constants & Geometry Settings ---
load_factor = 1.0

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
Fn_eng = m_eng * 9.81  # Weight force acting DOWN (Positive in local frame)
Ft_eng = -thrust_eng  # Thrust force acting FORWARD (Negative in local frame)

@dataclass
class WingResults:
    total_mass: float
    load_factor: float
    mass_webs: float
    mass_skins: float
    mass_str: float
    bending_n_val: np.ndarray
    bending_t_val: np.ndarray
    twist_val: np.ndarray
    y_stations: np.ndarray


def main(
    S, AR, taper,
    xflr5_file, airfoil_file, xfs_pct, xrs_pct, tskin_m, tspar_m, Astr_one_m, n_str,
    materials, E_ref, G_ref, load_factor, y_eng, Fn_eng, Ft_eng, dz_eng, dx_eng,
):
    # 1. Inputs
    b, y_tip, c_root, c_tip, C_mac = calculate_wing_geometry(S, AR, taper)

    # Rename variables to clarify they are aerodynamic-only
    y_stations, T_c4, Fn_aero, Ft_aero, xflr5_chords = load_aerodynamic_data(xflr5_file)

    analytical_chords = c_root - (c_root - c_tip) * (y_stations / y_tip)

    # Chord length check
    max_error = np.max(np.abs(xflr5_chords - analytical_chords))
    if max_error > 0.01:
        print(f"WARNING: XFLR5 geometry deviates from analytical by max {max_error:.4f} m!")

    chords = analytical_chords

    # 2. Structural Properties
    Ixx_eq, Izz_eq, J_eq, x_sc_stations = compute_wingbox_properties(
        chords, airfoil_file, xfs_pct, xrs_pct, tskin_m, tspar_m, Astr_one_m, n_str,
        materials, E_ref, G_ref
    )

    # 3. Mass Calculation (Must happen before internal loads to get m_prime)
    half_wing_mass, mass_webs, mass_skins, mass_str, m_prime = compute_wingbox_mass(
        y_stations, chords, airfoil_file, xfs_pct, xrs_pct,
        tskin_m, tspar_m, Astr_one_m, n_str, materials
    )

    # 4. Internal Load Distribution
    V_stations, Mx_stations, Mz_stations, T_stations = compute_internal_loads(
        y_stations, chords, T_c4, Fn_aero, Ft_aero, m_prime, x_sc_stations,
        load_factor, y_eng, Fn_eng, Ft_eng, dz_eng, dx_eng
    )

    # 5. Deflection Integration Using Reference Modulus
    twist_val = calculate_torsional_deflection(y_stations, T_stations, G_ref, J_eq)
    bending_n_val = calculate_bending_deflection(y_stations, Mx_stations, E_ref, Ixx_eq)
    bending_t_val = calculate_bending_deflection(y_stations, Mz_stations, E_ref, Izz_eq)

    total_mass = half_wing_mass * 2

    return WingResults(
        total_mass=total_mass,
        load_factor=load_factor,

        mass_webs=mass_webs,
        mass_skins=mass_skins,
        mass_str=mass_str,

        bending_n_val=bending_n_val,
        bending_t_val=bending_t_val,
        twist_val=twist_val,

        y_stations=y_stations
    )

def output_results(results):
    bending_normal = -results.bending_n_val
    bending_tangential = -results.bending_t_val
    twist = results.twist_val

    print("--- STRUCTURAL MASS ---")
    print(f"Skin Mass:      {results.mass_skins:.2f} kg")
    print(f"Spars Mass:     {results.mass_webs:.2f} kg")
    print(f"Stringers Mass: {results.mass_str:.2f} kg")
    print(f"TOTAL AIRCRAFT WINGBOX MASS: {results.total_mass:.2f} kg\n")

    print("--- DEFLECTIONS (Positive = Up / Forward / Nose-Up) ---")
    print(f"Load Factor: {results.load_factor}")
    print(f"Total Tip Twist: {twist[-1]:.4f} deg")
    print(f"Total Tip Vertical Deflection: {bending_normal[-1]:.4f} mm")
    print(f"Total Tip Chordwise Deflection: {bending_tangential[-1]:.4f} mm")

    fig, axes = plt.subplots(2, 1, figsize=(8, 8))

    # Vertical bending plot
    axes[0].plot(results.y_stations, bending_normal, linewidth=2)
    axes[0].set_title("Vertical Bending Deflection")
    axes[0].set_xlabel("Spanwise Position y [m]")
    axes[0].set_ylabel("Vertical Deflection [mm]")
    axes[0].grid(True)

    # Twist angle plot
    axes[1].plot(results.y_stations, twist, linewidth=2, color='orange')
    axes[1].set_title("Twist Angle Distribution")
    axes[1].set_xlabel("Spanwise Position y [m]")
    axes[1].set_ylabel("Twist Angle [deg]")
    axes[1].grid(True)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    results = main(
        S=S,
        AR=AR,
        taper=taper,
        xflr5_file="MainWing_a=3.50_v=72.00ms.txt",
        airfoil_file="onze_airfoil.dat",
        xfs_pct=xfs_pct,
        xrs_pct=xrs_pct,
        tskin_m=tskin_m,
        tspar_m=tspar_m,
        Astr_one_m=Astr_one_m,
        n_str=n_str,
        materials=materials,
        E_ref=E_ref,
        G_ref=G_ref,
        load_factor=load_factor,
        y_eng=y_eng,
        Fn_eng=Fn_eng,
        Ft_eng=Ft_eng,
        dz_eng=dz_eng,
        dx_eng=dx_eng
    )

    output_results(results)