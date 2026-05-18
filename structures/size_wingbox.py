import numpy as np
from scipy.integrate import simpson, cumulative_trapezoid
from geom_wingbox import *

# --- Constants & Geometry Settings ---
xfs_pct = 0.15 # 15-20
xrs_pct = 0.65 # 55-65

tskin_m = 0.002
tspar_m = 0.003
Astr_one_m = 0.0005 # Cross-sectional area of ONE stringer (m^2)
n_str = 16  # Total number of stringers (MUST be even)

# Material Properties
E_modulus = 70e9  # Young's Modulus in Pa
G_modulus = 26e9  # Shear Modulus in Pa

# Wing Geometry
S = 26  # Surface Area (m^2)
AR = 9.0  # Aspect Ratio
taper = 0.4  # Taper Ratio (tip_chord / root_chord)

# Calculate semi-span and chords
b = np.sqrt(S * AR)  # Total span (m)
y_tip = b / 2  # Semi-span (m)
c_root = (2 * S) / (b * (1 + taper))
c_tip = taper * c_root




data = np.genfromtxt("xflr5_parsed.csv", delimiter=",", names=True)

# 1. Enforce Root-to-Tip Sorting
# This is required for integration starting at the root (y=0)
sort_idx = np.argsort(data['y'])

y_stations = data['y'][sort_idx]
chords = data['chord'][sort_idx]

# Global Forces (Discrete Panel Point Loads in N)
Fx = data['Fx'][sort_idx]
Fy = data['Fy'][sort_idx]
Fz = data['Fz'][sort_idx]

# Local Normal Unit Vectors
nx = data['nx'][sort_idx]
ny = data['ny'][sort_idx]
nz = data['nz'][sort_idx]

# Local Tangential Unit Vectors
cx = data['cx'][sort_idx]
cy = data['cy'][sort_idx]
cz = data['cz'][sort_idx]

# Total Panel Torsion at Quarter Chord (Discrete Panel Moment in N*m)
T_c4 = data['Torsion'][sort_idx]

# Transform Forces to Wingbox Local Frame
Fn_stations = (Fx * nx) + (Fy * ny) + (Fz * nz)  # Positive DOWN
Ft_stations = (Fx * cx) + (Fy * cy) + (Fz * cz)  # Positive AFT




def calculate_torsional_deflection(y, T, G, J):
    twist_rate = T / (G * J)
    total_twist_rad = simpson(twist_rate, x=y)
    return np.rad2deg(total_twist_rad)


def calculate_bending_deflection(y, M, E, I_inertia):
    curvature = M / (E * I_inertia)
    # Starts integration at index 0 (Root, y=0) so boundary condition slope=0 is respected
    slope = cumulative_trapezoid(curvature, x=y, initial=0)
    tip_deflection_m = simpson(slope, x=y)
    return tip_deflection_m * 1000  # Return in mm

def wingbox_volume(
    c_root,
    c_tip,
    y_tip,
    xfs_pct,
    xrs_pct,
    hfs_norm,
    hrs_norm,
    n=500
):

    y = np.linspace(0.0, y_tip, n)

    chord = c_root + (c_tip - c_root) * (y / y_tip)

    width = (xrs_pct - xfs_pct) * chord

    h_front = hfs_norm * chord
    h_rear  = hrs_norm * chord

    area = 0.5 * (h_front + h_rear) * width

    volume = simpson(area, y)

    return volume


# --- Main Calculations ---

hfs_norm, hrs_norm = get_airfoil_heights("onze_airfoil.dat", xfs_pct, xrs_pct)

Ixx_stations = []
Izz_stations = []
J_stations = []
x_sc_stations = []

for c in chords:
    local_hfs = hfs_norm * c
    local_hrs = hrs_norm * c
    local_xfs = xfs_pct * c
    local_xrs = xrs_pct * c

    section_props = calculate_section_inertia(
        tskin_m, tspar_m, Astr_one_m, n_str,
        local_xfs, local_xrs, local_hfs, local_hrs
    )

    Ixx_stations.append(section_props["Ixx"])
    Izz_stations.append(section_props["Izz"])
    J_stations.append(section_props["J"])
    x_sc_stations.append(section_props["x_sc_midpoint"])

Ixx_stations = np.array(Ixx_stations)
Izz_stations = np.array(Izz_stations)
J_stations = np.array(J_stations)
x_sc_stations = np.array(x_sc_stations)

# Move Torsion to Shear Center
dx_sc = x_sc_stations - (0.25 * chords)
T_sc_panel = T_c4 - (dx_sc * Fn_stations)

# Summation for Internal Loads (Tip to Root)
V_stations = np.zeros_like(y_stations)
Mx_stations = np.zeros_like(y_stations)
Mz_stations = np.zeros_like(y_stations)
T_stations = np.zeros_like(y_stations)

for i in range(len(y_stations)):
    # Find all panels strictly outboard (or at) the current station
    outboard_mask = y_stations >= y_stations[i]
    lever_arms = y_stations[outboard_mask] - y_stations[i]

    V_stations[i] = np.sum(Fn_stations[outboard_mask])
    Mx_stations[i] = np.sum(Fn_stations[outboard_mask] * lever_arms)
    Mz_stations[i] = np.sum(Ft_stations[outboard_mask] * lever_arms)
    T_stations[i] = np.sum(T_sc_panel[outboard_mask])

# Deflection Integration
twist_val = calculate_torsional_deflection(y_stations, T_stations, G_modulus, J_stations)
bending_n_val = calculate_bending_deflection(y_stations, Mx_stations, E_modulus, Ixx_stations)
bending_t_val = calculate_bending_deflection(y_stations, Mz_stations, E_modulus, Izz_stations)

# Transform integral signs to physical directions
display_bending_normal = -1 * bending_n_val  # Up = Positive
display_bending_tangential = -1 * bending_t_val  # Forward = Positive
display_twist = twist_val  # Nose-Up = Positive

print("--- DEFLECTIONS (Positive = Up / Forward / Nose-Up) ---")
print(f"Total Tip Twist: {display_twist:.4f} deg")
print(f"Total Tip Vertical Deflection: {display_bending_normal:.4f} mm")
print(f"Total Tip Chordwise Deflection: {display_bending_tangential:.4f} mm")


V_wingbox = wingbox_volume(
    c_root,
    c_tip,
    y_tip,
    xfs_pct,
    xrs_pct,
    hfs_norm,
    hrs_norm
)

print(V_wingbox)