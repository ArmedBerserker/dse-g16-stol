import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import simpson, cumulative_trapezoid

# --- Constants & Geometry Settings ---
xfs_pct = 0.15
xrs_pct = 0.65

tskin_m = 0.001
tspar_m = 0.003
Astr_m = 0.001

# --- Wing Geometry Inputs ---
S = 23.0  # Surface Area (m^2)
AR = 9.0  # Aspect Ratio
taper = 0.4  # Taper Ratio (tip_chord / root_chord)

# Calculate semi-span and chords
b = np.sqrt(S * AR)  # Total span (m)
y_tip = b / 2  # Semi-span (m)
c_root = (2 * S) / (b * (1 + taper))
c_tip = taper * c_root

# --- Load Internal Loads Data ---
# CSV columns: y, Vz, Mx, Vy, Mz, T
# skiprows=1 handles the "y_m" header error
data = np.loadtxt("internal_loads.csv", delimiter=",", skiprows=1)

y_stations = data[:, 0]      # Column 0: Spanwise position (y)
M_stations = data[:, 2]      # Column 2: Bending Moment (Mx)
T_stations = data[:, 5]      # Column 5: Torsion (T)

# Calculate chords at the specific CSV stations
chords = c_root - (c_root - c_tip) * (y_stations / y_tip)

# --- Your Functions (Untouched) ---
def get_airfoil_heights(dat_file_path, xfs_pct, xrs_pct):
    # Load coordinates, skipping the header line
    coords = np.loadtxt(dat_file_path, skiprows=1)
    x, y = coords[:, 0], coords[:, 1]

    # Split into Upper and Lower surfaces
    le_idx = np.argmin(x)

    # Ensure x is strictly increasing for interpolation
    x_upper, y_upper = x[:le_idx + 1][::-1], y[:le_idx + 1][::-1]
    x_lower, y_lower = x[le_idx:], y[le_idx:]

    f_upper = interp1d(x_upper, y_upper, kind='linear', bounds_error=False, fill_value="extrapolate")
    f_lower = interp1d(x_lower, y_lower, kind='linear', bounds_error=False, fill_value="extrapolate")

    # Normalized heights (thickness/chord)
    hfs_norm = float(f_upper(xfs_pct) - f_lower(xfs_pct))
    hrs_norm = float(f_upper(xrs_pct) - f_lower(xrs_pct))

    return hfs_norm, hrs_norm

def calculate_section_inertia(tskin, tspar, Astr, xfs, xrs, hfs, hrs):
    # Width of the box
    w = xrs - xfs

    # Ixx Calculation
    I_spars = (1 / 12) * tspar * (hfs ** 3 + hrs ** 3)

    # Consider the sloped skins on top and bottom and take integral of y^2
    int_factor = (w / 12) * (hfs ** 2 + hfs * hrs + hrs ** 2)

    I_skins = 2 * tskin * int_factor
    I_stringers = 2 * (Astr / w) * int_factor

    Ixx = I_spars + I_skins + I_stringers

    # J Calculation
    # Enclosed Area (Trapezoid)
    Ae = 0.5 * (hfs + hrs) * w

    # Length of the sloping skins (Pythagoras)
    slope_length = np.sqrt(w ** 2 + (0.5 * (hfs - hrs)) ** 2)

    # Perimeter Integral (sum of length/thickness for each wall)
    integral = (hfs / tspar) + (hrs / tspar) + (2 * slope_length / tskin)

    J = (4 * Ae ** 2) / integral  # Bredt-Batho

    return Ixx, J

def calculate_torsional_deflection(y, T, G, J):
    # Twist rate (d_theta/dy)
    twist_rate = T / (G * J)

    # Integrate twist rate along the span using Simpson's rule
    total_twist_rad = simpson(twist_rate, y)

    return np.rad2deg(total_twist_rad)

def calculate_bending_deflection(y, M, E, Ixx):
    # Curvature (d2v/dy2)
    curvature = M / (E * Ixx)

    # First integration: Slope (phi)
    slope = cumulative_trapezoid(curvature, y, initial=0)

    # Second integration: Deflection (v)
    tip_deflection_m = simpson(y=slope, x=y)

    return tip_deflection_m * 1000  # Return in mm



# --- Calculations ---

hfs_norm, hrs_norm = get_airfoil_heights("onze_airfoil.dat", xfs_pct, xrs_pct)

Ixx_stations = []
J_stations = []

for c in chords:
    # Scale geometry for the local chord
    local_hfs = hfs_norm * c
    local_hrs = hrs_norm * c
    local_xfs = xfs_pct * c
    local_xrs = xrs_pct * c

    # Calculate inertia for this specific station
    I_val, J_val = calculate_section_inertia(
        tskin_m, tspar_m, Astr_m,
        local_xfs, local_xrs, local_hfs, local_hrs
    )

    Ixx_stations.append(I_val)
    J_stations.append(J_val)

Ixx_stations = np.array(Ixx_stations)
J_stations = np.array(J_stations)

# Final Deflection Results
# Aluminum Properties: E = 70 GPa, G = 26 GPa
twist = calculate_torsional_deflection(y_stations, T_stations, 26e9, J_stations)
bending = calculate_bending_deflection(y_stations, M_stations, 70e9, Ixx_stations)

print(f"Total Tip Twist: {twist:.4f} deg")
print(f"Total Tip Bending: {bending:.4f} mm")