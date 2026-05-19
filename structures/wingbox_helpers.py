import numpy as np
from scipy.integrate import cumulative_simpson


def calculate_wing_geometry(S, AR, taper):
    b = np.sqrt(S * AR)  # Total span (m)
    y_tip = b / 2  # Semi-span (m)
    c_root = (2 * S) / (b * (1 + taper))
    c_tip = taper * c_root
    C_mac = (2 / 3) * c_root * ((1 + taper + taper ** 2) / (1 + taper))
    return b, y_tip, c_root, c_tip, C_mac

def calculate_torsional_deflection(y, T, G, J):
    twist_rate = T / (G * J)
    # Integrate once to get twist angle over the span
    twist_rad_array = cumulative_simpson(twist_rate, x=y, initial=0)
    return np.rad2deg(twist_rad_array)

def calculate_bending_deflection(y, M, E, I_inertia):
    curvature = M / (E * I_inertia)
    # Integrate once to get slope
    slope = cumulative_simpson(curvature, x=y, initial=0)
    # Integrate again to get deflection over the span
    deflection_m_array = cumulative_simpson(slope, x=y, initial=0)
    return deflection_m_array * 1000  # Return array in mm

def load_aerodynamic_data(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    v_inf = 72.0  # Fallback default
    for line in lines:
        if line.startswith("QInf"):
            parts = line.replace('=', ' ').split()
            if len(parts) >= 2:
                v_inf = float(parts[1])
            break

    rho = 1.225  # kg/m^3
    q_inf = 0.5 * rho * (v_inf ** 2)

    # Locate the data block
    data_start = -1
    for i, line in enumerate(lines):
        if "y-span" in line and "Chord" in line and "Cl" in line:
            data_start = i + 1
            break

    if data_start == -1:
        raise ValueError("Could not find data table header in the input file.")

    # Parse spanwise coefficients
    y_list, c_list, cl_list, cd_list, cm_list = [], [], [], [], []
    for line in lines[data_start:]:
        parts = line.strip().split()

        if len(parts) < 12:
            if len(y_list) > 0:
                break
            continue

        try:
            y = float(parts[0])
            if y < 0:
                continue
            if y_list and abs(y - y_list[-1]) < 1e-5:
                continue

            chord = float(parts[1])
            cl = float(parts[3])
            pcd = float(parts[4])
            icd = float(parts[5])
            cm = float(parts[7])  # CmAirf@chord/4

            y_list.append(y)
            c_list.append(chord)
            cl_list.append(cl)
            cd_list.append(pcd + icd)
            cm_list.append(cm)
        except ValueError:
            pass

    # Convert to numpy arrays
    y_stations = np.array(y_list)
    chords = np.array(c_list)
    cls = np.array(cl_list)
    cds = np.array(cd_list)
    cms = np.array(cm_list)

    # Calculate Span (dy) for each station
    dy = np.zeros_like(y_stations)
    if len(y_stations) > 1:
        dy[1:-1] = (y_stations[2:] - y_stations[:-2]) / 2  # Interior points
        dy[0] = (y_stations[1] - y_stations[0]) / 2        # Root edge
        dy[-1] = (y_stations[-1] - y_stations[-2]) / 2     # Tip edge

    # Calculate final point loads per strip, small angle approx
    Fn_stations = -(q_inf * chords * cls) * dy  # Lift - Normal force (Positive Down)
    Ft_stations = (q_inf * chords * cds) * dy  # Drag - Tangential force (Positive Aft)
    T_c4 = (q_inf * (chords ** 2) * cms) * dy  # Pitching Moment - Torsion

    return y_stations, T_c4, Fn_stations, Ft_stations, chords


def compute_internal_loads(y_stations, chords, T_c4, Fn_aero, Ft_aero, m_prime, x_sc_stations,
    n, y_eng, Fn_eng, Ft_eng, dz_eng, dx_eng):
    # 1. Calculate span (dy) for mass to point load conversion
    dy = np.zeros_like(y_stations)
    if len(y_stations) > 1:
        dy[1:-1] = (y_stations[2:] - y_stations[:-2]) / 2
        dy[0] = (y_stations[1] - y_stations[0]) / 2
        dy[-1] = (y_stations[-1] - y_stations[-2]) / 2
    Fn_weight = m_prime * 9.81 * dy

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