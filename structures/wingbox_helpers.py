import numpy as np
from scipy.integrate import cumulative_simpson


def calculate_wing_geometry(S, AR, taper):
    b = np.sqrt(S * AR)  # Total span (m)
    y_tip = b / 2  # Semi-span (m)
    c_root = (2 * S) / (b * (1 + taper))
    c_tip = taper * c_root
    C_mac = (2 / 3) * c_root * ((1 + taper + taper ** 2) / (1 + taper))
    return b, y_tip, c_root, c_tip, C_mac


def calculate_torsional_deflection(y, T, G, J, y_fus=0.0):
    # Slice arrays to only include stations from the fuselage wall outward
    struct_mask = y >= y_fus
    y_s, T_s, J_s = y[struct_mask], T[struct_mask], J[struct_mask]

    twist_rate = T_s / (G * J_s)

    # Cumulative integration naturally initializes at 0 for the first element (y_fus)
    twist_rad = cumulative_simpson(twist_rate, x=y_s, initial=0)

    # Reconstruct full array so it matches the original shape expected by output scripts
    full_twist_deg = np.zeros_like(y)
    full_twist_deg[struct_mask] = np.rad2deg(twist_rad)
    return full_twist_deg


def calculate_bending_deflection(y, M, E, I, y_fus=0.0):
    # Slice arrays to only include stations from the fuselage wall outward
    struct_mask = y >= y_fus
    y_s, M_s, I_s = y[struct_mask], M[struct_mask], I[struct_mask]

    curvature = M_s / (E * I_s)

    # First integration: Slope (Starts at 0 at y_fus)
    slope = cumulative_simpson(curvature, x=y_s, initial=0)

    # Second integration: Deflection (Starts at 0 at y_fus)
    deflection_m = cumulative_simpson(slope, x=y_s, initial=0)

    # Reconstruct full array so inboard is 0 and outboard is elastically deformed
    full_deflection_mm = np.zeros_like(y)
    full_deflection_mm[struct_mask] = deflection_m * 1000
    return full_deflection_mm



def load_aerodynamic_data(filepath, verbose=False):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    v_inf = 72.0  # Fallback default

    # Initialize variables for XFLR5 global parameters
    header_CL = 0.0
    header_Cd = 0.0
    header_Bending = 0.0

    # Parse header parameters
    for line in lines:
        line_stripped = line.strip()
        if line_stripped.startswith("QInf"):
            parts = line_stripped.replace('=', ' ').split()
            if len(parts) >= 2:
                v_inf = float(parts[1])
        elif line_stripped.startswith("CL"):
            parts = line_stripped.replace('=', ' ').split()
            if len(parts) >= 2:
                header_CL = float(parts[1])
        elif line_stripped.startswith("Cd"):
            parts = line_stripped.replace('=', ' ').split()
            if len(parts) >= 2:
                header_Cd = float(parts[1])
        elif line_stripped.startswith("Bending"):
            parts = line_stripped.replace('=', ' ').split()
            if len(parts) >= 2:
                header_Bending = float(parts[1])
        elif "y-span" in line and "Chord" in line and "Cl" in line:
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
            # Keep both sides for an accurate total area, lift, and drag calculation
            # if y < 0: continue  <- Removed this so we integrate the whole wing
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
        dy[0] = (y_stations[1] - y_stations[0]) / 2  # Root edge
        dy[-1] = (y_stations[-1] - y_stations[-2]) / 2  # Tip edge

    # Calculate final point loads per strip, small angle approx
    Fn_stations = -(q_inf * chords * cls) * dy  # Normal force (Positive Down)
    Ft_stations = (q_inf * chords * cds) * dy  # Tangential force (Positive Aft)
    T_c4 = (q_inf * (chords ** 2) * cms) * dy  # Pitching Moment - Torsion

    # UNIT TEST: Calculate Totals and Compare to XFLR5

    # Calculate Total Wing Area (S) from strips
    S_ref = np.sum(chords * dy)

    # Calculated Total Forces
    # Note: Fn_stations is negative for upward lift, so we negate it for total lift sum
    total_lift_calc = np.sum(-Fn_stations)
    total_drag_calc = np.sum(Ft_stations)

    # Calculate Root Bending Moment (Integral of L' * y dy over the half-wing)
    mask_half = y_stations >= 0
    root_bm_calc = np.sum(-Fn_stations[mask_half] * y_stations[mask_half])

    # Expected Totals from XFLR5 header
    total_lift_xflr = q_inf * S_ref * header_CL
    total_drag_xflr = q_inf * S_ref * header_Cd
    root_bm_xflr = header_Bending

    if verbose==True:
        # Output Comparison
        print(f"Wing Area (S_ref) : {S_ref:.4f} m^2")
        print("-" * 50)
        print(f"Total Lift        | Calc: {total_lift_calc:9.2f} N | XFLR5: {total_lift_xflr:9.2f} N")
        print(f"Total Drag        | Calc: {total_drag_calc:9.2f} N | XFLR5: {total_drag_xflr:9.2f} N")
        print(f"Root Bending Mmt  | Calc: {root_bm_calc:9.2f} Nm| XFLR5: {root_bm_xflr:9.2f} Nm")
        print("=" * 50 + "\n")

    # If your downstream script ONLY expects the right side of the wing,
    # you can slice the arrays here before returning:
    mask_right = y_stations >= 0

    return (
        y_stations[mask_right],
        T_c4[mask_right],
        Fn_stations[mask_right],
        Ft_stations[mask_right],
        chords[mask_right]
    )

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

def check_buckle(h_web, t_web, V_stations, E_spar, b_spacing):
    nu = 0.33  # Poisson's ratio

    # Shear buckling coefficient for web panel
    ratio = np.minimum(h_web, b_spacing) / np.maximum(h_web, b_spacing)
    ks = 5.35 + 4 * (ratio ** 2)

    # Actual shear stress per panel assumming about the same spar height
    A_web_panel = 2 * t_web * h_web
    tau_actual = np.abs(V_stations) / A_web_panel


    # Critical shear stress for buckling along web
    tau_critical = ks * ((np.pi ** 2 * E_spar) / (12 * (1 - nu ** 2))) * (t_web / h_web) ** 2

    # Margin of Safety
    with np.errstate(divide='ignore', invalid='ignore'):
        spar_margins = (tau_critical / tau_actual) - 1
        spar_margins = np.nan_to_num(spar_margins, nan=10.0, posinf=10.0)

    return float(np.min(spar_margins))