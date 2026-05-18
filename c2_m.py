from classes.aircraft_2 import Aircraft
from lookups.consts import *
import pandas as pd
import numpy as np
from scipy.interpolate import RegularGridInterpolator, interp1d
from scipy.optimize import brentq
from classes.isa import Atmosphere
import matplotlib.pyplot as plt

# NOTE: Check how to calculate W_ops for medivac, check fti mass, check if ballasts needed


# Helper functions for wing geometry:
def interpolate_csv_z(x_query, y_query, filepath):
    """
    Interpolate z value from a CSV lookup table.

    Parameters
    ----------
    x_query : float
        Desired x value.
    y_query : float
        Desired y value.
    filepath : str
        Path to CSV file.

    Returns
    -------
    float
        Interpolated z value.
    """

    # Read CSV
    df = pd.read_csv(filepath, index_col=0)

    # Extract x, y, z
    x_vals = df.columns.astype(float).to_numpy()
    y_vals = df.index.astype(float).to_numpy()
    z_vals = df.to_numpy(dtype=float)

    # Clamp to boundaries
    x_query = np.clip(x_query, x_vals.min(), x_vals.max())
    y_query = np.clip(y_query, y_vals.min(), y_vals.max())

    # Create interpolator
    interpolator = RegularGridInterpolator(
        (y_vals, x_vals),  # order matters: rows=y, cols=x
        z_vals,
        method='linear',
        bounds_error=False
    )

    # Query point
    z_interp = interpolator([[y_query, x_query]])[0]

    return float(z_interp)

def scissor_plot_intersection_points(x1, x2, y1, y2):
    # Interpolation functions
    # -----------------------------
    f1 = interp1d(x1, y1, kind='linear')
    f2 = interp1d(x2, y2, kind='linear')

    # Function for difference
    def diff(x):
        return f1(x) - f2(x)

    # -----------------------------
    # Find sign changes
    # -----------------------------
    x_search = np.linspace(
        max(min(x1), min(x2)),
        min(max(x1), max(x2)),
        1000
    )

    d = diff(x_search)

    # Indices where sign changes occur
    crossing_indices = np.where(np.diff(np.sign(d)))[0]

    intersections = []

    # -----------------------------
    # Refine intersections
    # -----------------------------
    for idx in crossing_indices:
        x_left = x_search[idx]
        x_right = x_search[idx + 1]

        # Exact root via interpolation
        x_root = brentq(diff, x_left, x_right)

        # Corresponding y-value
        y_root = f1(x_root)

        intersections.append((x_root, y_root))

    # -----------------------------
    # Choose lowest y-value
    # -----------------------------
    if intersections:
        x_best, y_best = min(intersections, key=lambda p: p[1])

        print(f"Selected intersection:")
        print(f"x = {x_best:.6f}")
        print(f"y = {y_best:.6f}")
        return x_best, y_best
    else:
        print("No intersections found")
        return None

def closest_value(x, values = [2, 4, 6, 8, 10]):
    return int(min(values, key=lambda v: abs(v - x)))

def LE_sweep_deg(sweep_c_4: float, # degrees
             c_r: float, #root chord length [m]
             b: float, # span [m]
             taper_ratio: float) -> float:
    return np.rad2deg(np.arctan(np.tan(np.deg2rad(sweep_c_4)) + 0.5 * c_r / b * (1 - taper_ratio)))

def sweep_at_x_c_deg(LE_sweep: float, # degrees
             c_r: float, #root chord length [m]
             b: float, # span [m]
             taper_ratio: float,
             x_c: float) -> float:
    return np.rad2deg(np.arctan(np.tan(np.deg2rad(LE_sweep)) - x_c * 2 * c_r / b * (1 - taper_ratio)))

def mac(c_r: float, # root chord [m]
        taper_ratio) -> float:
    return 2 * c_r / 3 * (1 + taper_ratio + taper_ratio**2) / (1 + taper_ratio)

def chord_at_y_span(c_r, taper_ratio, y, b):
    return c_r - 2 * y / b * c_r * (1 - taper_ratio)

def y_pos_at_chord_length(c_r, taper_ratio, chord, b):
    return (c_r - chord) / (2 / b * c_r * (1 - taper_ratio))

def x_pos_le_along_span_from_nose(le_sweep_deg: float, # Leading edge sweep in degrees
                                  y: float, # Location along semi-span
                                  x_le # Distance of leading edge of wing/tail surface from nose of aircraft
                                  ):
    return x_le + y * np.tan(np.deg2rad(le_sweep_deg))

def Snet(S, b_f, taper, b, c_r):
    c_fus_int = chord_at_y_span(c_r, taper, b_f/2, b)
    return S - (c_r + c_fus_int) * b_f / 2

def S_wf(y_start_f, y_end_f, taper, c_r, b):
    c_start_f = chord_at_y_span(c_r, taper, y_start_f, b)
    c_end_f = chord_at_y_span(c_r, taper, y_end_f, b)
    return (c_end_f + c_start_f) * np.abs(y_end_f - y_start_f)

def beta(V, altitude, temp_shift):
    M = V / np.sqrt(287 * 1.4 * Atmosphere(altitude, temp_shift)).temp
    return np.sqrt(1 - M**2)

def lift_slope(A, beta, sweep_c_2_deg, eta=0.95):
    return 2 * np.pi * A / (2 + np.sqrt(4 + (A * beta / eta)**2 * (1 + np.tan(np.deg2rad(sweep_c_2_deg))**2 / beta**2)))

def W_to(ac: Aircraft, w_oe, w_f, w_pl, w_crew = 0, update_ac: bool = False):
    W_to = sum(w_oe, w_f, w_pl, w_crew) * (1 + 0.005 / 0.995)  # 0.5% trapped fuel and oil
    if update_ac:
        ac.weights.m_takeoff = W_to
    return W_to

def W_oe_and_cg_from_nose(ac: Aircraft, x_le_w, x_le_ht, x_le_vt, update_ac: bool = False, 
                          pie_chart_output_path: str = None, show_pie_chart: bool = False, 
                          struc_pie_chart_output_path: str = None, struc_show_pie_chart: bool = False) -> tuple:
    w_power = W_pwr(ac)
    w_mlg, w_nlg = W_gear(ac)
    w_ht, w_vt = W_emp(ac)
    wwing = W_wing(ac)
    wfus = W_fus(ac)
    wnac = W_nac(ac)
    w_structure = sum(wwing, w_ht, w_vt, wfus, wnac, w_mlg, w_nlg)
    w_fxeq, x_cg_fxeq = W_feq_and_cg_from_nose(ac)
    W_oe = w_structure + w_power + w_fxeq
    x_cg_oe = (w_structure * x_cg_structural_from_nose(ac, x_le_w, x_le_vt, x_le_ht, update_ac=False)[0] + w_fxeq * x_cg_fxeq + w_power * x_cg_pwr_from_nose(ac))
    if update_ac:
        ac.weights.m_empty = W_oe
        ac.weights.x_cg_oew = x_cg_oe
    if pie_chart_output_path is not None:
        categories = ['Structural', 'Power', 'Fixed equipment']
        values = [w_structure, w_power, w_fxeq] / W_oe * 100
        fig, ax = plt.subplots(figsize=(8, 6))
        wedges, _ = ax.pie(values, startangle=90)
        total = sum(values)
        for i, wedge in enumerate(wedges):
            angle = (wedge.theta2 + wedge.theta1) / 2
            x = np.cos(np.radians(angle))
            y = np.sin(np.radians(angle))
            percentage = values[i] / total * 100

            ax.annotate(
                f'{categories[i]}\n{percentage:.1f}%',
                xy=(x, y),
                xytext=(1.3 * np.sign(x), 1.3 * y),
                arrowprops=dict(arrowstyle='->'),
                ha='left' if x > 0 else 'right'
            )

        ax.axis('equal')
        plt.title('Distribution of OEW')
        plt.savefig(pie_chart_output_path)
        if show_pie_chart:
            plt.show()
    if struc_pie_chart_output_path is not None:
        categories = ['Wing', 'Horizontal tail', 'Vertical tail', 'Fuselage', 'Nacelles', 'Landing gear']
        values = [wwing, w_ht, w_vt, wfus, wnac, w_mlg+w_nlg] / w_structure * 100
        fig, ax = plt.subplots(figsize=(8, 6))
        wedges, _ = ax.pie(values, startangle=90)
        total = sum(values)
        for i, wedge in enumerate(wedges):
            angle = (wedge.theta2 + wedge.theta1) / 2
            x = np.cos(np.radians(angle))
            y = np.sin(np.radians(angle))
            percentage = values[i] / total * 100

            ax.annotate(
                f'{categories[i]}\n{percentage:.1f}%',
                xy=(x, y),
                xytext=(1.3 * np.sign(x), 1.3 * y),
                arrowprops=dict(arrowstyle='->'),
                ha='left' if x > 0 else 'right'
            )

        ax.axis('equal')
        plt.title('Distribution of structural weight')
        plt.savefig(struc_pie_chart_output_path)
        if struc_show_pie_chart:
            plt.show()
    return W_oe, x_cg_oe, ac

def W_to_new(ac: Aircraft, x_le_w, x_le_ht, x_le_vt,
             m_ff, # from class I
             m_res, # from class I
             m_tfo, # from class I
             W_crew = 0.0, # included in PL probably
             update_ac: bool = False,
             pie_chart_output_path: str = None,
             show_pie_chart: bool = False
             ):
    W_PL = ac.weights.m_payload
    W_e = W_oe_and_cg_from_nose(ac, x_le_w, x_le_ht, x_le_vt)[0]
    m_ff = ... # Insert class I method called
    m_res = ... # Assume value NOTE: check if Shubhankar used the whole range + diversion for fuel mass est, else add here
    m_tfo = ... 
    W_to = (W_e + W_PL + W_crew) / (m_ff * (1 + m_res) - m_res - m_tfo)
    W_F = W_to-W_e-m_tfo-W_PL
    if update_ac:
        ac.weights.m_takeoff = W_to
        ac.weights.m_fuel = W_F

    if pie_chart_output_path is not None:
        categories = ['Empty weight', 'Payload', 'Trapped fuel and oil', 'Fuel']
        values = [W_e, W_PL, m_tfo, W_to-W_e-m_tfo-W_PL] / W_to * 100
        fig, ax = plt.subplots(figsize=(8, 6))
        wedges, _ = ax.pie(values, startangle=90)
        total = sum(values)
        for i, wedge in enumerate(wedges):
            angle = (wedge.theta2 + wedge.theta1) / 2
            x = np.cos(np.radians(angle))
            y = np.sin(np.radians(angle))
            percentage = values[i] / total * 100

            ax.annotate(
                f'{categories[i]}\n{percentage:.1f}%',
                xy=(x, y),
                xytext=(1.3 * np.sign(x), 1.3 * y),
                arrowprops=dict(arrowstyle='->'),
                ha='left' if x > 0 else 'right'
            )

        ax.axis('equal')
        plt.title('Distribution of MTOW')
        plt.savefig(pie_chart_output_path)
        if show_pie_chart:
            plt.show()
    return W_to, W_F, ac

def W_pwr(ac: Aircraft):
    P_to = ac.engine.power_to / HP_TO_W
    W_eng = ... / LBS_TO_KG  # Weight of one engine in lbs (from manufacturer)
    N_e = ac.engine.count
    alpha_p_id = ac.engine.alpha_p_id
    W_fuel = ac.weights.m_fuel
    Fuel_type = ...  # Check if avgas or Jet-A1 or other
    N_t = ... # Number of separate fuel tanks
    INT = ...  # Fraction of fuel tanks that are intergral
    # W_pwr = W_eng + W_ai + W_prop + W_fs + W_p + W_batt = W_batt + W_pwr1 + W_fs (fuel system)
    
    # USAF - Roskam eqn 6.3
    W_pwr1 = 2.575 * W_eng**0.922 * N_e

    if alpha_p_id != 'hydrogen':
        if Fuel_type == 'avgas':
            K_fsp = 5.87  # lbs/gal
        elif Fuel_type == 'Jet-A1':
            K_fsp = 6.55  # lbs/gal
        else:
            raise ValueError(f"Fuel type given: {Fuel_type}, should be either 'avgas' or 'Jet-A1'")
        W_fs = 2.49 * ((W_fuel / K_fsp)**0.6 * (1 / (1 + INT))**0.3 * N_t**0.2 * N_e**0.13)**1.21
    
    elif alpha_p_id == 'hydrogen':  # NOTE: insert method of fixed equipment weight
        W_fs = ...
    return W_pwr1 + W_fs + ac.weights.m_battery

def W_feq_and_cg_from_nose(ac: Aircraft):
    W_to = ac.weights.m_takeoff / LBS_TO_KG
    W_e = ac.weights.m_empty
    N_pax = ...  # Including crew
    M_D = ...  # Design dive Mach number
    no_pressurization_const = ...  # Fraction to take into account api does not have p (pressurization) in our case
    V_pax_cargo = ... * M2_TO_F2/FT_TO_M  # Volume of passenger cabin and cargo [ft3]

    W_fc = 0.5 * (0.0168 * W_to + 1.066 * W_to**0.626)

    # Cessna
    W_hps = 0.009 * W_to
    W_els = 0.0268 * W_to
    # # Torenbeek (W_hps + W_els = W_hps_els)
    # W_hps_els = 0.0078 * W_to**1.2
    W_iae = 40 + 0.008 * W_to
    W_api = (0.265 * W_to**0.52 * N_pax**0.68 * W_iae**0.17 * M_D**0.08) * no_pressurization_const
    W_fur = 0.5 * (0.412 * N_pax**1.145 * W_to**0.489 + 15 * N_pax + V_pax_cargo)
    W_ops = 0
    W_fti = 0.5 * (155 / 9980 * W_to + 708 / 24912 * W_to)
    W_aux = 0.01 * W_e
    W_bal = 0
    W_pt = 0.0045 * W_to
    W_etc = 0  # NOTE: check if there are other things to include
    Weights = np.array([W_fc, W_hps, W_els, W_iae, W_api, W_fur, W_ops, W_fti, W_aux, W_pt, W_etc])
    W_feq = np.sum(Weights)

    # CG from nose - *indicated source = https://archive.aoe.vt.edu/mason/Mason_f/M96SC02.pdf 
    c_r_w = ...
    taper_w = ...
    b_w = ...
    b_ht = ...
    b_vt = ...
    w_sweep_c_4_deg = ...
    x_le_w = ...  # Wing root chord LE distance from nose
    x_le_ht = ...  # HT root chord LE distance from nose
    x_le_vt = ...  # VT root chord LE distance from nose
    x_c_front_spar = ...
    x_c_rear_spar = ...
    start_flap_along_span = ...
    end_aileron_along_span = ...
    l_fus = ac.fuselage.length
    l_tc = ...  # length of tailcone
    l_cabin = ...
    start_cabin = ...
    x_nlg = ...

    chord_fc = (chord_at_y_span(c_r_w, taper_w, start_flap_along_span, b_w) + chord_at_y_span(c_r_w, taper_w, end_aileron_along_span, b_w)) / 2
    x_cg_fc = (x_c_rear_spar + 1) / 2 * chord_fc + x_pos_le_along_span_from_nose(LE_sweep_deg(w_sweep_c_4_deg, c_r_w, b_w, taper_w), (start_flap_along_span+end_aileron_along_span)/2, x_le_w)  # *Between aft spar and trailing-edge
    x_cg_hps = ((x_c_rear_spar - x_c_front_spar) / 2 * c_r_w + x_le_w) * 0.7 + (l_fus - 0.5 * l_tc) * 0.3  # *
    x_cg_els = x_cg_hps * 0.7  # * Battery cables between 
    x_cg_iae = x_nlg * 0.65 # *
    x_cg_api = (b_w * x_le_w + b_ht * x_le_ht + b_vt * x_le_vt) / (b_w + b_ht + b_vt)  # De-icing on leading edges of wing, ht and vt
    x_cg_fur = start_cabin + 0.525 * l_cabin  # *
    x_cg_ops = x_cg_fur  # *
    x_cg_fti = start_cabin + l_cabin  # Source: we saw it was in the back for the PH-LAB
    x_cg_aux = start_cabin  # Assuming additional stuff like fire axes near cockpit
    x_cg_pt = start_cabin + 0.5 * l_cabin  # Guestimation
    x_cg_etc = 0  # NOTE: check if there was anything to add
    cgs = np.array([x_cg_fc, x_cg_hps, x_cg_els, x_cg_iae, x_cg_api, x_cg_fur, x_cg_ops, x_cg_fti, x_cg_aux, x_cg_pt, x_cg_etc])
    x_cg_feq = (Weights@cgs) / W_feq

    return W_feq, x_cg_feq

def W_wing(ac: Aircraft, update_ac: bool = False):

    W_to = ac.weights.m_takeoff / LBS_TO_KG
    S = ac.wing.area * M2_TO_F2
    n_ult = ...  # NOTE: add later
    A = ac.wing.aspect_ratio
    t_c_max = ...
    sweep_c_4_deg = ...
    taper = ac.wing.taper_ratio
    V_h = ...  # NOTE: add later: maximum level speed at sealevel in kts
    wing_type = ...  # 'cantilever' or not
    b = ac.wing.span / FT_TO_M
    c_r = ... / FT_TO_M
    t_r = t_c_max * c_r  # max thickness of wing root chord in ft
    sweep_le_rad = np.arctan(np.tan(np.rad2deg(sweep_c_4_deg)) + c_r / 2 / b * (1 + taper))
    sweep_c_2_rad = np.arctan(np.tan(sweep_le_rad) - c_r / b * (1 + taper))

    # Cessna
    if wing_type == 'cantilever':
        W_wing_c = 0.04674 * (W_to ** 0.397) * (S ** 0.36) * (n_ult ** 0.397) * (A ** 1.712)

    else:
        W_wing_c = 0.002933 * (S ** 1.018) * (A ** 2.473) * (n_ult ** 0.611)

    # USAF
    W_wing_u = 96.948 * (((W_to * n_ult * 1e-5) ** 0.65) * ((A / np.cos(np.rad2deg(sweep_c_4_deg))) ** 0.57) * ((S / 100) ** 0.61) * (((1 + taper) / (2 * t_c_max)) ** 0.36) * ((1 + V_h / 500) ** 0.5)) ** 0.993

    # Torenbeek
    W_wing_t = 0.00125 * W_to * ((b / np.cos(sweep_c_2_rad)) ** 0.75) * (1 + (6.3 * np.cos(sweep_c_2_rad) / b) ** 0.5) * (n_ult ** 0.55) * (b * S / (t_r * W_to * np.cos(sweep_c_2_rad)))**0.3

    return (W_wing_c + W_wing_u + W_wing_t) / 3

def W_emp(ac: Aircraft, update_ac: bool = False):
    n_ult = ...  # NOTE: add later
    S_h = ... * M2_TO_F2
    S_v = ... * M2_TO_F2
    l_h = ... / FT_TO_M
    b_h = ... / FT_TO_M
    b_v = ... / FT_TO_M
    t_c_max_h = ...
    t_c_max_v = ...
    c_r_h = ... / FT_TO_M
    c_r_v = ... / FT_TO_M
    t_r_h = t_c_max_h * c_r_h
    t_r_v = t_c_max_v * c_r_v
    A_h = ...
    A_v = ...
    vt_sweep_c_4_deg = ...

    # Cessna
    W_h_c = (3.184 * (W_to ** 0.887) * (S_h ** 0.101) * (A_h ** 0.138)) / (174.04 * (t_r_h **0.223))
    W_v_c = (1.68 * (W_to ** 0.567) * (S_v ** 1.249) * (A_v ** 0.482)) / (639.95 * (t_r_v ** 0.747) * (np.cos(np.rad2deg(vt_sweep_c_4_deg)) ** 0.882))
    W_emp_c = W_h_c + W_v_c

    # USAF
    W_h_u = 127 * ((W_to * n_ult * 1e-5)**0.87 * (S_h / 100)**1.2 * 0.289 * (l_h / 10)**0.483 * (b_h / t_r_h)**0.5)**0.458
    W_v_u = 98.5 * ((W_to * n_ult * 1e-5)**0.87 * (S_v / 100)**1.2 * 0.289 * (b_v / t_r_v)**0.5)**0.458
    W_emp_u = W_h_u + W_v_u

    # # Torenbeek
    # W_emp_t = 0.04 * (n_ult * (S_v + S_h)**2)**0.75

    W_ht = (W_h_c + W_h_u) / 2
    W_vt = (W_v_c + W_v_u) / 2

    return W_ht, W_vt

def W_fus(ac: Aircraft, update_ac: bool = False):
    fus = ac.fuselage
    l_f = fus.length / FT_TO_M  # Fuselage length in ft
    w_f = fus.width/ FT_TO_M  # Max fuselage width in ft
    h_f = fus.height / FT_TO_M  # Max fuselage height in ft
    n_ult = ...
    alt_c = ac.requirements.cruise['cr_altitude'] * FT_TO_M  # m
    rho_c = Atmosphere(height=alt_c) # kg/m3
    V_c = ac.requirements.cruise['cr_speed'] * np.sqrt(rho_c / 1.225) # KEAS
    l_f_n = l_f # Fuselage length - nacelle in ft
    P_max = ...  # Maximum fuselage perimeter in ft
    N_pax = ...  # Number of passengers including pilot(s)

    # Cessna
    W_fus_c = 14.86 * (W_to**0.144) * (l_f_n / P_max)**0.778 * (l_f_n**0.383) * (N_pax**0.455)
    
    # USAF
    W_fus_u = 200 * ((W_to * n_ult * 1e-5)**0.286 * (l_f / 10)**0.857 * ((w_f + h_f) / 10) * (V_c / 100)**0.338)**1.1

    return (W_fus_c + W_fus_u) / 2

def W_nac(ac: Aircraft, update_ac: bool = False):
    alpha_p_id = ac.engine.alpha_p_id  # turboprop or hydrogen or piston
    engine_over_wing: bool = ...
    P_to = ac.engine.power_to / HP_TO_W  # HP

    # Cessna NOTE: fill in methods
    if alpha_p_id == 'piston':
        engine_type = ...
        if engine_type == 'radial':
            K_n = 0.37   # lbs/hp 
        elif engine_type == 'horizontally opposed':
            K_n = 0.24   # lbs/hp 

        W_n_c = K_n * P_to
    elif alpha_p_id == 'hydrogen':
        W_n_c = ...
    elif alpha_p_id =='turboprop':
        W_n_c = ...

    # Torenbeek
    if alpha_p_id == 'piston':
        engine_type = ...
        if engine_type == 'radial':
            N_e = ac.engine.count
            W_n_t = 0.045 * P_to**1.25 * N_e**(-0.25)
        elif engine_type == 'horizontally opposed':
            W_n_t = 0.32 * P_to 

    elif alpha_p_id == 'hydrogen':
        W_n_t = ...
    elif alpha_p_id =='turboprop':
        W_n_t = 0.14 * P_to 

    if engine_over_wing:
        W_n_t += 0.11 * P_to

    return (W_n_c + W_n_t) / 2

def W_gear(ac: Aircraft, update_ac: bool = False):
    W_to = ac.weights.m_takeoff / LBS_TO_KG
    W_L = ac.requirements.landing['la_mass_frac'] * W_to
    shock_strut_frac_whole = 
    l_s_m = shock_strut_frac_whole * (np.abs(ac.landing_gear.height_mlg) - ac.landing_gear.selected_mlg_tire['Tire Radius (In)'] * 2.54 / 100) / FT_TO_M # Shock strut length main gear [ft]
    l_s_n = shock_strut_frac_whole * (np.abs(ac.landing_gear.height_nlg) - ac.landing_gear.selected_nlg_tire['Tire Radius (In)'] * 2.54 / 100) / FT_TO_M # Shock strut length nose gear [ft]
    n_ult = ...

    # Cessna
    W_mlg = (0.013 * W_to + 0.362 * (W_L**0.417) * (n_ult**0.950) * (l_s_m**0.183))
    W_nlg = (6.2 + 0.0013 * W_to + 0.007157 * (W_L**0.749) * (n_ult) * (l_s_n**0.788))
    
    # # USAF
    # W_g_u = 0.054 * (l_s_m**0.501) * (W_L * n_ult)**0.684

    return W_mlg, W_nlg

def x_cg_structural_from_nose(ac: Aircraft,
                            x_le_w: float,
                            x_le_vt: float,
                            x_le_ht: float,
                            update_ac=False):
    l_fus = ac.fuselage.length
    w_sweep_c_4_deg = ac.wing.sweep
    ht_sweep_c_4_deg = ...
    vt_sweep_c_4_deg = ...
    w_c_r = ...
    ht_c_r = ...
    vt_c_r = ...
    w_taper = ac.wing.taper_ratio
    ht_taper = ...
    vt_taper = ...
    b_w = ac.wing.span
    b_ht = ...
    b_vt = ...
    x_c_front_spar = ...
    x_c_rear_spar = ...
    t_tail_condition: bool = ...  # is it a t-tail or not
    l_nacelle = ...
    nac_mount_dist = ...  # distance nacelle is mounted behind le of wing (front of nacelle)
    x_nlg = ...
    x_mlg = ...
    
    # Wing:
    if ac.wing.sweep == 0:
        c_y = chord_at_y_span(w_c_r, w_taper, y=0.4*0.5*b_w, b=b_w)
        x_along_chord = 0.4 * c_y
        x_cg_w = x_pos_le_along_span_from_nose(LE_sweep_deg(sweep_c_4=0, c_r=w_c_r, b=b_w, taper_ratio=w_taper), 0.5 * 0.4 * b_w, x_le_w) + x_along_chord
    else:
        c_y = chord_at_y_span(w_c_r, w_taper, y=0.35*0.5*b_w, b=b_w)
        x_along_chord = c_y * (x_c_front_spar + 0.7 * (x_c_rear_spar - x_c_front_spar))
        x_cg_w = x_pos_le_along_span_from_nose(LE_sweep_deg(w_sweep_c_4_deg, w_c_r, b_w, w_taper), 0.5 * 0.35 * b_w, x_le_w) + x_along_chord

    # HT:
    c_y = chord_at_y_span(ht_c_r, ht_taper, y=0.38*0.5*b_ht, b=b_ht)
    x_along_chord = 0.42 * c_y
    x_cg_ht = x_pos_le_along_span_from_nose(LE_sweep_deg(ht_sweep_c_4_deg, ht_c_r, b_ht, ht_taper), 0.38*0.5*b_ht, x_le_ht) + x_along_chord
    # VT:
    if t_tail_condition:
        c_y = chord_at_y_span(vt_c_r, vt_taper, y=0.55*0.5*b_vt, b=b_vt)
        x_along_chord = 0.42 * c_y
        x_cg_vt = x_pos_le_along_span_from_nose(LE_sweep_deg(vt_sweep_c_4_deg, vt_c_r, b_vt, vt_taper), 0.55*0.5*b_vt, x_le_vt) + x_along_chord
    else:
        c_y = chord_at_y_span(vt_c_r, vt_taper, y=0.38*0.5*b_vt, b=b_vt)
        x_along_chord = 0.42 * c_y
        x_cg_vt = x_pos_le_along_span_from_nose(LE_sweep_deg(vt_sweep_c_4_deg, vt_c_r, b_vt, vt_taper), 0.38*0.5*b_vt, x_le_vt) + x_along_chord
    # Fuselage:
    x_cg_fus = 0.45 * l_fus  # https://archive.aoe.vt.edu/mason/Mason_f/M96SC02.pdf

    # Nacelle:
    n_eng = ac.engine.count
    if n_eng % 2 !=0:
        raise ValueError(f"Number of engines = {n_eng}, need an even number")
    if n_eng == 2:
        y_pos_eng: float = ...
        x_cg_nac = 0.4 * l_nacelle + x_pos_le_along_span_from_nose(LE_sweep_deg(w_sweep_c_4_deg, w_c_r, b_w, w_taper), y_pos_eng, x_le_w) + nac_mount_dist
    elif n_eng > 2:
        y_pos_eng: list = ...  # NOTE: get a list/array or someting with engine y-positions for one side
        if len(y_pos_eng) != n_eng // 2:
            raise ValueError(
                f"Error: len(y_pos_eng) = {len(y_pos_eng)}, "
                f"but expected {n_eng // 2}"
            )
        x_cg_nac = 0
        for i, y_pos in enumerate(y_pos_eng):
            x_cg_nac_elem = 0.4 * l_nacelle + x_pos_le_along_span_from_nose(LE_sweep_deg(w_sweep_c_4_deg, w_c_r, b_w, w_taper), y_pos, x_le_w) + nac_mount_dist
            x_cg_nac += 2 * x_cg_nac_elem / n_eng

    # Landing gear:
    x_cg_nlg = x_nlg
    x_cg_mlg = x_mlg

    x_cg_data = {'wing': x_cg_w, 'ht': x_cg_ht, 'vt': x_cg_vt, 'fuselage': x_cg_fus, 
                 'nose landing gear': x_cg_nlg, 'main landing gear': x_cg_mlg, 'nacelle': x_cg_nac}
    W_ht, W_vt = W_emp(ac)
    W_mlg, W_nlg = W_gear(ac)
    Weights = np.array([W_wing(ac), W_ht, W_vt, W_fus(ac), W_mlg, W_nlg, W_nac(ac)])
    cgs = np.array([x_cg_w, x_cg_ht, x_cg_vt, x_cg_fus, x_cg_mlg, x_cg_nlg, x_cg_nac])
    x_cg_struc = (Weights @ cgs) / np.sum(Weights)
    if update_ac:
        ac.weights.x_cg_structural = x_cg_data
    return x_cg_struc, x_cg_data, ac

def x_cg_pwr_from_nose(ac: Aircraft):
    ...

def convert_x_cg_from_nose_to_lemac_frac_mac(x_cg_nose, ac: Aircraft):
    c_r = ...
    taper = ac.wing.taper_ratio
    sweep_c_4 = ac.wing.sweep
    x_le = ...
    b = ac.wing.span
    mac = 2 * c_r / 3 * (1 + taper + taper**2) / (1 + taper)
    le_sweep_deg = LE_sweep_deg(sweep_c_4, c_r, b, taper)
    y_span_mac = (c_r - mac) / (2 / b * c_r * (1 - taper))
    x_lemac = x_pos_le_along_span_from_nose(le_sweep_deg, y_span_mac, x_le)
    return (x_cg_nose - x_lemac) / mac

def update_m_and_cg(m_current, cg_current, added_m, added_cg):
    new_cg = (m_current * cg_current + added_m * added_cg) / (m_current + added_cg)
    new_m  =(m_current + added_cg)
    return new_m, new_cg

def loading_diagram(x_le_w, ac: Aircraft, show_plot: bool=False, output_filepath=None, update_ac_cgs = False):
    m_cargo = ac.weights.m_cargo
    x_le_ht = ...
    x_le_vt = ...
    x_cg_seats = ...  # Front to back cg positions of seats (length 3)
    assert len(x_cg_seats) == 3, f"There must be 3 seat cg positions, {len(x_cg_seats)} were given"
    x_cg_seats = convert_x_cg_from_nose_to_lemac_frac_mac(x_cg_seats, ac)
    n_pax = ...
    m_pax = ac.weights.m_pax
    n_rows = n_pax / len(x_cg_seats)
    n_window = ...  # should be 2
    n_aisle = ...
    n_middle = ...
    m_fuel = ac.weights.m_fuel
    x_cg_fuel_tanks = ...
    x_cg_fuel_tanks = convert_x_cg_from_nose_to_lemac_frac_mac(x_cg_fuel_tanks, ac)
    fuel_tanks_mass_fracs = ...
    n_fuel_tanks = len(x_cg_fuel_tanks)
    m_fuel_tanks = fuel_tanks_mass_fracs * m_fuel

    x_cg_cargo_holds = ...
    x_cg_cargo_holds = convert_x_cg_from_nose_to_lemac_frac_mac(x_cg_cargo_holds, ac)
    cargo_mass_frac = ...  # list same length as x_cg_cargo_holds with corresponding fractions of cargo mass held at specific cg. point


    # Initiate lists to plot:
    x_cg_lemac_plot_ftb = []
    W_plot_ftb = []
    x_cg_lemac_plot_btf = []
    W_plot_btf = []

    # OEW
    W_oe, x_cg_oe = W_oe_and_cg_from_nose(ac, x_le_w, x_le_ht, x_le_vt)
    W_plot_ftb.append(W_oe)
    x_cg_lemac_plot_ftb.append(convert_x_cg_from_nose_to_lemac_frac_mac(x_cg_oe, ac))
    W_plot_btf.append(W_oe)
    x_cg_lemac_plot_btf.append(convert_x_cg_from_nose_to_lemac_frac_mac(x_cg_oe, ac))

    # Cargo/PL
    cargo_masses = cargo_mass_frac * m_cargo
    cargo_masses_reverse = cargo_masses[::-1]
    x_cg_cargo_reverse = x_cg_cargo_holds[::-1]

    for i, x_cg in enumerate(x_cg_cargo_holds):
        new_m, new_cg = update_m_and_cg(W_plot_ftb[-1], x_cg_lemac_plot_ftb[-1], cargo_masses[i], x_cg)
        W_plot_ftb.append(new_m)
        x_cg_lemac_plot_ftb.append(new_cg)
        new_m, new_cg = update_m_and_cg(W_plot_ftb[-1], x_cg_lemac_plot_ftb[-1], cargo_masses_reverse[i], x_cg_cargo_reverse[i])
        W_plot_btf.append(new_m)
        x_cg_lemac_plot_btf.append(new_cg)
    assert W_plot_btf[-1] == W_plot_ftb[-1], f"After loading all cargo holds total weights should match. \n ftb weight = {W_plot_ftb[-1]} \t btf weight = {W_plot_btf[-1]}"
    assert x_cg_lemac_plot_btf[-1] == x_cg_lemac_plot_ftb[-1], f"After loading all cargo holds total cg should match. \n ftb cg = {x_cg_lemac_plot_ftb[-1]} \t btf cg = {x_cg_lemac_plot_btf[-1]}"
    
    # Passengers
    if n_aisle == 0 and n_middle == 0:
        pax_row_mass = m_pax / n_rows
        x_cg_seats_reverse = x_cg_seats[::-1]

        for i, x_cg in enumerate(x_cg_seats):
            new_m, new_cg = update_m_and_cg(W_plot_ftb[-1], x_cg_lemac_plot_ftb[-1], pax_row_mass, x_cg)
            W_plot_ftb.append(new_m)
            x_cg_lemac_plot_ftb.append(new_cg)
            new_m, new_cg = update_m_and_cg(W_plot_ftb[-1], x_cg_lemac_plot_ftb[-1], pax_row_mass, x_cg_seats_reverse[i])
            W_plot_btf.append(new_m)
            x_cg_lemac_plot_btf.append(new_cg)
        assert W_plot_btf[-1] == W_plot_ftb[-1], f"After loading all passengers total weights should match. \n ftb weight = {W_plot_ftb[-1]} \t btf weight = {W_plot_btf[-1]}"
        assert x_cg_lemac_plot_btf[-1] == x_cg_lemac_plot_ftb[-1], f"After loading all passengers total cg should match. \n ftb cg = {x_cg_lemac_plot_ftb[-1]} \t btf cg = {x_cg_lemac_plot_btf[-1]}"
    
    # Fuel
    m_fuel_tanks_reverse = m_fuel_tanks[::-1]
    x_cg_fuel_tanks_reverse = x_cg_fuel_tanks[::-1]
    for i, x_cg in enumerate(x_cg_fuel_tanks):
        new_m, new_cg = update_m_and_cg(W_plot_ftb[-1], x_cg_lemac_plot_ftb[-1], m_fuel_tanks[i], x_cg)
        W_plot_ftb.append(new_m)
        x_cg_lemac_plot_ftb.append(new_cg)
        new_m, new_cg = update_m_and_cg(W_plot_ftb[-1], x_cg_lemac_plot_ftb[-1], m_fuel_tanks_reverse[i], x_cg_fuel_tanks_reverse[i])
        W_plot_btf.append(new_m)
        x_cg_lemac_plot_btf.append(new_cg)
    assert W_plot_btf[-1] == W_plot_ftb[-1], f"After loading all fuel tanks total weights should match. \n ftb weight = {W_plot_ftb[-1]} \t btf weight = {W_plot_btf[-1]}"
    assert x_cg_lemac_plot_btf[-1] == x_cg_lemac_plot_ftb[-1], f"After loading all fuel tanks total cg should match. \n ftb cg = {x_cg_lemac_plot_ftb[-1]} \t btf cg = {x_cg_lemac_plot_btf[-1]}"
    
    # Plotting
    min_x = min(min(x_cg_lemac_plot_btf), min(x_cg_lemac_plot_ftb))
    max_x = max(max(x_cg_lemac_plot_btf), max(x_cg_lemac_plot_ftb))
    range_x = np.abs(max_x - min_x)
    min_y = min(min(W_plot_btf), min(W_plot_ftb))
    max_y = max(max(W_plot_btf), max(W_plot_ftb))
    range_y = np.abs(max_y - min_y)

    if output_filepath is not None:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.set_xlim(min_x - 0.1 * range_x, max_x + 0.1 * range_x)
        ax.set_ylim(min_y, max_y + 0.1 * range_y)
        ax.scatter(x_cg_lemac_plot_btf, W_plot_btf, color='green', label='loading back to front', marker = 'x', linewidth=2)
        ax.scatter(x_cg_lemac_plot_ftb, W_plot_ftb, color='purple', label='loading front to back', marker = 'x', linewidth=2)
        for i, val in enumerate((min_x, max_x)):
            ax.vline(val, color='grey', linestyle='--')
        for i, val in enumerate((min_x-0.05*range_x, max_x+0.05*range_x)):
            ax.vline(val, color='black', label='cg range (including 5% margin)')

        # Labels and shit
        ax.set_xlabel(f"x$_{{cg}}$/mac from lemac")
        ax.set_ylabel(f"Weight [kg]")
        ax.set_title("Loading diagram")
        ax.legend()
        ax.grid(True, linestyle="--", alpha=0.5)

        plt.tight_layout()
        plt.savefig(output_filepath, dpi=300)
        if show_plot:
            plt.show()
        print(f'Loading diagram saved to {output_filepath}')
    
    fwd_cg = min_x - 0.05 * range_x
    aft_cg = max_x + 0.05 * range_x
    if update_ac_cgs:
        ac.weights.x_cg_fwd = fwd_cg
        ac.weights.x_cg_aft = aft_cg
    return fwd_cg, aft_cg, ac

def scissor_plot(ac: Aircraft, x_cg_lemac_mac: np.ndarray, SM: float = 0.05, output_filepath=None, show_plot: bool = False)->np.ndarray:
    n_eng = ac.engine.count
    t_tail_condition: bool = ...  # T-tail or not
    tail_type: str =  # 'fully moving' or 'adjustable' or 'fixed'

    # Flight conditions
    V_cruise = ac.requirements.cruise['cr_speed'] * KTS_TO_MS
    alt_cr = ac.requirements.cruise['cr_altitude'] * FT_TO_M
    V_app = ac.requirements.approach['ap_speed'] * KTS_TO_MS
    alt_LD = ac.requirements.cruise['la_altitude'] * FT_TO_M
    LD_temp_shift = ac.requirements.landing['la_temp_shift']

    # Fuselage params
    b_f = ac.fuselage.width
    l_f = ac.fuselage.length
    h_f = ac.fuselage.height

    # Wing params
    x_le =  # wing le distance from nose
    A = ac.wing.aspect_ratio
    b = ac.wing.span
    S = ac.wing.area
    c_r = 
    taper = ac.wing.taper_ratio
    S_net = Snet(S, b_f, taper, b, c_r)
    sweep_c_4_deg = ac.wing.sweep
    sweep_c_2_deg = sweep_at_x_c_deg(LE_sweep_deg(sweep_c_4_deg, c_r, b, taper), c_r, b, taper, x_c=0.5)
    sweep_LE_deg = LE_sweep_deg(sweep_c_4_deg, c_r, b, taper)
    C_m0_airfoil = 
    C_L_0 =  # CL of flapped wing at alpha = 0
    Delta_Cl_max =  # the airfoil lift coefficient increase due to flap extension at landing condition (estimated in the wing design module)
    C_L_LD = 
    mac = 2 * c_r / 3 * (1 + taper + taper**2) / (1 + taper)
    mgc = S/b

    # HLD
    flap_type: str =  # containing fowler or slotter or none for Cm calculations
    flap_defl_ld_deg = 
    ext_flap_chord_ratio =  # c’/c = the ratio between the chord of the airfoil with extended flap and the chord in clean configuration
    cf_c =  # simple flaps: 0.25, highly efficient slotted flaps: 0.35-0.4
    cf_ext_flap_chord = cf_c / ext_flap_chord_ratio  # c_f/c'
    y_start_f = 
    y_end_f = 
    flap_span_wing_span = np.abs(y_start_f - y_end_f) / b
    Swf = S_wf(y_start_f, y_end_f, taper, c_r, b)

    # Tail params
    height_ht =  # height of ht above wing (perpendicular to chord lines)
    l_h = 
    A_h = 
    r = 2 * l_h / b
    m_tv = height_ht * 2 / b

    # Nacelle params
    nac_mount_dist =  # Distance front of nacelle from nose
    l_fn = x_pos_le_along_span_from_nose(sweep_LE_deg, b_f/2, x_le)
    b_n =  # nacelle diameter (max)
    x_mac_c_4 = x_pos_le_along_span_from_nose(sweep_LE_deg, y_pos_at_chord_length(c_r, taper, mac, b), x_le) + 0.25 * mac
    l_n = x_mac_c_4 - nac_mount_dist
    if any(word in flap_type.lower() for word in ['fowler', 'slotted']):
        fowler_data = np.array([[0.15, 0.2, 0.25, 0.3, 0.35, 0.4],
                                [0.235, 0.222, 0.213, 0.207, 0.202, 0.2]])
        cf_ext_flap_chord = np.clip(cf_ext_flap_chord, fowler_data[0,:].min(), fowler_data[0,:].max())
        mu1 = float(np.interp(cf_ext_flap_chord, fowler_data[0,:], fowler_data[1,:]))
    else: 
        mu1 = interpolate_csv_z(flap_defl_ld_deg, cf_ext_flap_chord, "lookups/mu_1.csv")
    mu2 = interpolate_csv_z(taper, flap_span_wing_span, 'lookups/mu_2.csv')
    mu3 = interpolate_csv_z(taper, flap_span_wing_span, 'lookups/mu_3.csv')
    
    # Lift coefficients and Lift curve slopes
    if tail_type == 'fully moving':
        C_L_h = -1.0
    elif tail_type == 'adjustable':
        C_L_h = -0.8
    elif tail_type == 'fixed':
        C_L_h = -0.35 * A_h**(1 / 3)
    C_L_A_less_h_ld = C_L_LD
    C_L_alpha_h_cr = lift_slope(A_h, beta(V_cruise, alt_cr, temp_shift=0), sweep_c_2_deg)  # /rad
    C_L_alpha_w_cr = lift_slope(A, beta(V_cruise, alt_cr, temp_shift=0), sweep_c_2_deg)  # /rad
    C_L_alpha_A_less_h_cr = C_L_alpha_w_cr * (1 + 2.15 * b_f / b) * S_net / S + np.pi / 2 * b_f**2 / S  # /rad
    # C_L_alpha_h_LD = lift_slope(A_h, beta(V_app, alt_LD, temp_shift=LD_temp_shift), sweep_c_2_deg)  # /rad
    C_L_alpha_w_LD = lift_slope(A, beta(V_app, alt_LD, temp_shift=LD_temp_shift), sweep_c_2_deg)  # /rad
    C_L_alpha_A_less_h_LD = C_L_alpha_w_LD * (1 + 2.15 * b_f / b) * S_net / S + np.pi / 2 * b_f**2 / S  # /rad

    # Aerodynamic center
    beta_A_cr = A * beta(V_cruise, alt_cr, temp_shift=0)
    beta_A_ld = A * beta(V_app, alt_LD, LD_temp_shift)
    sweep_beta_cr = np.rad2deg(np.arctan(np.tan(np.deg2rad(sweep_c_4_deg)) / beta_A_cr * A))
    sweep_beta_ld = np.rad2deg(np.arctan(np.tan(np.deg2rad(sweep_c_4_deg)) / beta_A_ld * A))
    data_val_cr = closest_value(beta_A_cr)
    data_val_ld = closest_value(beta_A_ld)
    x_ac_w_cr = interpolate_csv_z(taper, sweep_beta_cr, f'lookups/Lambda_beta_ABeta{data_val_cr}.csv')
    x_ac_w_ld = interpolate_csv_z(taper, sweep_beta_ld, f'lookups/Lambda_beta_ABeta{data_val_ld}.csv')
    x_ac_f1_cr = -1.8 * b_f * h_f * l_fn / (C_L_alpha_A_less_h_cr)
    x_ac_f1_ld = -1.8 * b_f * h_f * l_fn / (C_L_alpha_A_less_h_LD)
    x_ac_f2 = 0.273 * b_f * mgc * (b - b_f) * np.tan(np.deg2rad(sweep_c_4_deg)) / ((1 + taper) * mac**2 * (b + 2.15 * b_f))
    x_ac_n_cr = n_eng * (-4) * b_n**2 * l_n / (S * mac * C_L_alpha_A_less_h_cr)
    x_ac_n_ld = n_eng * (-4) * b_n**2 * l_n / (S * mac * C_L_alpha_A_less_h_LD)
    x_ac_cr = x_ac_w_cr + x_ac_f1_cr + x_ac_f2 + x_ac_n_cr
    x_ac_ld = x_ac_w_ld + x_ac_f1_ld + x_ac_f2 + x_ac_n_ld

    # depsilon/dalpha
    K_eL = (0.1124 + 0.1265 * np.deg2rad(sweep_c_4_deg) + 0.1766 * np.deg2rad(sweep_c_4_deg)**2) / (r**2) + 0.1024 / r + 2
    K_eL0 = 0.1124 / (r**2) + 0.1024 / r + 2
    depsilon_dalpha_cr = K_eL / K_eL0 * (r / (r**2 + m_tv**2) * 0.4876 / np.sqrt(r**2 + 0.6319 + m_tv**2) + (1 + (r**2 / (r**2 + 0.7915 + 5.0734 * m_tv**2))**0.3113) * (1 - np.sqrt(m_tv**2 / (m_tv**2 + 1)))) * C_L_alpha_w_cr / (np.pi * A)
    depsilon_dalpha_ld = K_eL / K_eL0 * (r / (r**2 + m_tv**2) * 0.4876 / np.sqrt(r**2 + 0.6319 + m_tv**2) + (1 + (r**2 / (r**2 + 0.7915 + 5.0734 * m_tv**2))**0.3113) * (1 - np.sqrt(m_tv**2 / (m_tv**2 + 1)))) * C_L_alpha_w_LD / (np.pi * A)
    V_h_V2 = 0.85
    if t_tail_condition:
        V_h_V2 = 1.0

    # Cm_ac
    C_m_ac_w = C_m0_airfoil * (A * np.cos(np.deg2rad(sweep_c_4_deg))**2 / (A + 2 * np.cos(np.deg2rad(sweep_c_4_deg))))
    C_m_ac_fus_ld = -1.8 * (1 - 2.5 * b_f / l_f) * np.pi * b_f * h_f * l_f / (4 * S * mac) * C_L_0 / C_L_alpha_A_less_h_LD
    C_m_ac_flap_ld = mu2 * (-mu1 * Delta_Cl_max * ext_flap_chord_ratio - (C_L_LD + Delta_Cl_max * (1 - Swf / S)) / 8 * ext_flap_chord_ratio * (ext_flap_chord_ratio - 1)) + 0.7 * A / (1 + 2 / A) * mu3 * Delta_Cl_max * np.tan(np.deg2rad(sweep_c_4_deg))
    C_m_ac_ld = C_m_ac_w + C_m_ac_flap_ld + C_m_ac_fus_ld

    # Tail sizes
    Sh_S_cont = (x_cg_lemac_mac - (x_ac_ld - C_m_ac_ld / C_L_A_less_h_ld)) / (C_L_h / C_L_A_less_h_ld * l_h / mac * V_h_V2)
    Sh_S_n_stab = (x_cg_lemac_mac - x_ac_cr) / (C_L_alpha_h_cr / C_L_alpha_A_less_h_cr * (1 - depsilon_dalpha_cr) * l_h / mac * V_h_V2)
    Sh_S_stab = (x_cg_lemac_mac - (x_ac_cr - SM)) / (C_L_alpha_h_cr / C_L_alpha_A_less_h_cr * (1 - depsilon_dalpha_cr) * l_h / mac * V_h_V2)

    # PLotting
    min_x = min(x_cg_lemac_mac)
    max_x = max(x_cg_lemac_mac)
    range_x = np.abs(max_x - min_x)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(min_x - 0.1 * range_x, max_x + 0.1 * range_x)
    ax.set_ylim(0, 0.7)
    ax.plot(x_cg_lemac_mac, Sh_S_cont, color='green', label='Controlability')
    ax.plot(x_cg_lemac_mac, Sh_S_n_stab, color='grey', label='Neutral stability', linestyle='--')
    ax.plot(x_cg_lemac_mac, Sh_S_stab, color='purple', label=f'Stability (SM = {SM})')
    ax.set_xlabel(f"x$_{{cg}}$/mac from lemac")
    ax.set_ylabel(f"$\\frac{{S_h}}{{S}}$")
    ax.set_title("Scissor plot")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    if output_filepath is not None:
        plt.savefig(output_filepath, dpi=300)
    if show_plot:
        plt.show()
    print(f'Scissor plot saved to {output_filepath}')
    return Sh_S_cont, Sh_S_n_stab, Sh_S_stab

def overlay_wing_pos_and_scissor_plot(ac: Aircraft, 
                                      x_le_w_fus_length_arr: np.ndarray,
                                      output_filepath: str = None,
                                      show_plot: bool = False):
    x_cg_lemac_mac_plot = convert_x_cg_from_nose_to_lemac_frac_mac(x_le_w_fus_length_arr * ac.fuselage.length, ac)

    fwd_cg, aft_cg, ac = loading_diagram(x_le_w_fus_length_arr * ac.fuselage.length, ac)
    x_le_w_l_fus = x_le_w_fus_length_arr
    Sh_S_cont, Sh_S_n_stab, Sh_S_stab = scissor_plot(ac, x_cg_lemac_mac_plot)

    # Plotting

    # Shared x-axis
    x = np.linspace(0, 1, 100)

    # Two different relationships
    y1 = -8 * (x - 0.3) + 0.5      # Example left-axis curve
    y2 = 2 * (x - 0.6)**2 + 0.4    # Example right-axis curve

    # Create figure and first axis
    fig, ax1 = plt.subplots(figsize=(8, 6))

    # Scissors 
    line1a = ax1.plot(x_cg_lemac_mac_plot, Sh_S_cont, color='navy', label='Controlability')
    line1b = ax1.plot(x_cg_lemac_mac_plot, Sh_S_n_stab, color='navy', label='Neutral stability')
    line1c = ax1.plot(x_cg_lemac_mac_plot, Sh_S_stab, color='navy', label=f'Stability (SM=0.05)')
    ax1.set_xlabel(r'$X_{cg}/MAC$')
    ax1.set_ylabel(r'$S_h/S$', color='navy')
    ax1.tick_params(axis='y', labelcolor='navy')

    # Create second y-axis sharing x-axis
    ax2 = ax1.twinx()

    # Right y-axis plot
    line2a = ax2.plot(fwd_cg, x_le_w_l_fus, color='steelblue', label='Forward cg')
    line2b = ax2.plot(aft_cg, x_le_w_l_fus, color='steelblue', label='Aft cg')
    ax2.set_ylabel(r'$X_{LE}/L_{fus}$', color='steelblue')
    ax2.tick_params(axis='y', labelcolor='steelblue')

    # Example intersection markers / annotations
    x_point = 0.4
    y1_point = np.interp(x_point, x, y1)
    y2_point = np.interp(x_point, x, y2)

    intersection1 = scissor_plot_intersection_points(x1=x_cg_lemac_mac_plot, x2=fwd_cg, y1=Sh_S_cont, y2=x_le_w_l_fus)
    intersection2 = scissor_plot_intersection_points(x1=x_cg_lemac_mac_plot, x2=aft_cg, y1=Sh_S_stab, y2=x_le_w_l_fus)
    if intersection1 is not None and intersection2 is not None:
        ax1.plot(intersection1[0], intersection1[1], 'ko')
        ax1.plot(intersection2[0], intersection2[1], 'ko')
        ax1.axhline(max(intersection1[1], intersection2[1]), linestyle='--', color='gray', alpha=0.7)

    plt.title('Tail sizing and wing position plot')
    plt.savefig(output_filepath)
    if show_plot:
        plt.show()
    
    if intersection1 is not None and intersection2 is not None:
        if intersection1[1] > intersection2[1]:
            Sh_S = intersection1[1]
            x_cg_lemac_mac = intersection1[0]
        else:
            Sh_S = intersection2[1]
            x_cg_lemac_mac = intersection2[0]
        wing_pos = float(np.interp(x_cg_lemac_mac, x_cg_lemac_mac_plot, x_le_w_l_fus))
        return Sh_S, wing_pos
    else: 
        return None





