from classes.aircraft_2 import Aircraft, loader, Requirements, Mission, Fuselage, Wing, Engine, Weights, Empennage, HLD_and_AIL, Landing_Gear
from lookups.consts import *
import pandas as pd
import numpy as np
from scipy.interpolate import RegularGridInterpolator, interp1d
from scipy.optimize import brentq
from classes.isa import Atmosphere
import matplotlib.pyplot as plt

# NOTE: Check how to calculate W_ops for medivac, check fti mass, check if ballasts needed
# NOTE: Check calculation for x_cg cargo hold 


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

def closest_value(x, values=[2, 4, 6, 8, 10]):
    return min(values, key=lambda v: abs(v - x))

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

def x_le_from_x_lemac(x_lemac, y_mac, le_sweep_deg):
    return x_lemac - y_mac * np.tan(np.deg2rad(le_sweep_deg))

def Snet(S, b_f, taper, b, c_r):
    c_fus_int = chord_at_y_span(c_r, taper, b_f/2, b)
    return S - (c_r + c_fus_int) * b_f / 2

def S_wf(y_start_f, y_end_f, taper, c_r, b):
    c_start_f = chord_at_y_span(c_r, taper, y_start_f, b)
    c_end_f = chord_at_y_span(c_r, taper, y_end_f, b)
    return (c_end_f + c_start_f) * np.abs(y_end_f - y_start_f)

def beta(V, altitude, temp_shift):
    M = V / np.sqrt(287 * 1.4 * float(Atmosphere(altitude, temp_shift).temp))
    return np.sqrt(1 - M**2)

def lift_slope(A, beta, sweep_c_2_deg, eta=0.95):
    return 2 * np.pi * A / (2 + np.sqrt(4 + (A * beta / eta)**2 * (1 + np.tan(np.deg2rad(sweep_c_2_deg))**2 / beta**2)))

def D_Cl_max(flap_type: str, # 'plain' or 'split' or 'slotted' or 'fowler' or 'double slotted'or 'triple slotted' or 'fowler'
             cdash_c: float) -> float:
    flap_values = {
        'plain': 0.9,
        'split': 0.9,
        'slotted': 1.3,
        'fowler': 1.3 * cdash_c,
        'double slotted': 1.6 * cdash_c,
        'triple slotted': 1.9 * cdash_c,
        'fixed slot': 0.2,
        'leading edge flap': 0.3,
        'kruger': 0.3,
        'slat': 0.4 * cdash_c
    }

    return flap_values.get(flap_type)

def W_to(ac: Aircraft, w_oe, w_f, w_pl, w_crew = 0, update_ac: bool = False):
    Wto = sum(w_oe, w_f, w_pl, w_crew) * (1 + 0.005 / 0.995)  # 0.5% trapped fuel and oil
    if update_ac:
        ac.weights.m_takeoff = Wto
    return Wto

def W_oe_and_cg_from_nose1(ac: Aircraft, update_ac: bool = False, 
                          pie_chart_output_path: str = None, show_pie_chart: bool = False, 
                          struc_pie_chart_output_path: str = None, struc_show_pie_chart: bool = False) -> tuple:
    w_power, x_cg_pwr = W_pwr_and_x_cg_from_nose(ac)
    w_mlg, w_nlg = W_gear(ac)
    w_ht, w_vt = W_emp(ac)
    wwing = W_wing(ac)
    wfus = W_fus(ac)
    wnac = W_nac(ac)
    w_structure = wwing + w_ht + w_vt + wfus + wnac + w_mlg + w_nlg
    w_fxeq, x_cg_fxeq = W_feq_and_cg_from_nose(ac)
    W_oe = w_structure + w_power + w_fxeq
    x_cg_oe = (w_structure * x_cg_structural_from_nose(ac, update_ac=False)[0] + w_fxeq * x_cg_fxeq + w_power * x_cg_pwr)
    ac.weights.oew_frac = W_oe / ac.weights.m_takeoff
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

def W_oe_and_cg_from_nose2(ac: Aircraft, update_ac: bool = False, 
                          pie_chart_output_path: str = None, show_pie_chart: bool = False, 
                          struc_pie_chart_output_path: str = None, struc_show_pie_chart: bool = False) -> tuple:
    w_power, x_cg_pwr = W_pwr_and_x_cg_from_nose(ac)
    w_mlg, w_nlg = W_gear(ac)
    w_ht, w_vt = W_emp(ac)
    wwing = W_wing(ac)
    wfus = W_fus(ac)
    wnac = W_nac(ac)
    w_structure = wwing + w_ht + w_vt + wfus + wnac + w_mlg + w_nlg
    w_fxeq, x_cg_fxeq = W_feq_and_cg_from_nose(ac)
    W_oe = w_structure + w_power + w_fxeq
    x_cg_oe = (w_structure * x_cg_structural_from_nose(ac, update_ac=False)[0] + w_fxeq * x_cg_fxeq + w_power * x_cg_pwr)
    ac.weights.oew_frac = W_oe / ac.weights.m_takeoff
    if update_ac:
        ac.weights.m_empty = W_oe
        ac.weights.x_cg_oew = x_cg_oe
    if pie_chart_output_path is not None:
        categories = ['Structural', 'Power', 'Fixed equipment']
        raw_values = [w_structure, w_power, w_fxeq]
        values = np.array(raw_values) / W_oe * 100
        fig, ax = plt.subplots(figsize=(8, 6))
        wedges, _ = ax.pie(values, startangle=90)
        total = sum(values)
        for i, wedge in enumerate(wedges):
            angle = (wedge.theta2 + wedge.theta1) / 2
            x = np.cos(np.radians(angle))
            y = np.sin(np.radians(angle))
            percentage = values[i] / total * 100

            ax.annotate(
                f'{categories[i]}\n{percentage:.1f}%\n({raw_values[i]:.1f} kg)',
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
        raw_values = [wwing, w_ht, w_vt, wfus, wnac, w_mlg+w_nlg]
        values = np.array(raw_values) / w_structure * 100
        fig, ax = plt.subplots(figsize=(8, 6))
        wedges, _ = ax.pie(values, startangle=90)
        total = sum(values)
        for i, wedge in enumerate(wedges):
            angle = (wedge.theta2 + wedge.theta1) / 2
            x = np.cos(np.radians(angle))
            y = np.sin(np.radians(angle))
            percentage = values[i] / total * 100

            ax.annotate(
                f'{categories[i]}\n{percentage:.1f}%\n({raw_values[i]:.1f} kg)',
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

def W_to_new2(ac: Aircraft,
             m_res = 0.0, # from class I
             W_crew = 0.0, # included in PL probably
             update_ac: bool = False,
             pie_chart_output_path: str = None,
             show_pie_chart: bool = False
             ):
    W_PL = ac.weights.m_payload
    W_e = ac.weights.m_empty
    print(f' New oew: {W_e}')
    m_ff = 1 - ac.weights.m_fuel / ac.weights.m_takeoff # Insert class I method called
    m_tfo = ac.weights.m_takeoff * 0.005
    Wto = W_e / ac.weights.oew_frac
    Wto = W_e + ac.weights.m_fuel + W_PL + m_tfo
    print(f'Wto: {Wto}, m_ff: {m_ff}, W_e: {W_e}, W_PL: {W_PL}, W_fuel: {ac.weights.m_fuel}')
    if update_ac:
        ac.weights.m_takeoff = Wto

    if pie_chart_output_path is not None:
        categories = ['Empty weight', 'Payload', 'Trapped fuel and oil', 'Fuel']
        raw_values = [W_e, W_PL, m_tfo, Wto-W_e-m_tfo-W_PL]
        values = np.array(raw_values) / Wto * 100
        print(f'Pie chart values: {values}')
        fig, ax = plt.subplots(figsize=(8, 6))
        wedges, _ = ax.pie(values, startangle=90)
        total = sum(values)
        for i, wedge in enumerate(wedges):
            angle = (wedge.theta2 + wedge.theta1) / 2
            x = np.cos(np.radians(angle))
            y = np.sin(np.radians(angle))
            percentage = values[i] / total * 100

            ax.annotate(
                f'{categories[i]}\n{percentage:.1f}%\n({raw_values[i]:.1f} kg)',
                xy=(x, y),
                xytext=(0.9 * np.sign(x), 0.9 * y),
                arrowprops=dict(arrowstyle='->'),
                ha='left' if x > 0 else 'right'
            )

        ax.axis('equal')
        plt.title('Distribution of MTOW')
        plt.savefig(pie_chart_output_path)
        if show_pie_chart:
            plt.show()
    return Wto, ac

def W_oe_and_cg_from_nose(ac: Aircraft, update_ac: bool = False,
                          pie_chart_output_path: str = None, show_pie_chart: bool = False,
                          struc_pie_chart_output_path: str = None, struc_show_pie_chart: bool = False) -> tuple:
    w_power, x_cg_pwr = W_pwr_and_x_cg_from_nose(ac)
    w_mlg, w_nlg = W_gear(ac)
    w_ht, w_vt = W_emp(ac)
    wwing = W_wing(ac)
    wfus = W_fus(ac)
    wnac = W_nac(ac)
    w_structure = wwing + w_ht + w_vt + wfus + wnac + w_mlg + w_nlg
    w_fxeq, x_cg_fxeq = W_feq_and_cg_from_nose(ac)
    W_oe = w_structure + w_power + w_fxeq
    x_cg_oe = (w_structure * x_cg_structural_from_nose(ac, update_ac=False)[0] + w_fxeq * x_cg_fxeq + w_power * x_cg_pwr)
    ac.weights.oew_frac = W_oe / ac.weights.m_takeoff
    if update_ac:
        ac.weights.m_empty = W_oe
        ac.weights.x_cg_oew = x_cg_oe
    if pie_chart_output_path is not None:
        categories = ['Structural', 'Power', 'Fixed equipment']
        raw_values = [w_structure, w_power, w_fxeq]
        values = np.array(raw_values) / W_oe * 100
        fig, ax = plt.subplots(figsize=(8, 6))
        wedges, _ = ax.pie(values, startangle=90)
        total = sum(values)
        for i, wedge in enumerate(wedges):
            angle = (wedge.theta2 + wedge.theta1) / 2
            x = np.cos(np.radians(angle)) * 0.7
            y = np.sin(np.radians(angle)) * 0.7
            percentage = values[i] / total * 100
            ax.text(x, y,
                    f'{categories[i]}\n{percentage:.1f}%\n({raw_values[i]:.1f} kg)',
                    ha='center', va='center', fontsize=8)
        ax.axis('equal')
        plt.title('Distribution of OEW')
        plt.savefig(pie_chart_output_path, dpi=400)
        if show_pie_chart:
            plt.show()
    if struc_pie_chart_output_path is not None:
        categories = ['Wing', 'Horizontal tail', 'Vertical tail', 'Fuselage', 'Nacelles', 'Landing gear']
        raw_values = [wwing, w_ht, w_vt, wfus, wnac, w_mlg+w_nlg]
        values = np.array(raw_values) / w_structure * 100
        fig, ax = plt.subplots(figsize=(8, 6))
        wedges, _ = ax.pie(values, startangle=90)
        total = sum(values)
        for i, wedge in enumerate(wedges):
            angle = (wedge.theta2 + wedge.theta1) / 2
            x = np.cos(np.radians(angle)) * 0.7
            y = np.sin(np.radians(angle)) * 0.7
            percentage = values[i] / total * 100
            ax.text(x, y,
                    f'{categories[i]}\n{percentage:.1f}%\n({raw_values[i]:.1f} kg)',
                    ha='center', va='center', fontsize=8)
        ax.axis('equal')
        plt.title('Distribution of structural weight')
        plt.savefig(struc_pie_chart_output_path, dpi=400)
        if struc_show_pie_chart:
            plt.show()
    return W_oe, x_cg_oe, ac

def W_to_new(ac: Aircraft,
             m_res = 0.0,
             W_crew = 0.0,
             update_ac: bool = False,
             pie_chart_output_path: str = None,
             show_pie_chart: bool = False
             ):
    W_PL = ac.weights.m_payload
    W_e = ac.weights.m_empty
    print(f' New oew: {W_e}')
    m_ff = 1 - ac.weights.m_fuel / ac.weights.m_takeoff
    m_tfo = ac.weights.m_takeoff * 0.005
    # Wto = W_e / ac.weights.oew_frac
    Wto = W_e + ac.weights.m_fuel + W_PL + m_tfo
    print(f'Wto: {Wto}, m_ff: {m_ff}, W_e: {W_e/Wto}, W_PL: {W_PL}, W_fuel: {ac.weights.m_fuel}, W_fuel/mtow: {ac.weights.m_fuel / Wto}')
    if update_ac:
        ac.weights.m_takeoff = Wto
        ac.weights.oew_frac = W_e / Wto
    if pie_chart_output_path is not None:
        categories = ['Empty weight', 'Payload', 'Trapped fuel and oil', 'Fuel']
        raw_values = [W_e, W_PL, m_tfo, ac.weights.m_fuel]
        values = np.array(raw_values) # / Wto * 100
        print(f'Pie chart values: {values}')
        fig, ax = plt.subplots(figsize=(8, 6))
        wedges, _ = ax.pie(values, startangle=90)
        total = sum(values)
        for i, wedge in enumerate(wedges):
            angle = (wedge.theta2 + wedge.theta1) / 2
            x = np.cos(np.radians(angle)) * 0.7
            y = np.sin(np.radians(angle)) * 0.7
            if categories[i]=='Trapped fuel and oil':
                x = np.cos(np.radians(angle))
            percentage = values[i] / total * 100
            ax.text(x, y,
                    f'{categories[i]}\n{percentage:.1f}%\n({raw_values[i]:.1f} kg)',
                    ha='center', va='center', fontsize=8)
        ax.axis('equal')
        plt.title('Distribution of MTOW')
        plt.savefig(pie_chart_output_path, dpi=400)
        if show_pie_chart:
            plt.show()
    return Wto, ac

def W_to_new1(ac: Aircraft,
             m_res = 0.0, # from class I
             W_crew = 0.0, # included in PL probably
             update_ac: bool = False,
             pie_chart_output_path: str = None,
             show_pie_chart: bool = False
             ):
    W_PL = ac.weights.m_payload
    W_e = ac.weights.m_empty
    print(f' New oew: {W_e}')
    m_ff = 1 - ac.weights.m_fuel / ac.weights.m_takeoff # Insert class I method called
    m_tfo = ac.weights.m_takeoff * 0.005
    # Wto = (W_e + W_PL + W_crew) / (m_ff * (1 + m_res) - m_res - m_tfo)
    Wto = W_e / ac.weights.oew_frac
    Wto = W_e + ac.weights.m_fuel + W_PL + m_tfo
    print(f'Wto: {Wto}, m_ff: {m_ff}, W_e: {W_e}, W_PL: {W_PL}, W_fuel: {ac.weights.m_fuel}')
    if update_ac:
        ac.weights.m_takeoff = Wto

    if pie_chart_output_path is not None:
        categories = ['Empty weight', 'Payload', 'Trapped fuel and oil', 'Fuel']
        values = [W_e, W_PL, m_tfo, Wto-W_e-m_tfo-W_PL] / Wto * 100
        print(f'Pie chart values: {values}')
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
    return Wto, ac

def W_pwr_and_x_cg_from_nose(ac: Aircraft):
    x_le = ac.wing.x_le
    alpha_p_id = ac.engine.alpha_p_id
    P_to = ac.engine.power_to / HP_TO_W
    W_piston = ac.weights.m_piston
    W_supercap = ac.weights.m_supercap
    W_eng = ac.weights.m_turboprop / LBS_TO_KG
    if alpha_p_id == 'piston':
        W_eng = W_piston / LBS_TO_KG
    N_e = ac.engine.count
    W_fuel = ac.weights.m_fuel / LBS_TO_KG
    Fuel_type = ac.engine.fuel_type  # Check if avgas or Jet-A1 or other
    N_t = ac.engine.n_fuel_tanks # Number of separate fuel tanks
    INT = 1.0  # Fraction of fuel tanks that are intergral
    # W_pwr = W_eng + W_ai + W_prop + W_fs + W_p + W_batt = W_batt + W_pwr1 + W_fs (fuel system)
    
    # USAF - Roskam eqn 6.3
    # W_pwr1a = 2.075 * W_eng**0.922 * N_e
    # W_pwr1b = N_e * ac.weights.m_propeller / LBS_TO_KG + W_eng + 1.03 * N_e**0.3 * P_to**0.7
    # print(f' \n Power system weights: USAF:{W_pwr1a} Torenbeek: {W_pwr1b}')
    # W_pwr1 = (W_pwr1a + W_pwr1b) / 2

    W_pwr1 = N_e * ac.weights.m_propeller / LBS_TO_KG + W_eng + 1.03 * N_e**0.3 * P_to**0.7

    if alpha_p_id != 'hydrogen':
        if Fuel_type == 'avgas':
            K_fsp = 5.87  # lbs/gal
        elif Fuel_type == 'Jet-A1':
            K_fsp = 6.55  # lbs/gal
        else:
            raise ValueError(f"Fuel type given: {Fuel_type}, should be either 'avgas' or 'Jet-A1'")
        W_fs = 2.49 * ((W_fuel / K_fsp)**0.6 * (1 / (1 + INT))**0.3 * N_t**0.2 * N_e**0.13)**1.21
    
    # elif alpha_p_id == 'hydrogen':  # NOTE: insert method of fixed equipment weight
    #     W_fs = ...
    
    if ac.engine.eng_x_pos == 'le':
        nac_y = ac.fuselage.width / 2 + ac.engine.eng_y_pos_fuselage
        x_cg_pwr1 = x_pos_le_along_span_from_nose(ac.wing.sweep_LE_deg, nac_y, x_le)
    x_cg = (W_pwr1 * x_cg_pwr1 + W_fs * (x_le + ac.wing.c_root * ac.engine.x_cg_fuel_tanks_c_r) + W_supercap * x_le) / (W_pwr1 + W_supercap + W_fs)
    ac.weights.power_system = (W_pwr1 + W_supercap + W_fs) * LBS_TO_KG
    # print(f' \nPower system weight: {(W_pwr1 + W_supercap + W_fs) * LBS_TO_KG}kg \n')
    return (W_pwr1 + W_supercap + W_fs) * LBS_TO_KG, x_cg

def W_feq_and_cg_from_nose(ac: Aircraft):
    Wto = ac.weights.m_takeoff / LBS_TO_KG
    W_e = ac.weights.m_empty
    N_pax = ac.fuselage.n_pax  # Including crew
    T_cr = float(Atmosphere(ac.requirements.cruise['cr_altitude'] * FT_TO_M).temp_isa)
    M_D = 1.5 * ac.requirements.cruise['cr_speed'] * KTS_TO_MS / np.sqrt(1.4 * 287 * T_cr)  # Design dive Mach number
    no_pressurization_const = 0.6 # Fraction to take into account api does not have p (pressurization) in our case
    V_pax_cargo = ac.fuselage.vol_cabin_and_cargo * M2_TO_F2/FT_TO_M  # Volume of passenger cabin and cargo [ft3]

    W_fc = 0.5 * (0.0168 * Wto + 1.066 * Wto**0.626)

    # Cessna
    W_hps = 0.009 * Wto
    W_els = 0.0268 * Wto
    # # Torenbeek (W_hps + W_els = W_hps_els)
    # W_hps_els = 0.0078 * Wto**1.2
    W_iae = 40 + 0.008 * Wto
    W_api = (0.265 * Wto**0.52 * N_pax**0.68 * W_iae**0.17 * M_D**0.08) * no_pressurization_const
    W_fur = 0.5 * (0.412 * N_pax**1.145 * Wto**0.489 + 15 * N_pax + V_pax_cargo)
    W_ops = 0
    W_fti = 0.5 * (155 / 9980 * Wto + 708 / 24912 * Wto)
    W_aux = 0.01 * W_e
    W_bal = 0
    W_pt = 0.0045 * Wto
    W_etc = 0  # NOTE: check if there are other things to include

    # CG from nose - *indicated source = https://archive.aoe.vt.edu/mason/Mason_f/M96SC02.pdf 
    c_r_w = ac.wing.c_root
    taper_w = ac.wing.taper_ratio
    b_w = ac.wing.span
    b_ht = ac.empennage.horizontal_tail['b_h']
    b_vt = ac.empennage.vertical_tail['b_v']
    w_sweep_c_4_deg = ac.wing.sweep
    ht_sweep_le_deg = ac.empennage.horizontal_tail['sweep_LE_deg']
    vt_sweep_le_deg = ac.empennage.vertical_tail['sweep_LE_deg']
    x_le_w = ac.wing.x_le  # Wing root chord LE distance from nose
    x_le_ht = ac.empennage.horizontal_tail['x_h_frac_lf'] * ac.fuselage.length - ac.empennage.horizontal_tail['MAC_h']* np.tan(np.deg2rad(ht_sweep_le_deg)) - 0.4 * ac.empennage.horizontal_tail['y_MAC_h'] # HT root chord LE distance from nose
    x_le_vt = ac.empennage.vertical_tail['x_v_frac_lf'] * ac.fuselage.length - ac.empennage.vertical_tail['MAC_v']* np.tan(np.deg2rad(vt_sweep_le_deg)) - 0.4 * ac.empennage.vertical_tail['y_MAC_v']  # VT root chord LE distance from nose
    x_c_front_spar = ac.wing.x_c_front_spar
    x_c_rear_spar = ac.wing.x_c_rear_spar
    start_flap_along_span = ac.hld_and_ailerons.flaps['y_flap_in']
    end_aileron_along_span = ac.hld_and_ailerons.ailerons['y_aileron_out']
    l_fus = ac.fuselage.length
    l_tc = ac.fuselage.tail_cone_length  # length of tailcone
    l_cabin = ac.fuselage.l_cabin
    start_cabin = ac.fuselage.start_cabin 
    x_nlg = ac.landing_gear.longitudinal_nlg

    chord_fc = (chord_at_y_span(c_r_w, taper_w, start_flap_along_span, b_w) + chord_at_y_span(c_r_w, taper_w, end_aileron_along_span, b_w)) / 2
    # print(f'Flight controls params: {start_flap_along_span, end_aileron_along_span, x_pos_le_along_span_from_nose(ac.wing.sweep_LE_deg, (start_flap_along_span+end_aileron_along_span)/2, x_le_w)}')
    y_fc = x_pos_le_along_span_from_nose(ac.wing.sweep_LE_deg, (start_flap_along_span+end_aileron_along_span)/2, x_le_w)
    x_cg_fc = (x_c_rear_spar + 1) / 2 * chord_fc + y_fc  # *Between aft spar and trailing-edge
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

    # print(f'Weights FEQ: {W_fc, W_hps, W_els, W_iae, W_api, W_fur, W_ops, W_fti, W_aux, W_pt, W_etc}')
    Weights = np.array([W_fc, W_hps, W_els, W_iae, W_api, W_fur, W_ops, W_fti, W_aux, W_pt, W_etc])
    W_feq = np.sum(Weights)
    cgs = np.array([x_cg_fc, x_cg_hps, x_cg_els, x_cg_iae, x_cg_api, x_cg_fur, x_cg_ops, x_cg_fti, x_cg_aux, x_cg_pt, x_cg_etc])
    x_cg_feq = (Weights@cgs) / W_feq

    return W_feq * LBS_TO_KG, x_cg_feq

def W_wing(ac: Aircraft, update_ac: bool = False):
    w = ac.wing
    Wto = ac.weights.m_takeoff / LBS_TO_KG
    S = w.area * M2_TO_F2
    n_ult = ac.requirements.general['n_ult']  # NOTE: add later
    A = w.aspect_ratio
    t_c_max = w.t_c_max
    # sweep_c_4_deg = w.sweep
    # taper = w.taper_ratio
    # V_h = ...  # NOTE: add later: maximum level speed at sealevel in kts
    wing_type = w.wing_type  # 'cantilever' or not
    b = w.span / FT_TO_M
    c_r = w.c_root / FT_TO_M
    t_r = t_c_max * c_r  # max thickness of wing root chord in ft
    
    # Cessna
    if wing_type == 'cantilever':
        W_wing_c = 0.04674 * (Wto ** 0.397) * (S ** 0.36) * (n_ult ** 0.397) * (A ** 1.712)

    else:
        W_wing_c = 0.002933 * (S ** 1.018) * (A ** 2.473) * (n_ult ** 0.611)

    # USAF
    # W_wing_u = 96.948 * (((Wto * n_ult * 1e-5) ** 0.65) * ((A / np.cos(np.deg2rad(sweep_c_4_deg))) ** 0.57) * ((S / 100) ** 0.61) * (((1 + taper) / (2 * t_c_max)) ** 0.36) * ((1 + V_h / 500) ** 0.5)) ** 0.993

    # Torenbeek
    W_wing_t = 0.00125 * Wto * ((b / np.cos(np.deg2rad(w.sweep_c_2_deg))) ** 0.75) * (1 + (6.3 * np.cos(np.deg2rad(w.sweep_c_2_deg)) / b) ** 0.5) * (n_ult ** 0.55) * (b * S / (t_r * Wto * np.cos(np.deg2rad(w.sweep_c_2_deg))))**0.3
    # print(f'\n Wing weight: \t Cessna: {W_wing_c * LBS_TO_KG} \t Torenbeek: {W_wing_t * LBS_TO_KG}')
    return W_wing_t * LBS_TO_KG
    # return (W_wing_c + W_wing_t) / 2 * LBS_TO_KG

def W_emp(ac: Aircraft, update_ac: bool = False):
    ht = ac.empennage.horizontal_tail
    vt = ac.empennage.vertical_tail

    n_ult = ac.requirements.general['n_ult']  # NOTE: add later
    S_h = ht['area'] * M2_TO_F2
    S_v = vt['area'] * M2_TO_F2
    l_h = ac.empennage.horizontal_tail['l_h'] / FT_TO_M
    b_h = ht['b_h'] / FT_TO_M
    b_v = vt['b_v'] / FT_TO_M
    t_c_max_h = ac.empennage.horizontal_tail['t_c_max']
    t_c_max_v = ac.empennage.vertical_tail['t_c_max']
    c_r_h = ht['c_r_h'] / FT_TO_M
    c_r_v = vt['c_r_v'] / FT_TO_M
    t_r_h = t_c_max_h * c_r_h
    t_r_v = t_c_max_v * c_r_v
    A_h = ht['aspect_ratio']
    A_v = vt['aspect_ratio']
    vt_sweep_c_4_deg = ac.empennage.vertical_tail['sweep']
    Wto = ac.weights.m_takeoff / LBS_TO_KG

    # Cessna
    W_h_c = (3.184 * (Wto ** 0.887) * (S_h ** 0.101) * (A_h ** 0.138)) / (174.04 * (t_r_h **0.223))
    W_v_c = (1.68 * (Wto ** 0.567) * (S_v ** 1.249) * (A_v ** 0.482)) / (639.95 * (t_r_v ** 0.747) * (np.cos(np.rad2deg(vt_sweep_c_4_deg)) ** 0.882))
    W_emp_c = W_h_c + W_v_c

    # USAF
    W_h_u = 127 * ((Wto * n_ult * 1e-5)**0.87 * (S_h / 100)**1.2 * 0.289 * (l_h / 10)**0.483 * (b_h / t_r_h)**0.5)**0.458
    W_v_u = 98.5 * ((Wto * n_ult * 1e-5)**0.87 * (S_v / 100)**1.2 * 0.289 * (b_v / t_r_v)**0.5)**0.458
    W_emp_u = W_h_u + W_v_u

    # # Torenbeek
    # W_emp_t = 0.04 * (n_ult * (S_v + S_h)**2)**0.75

    W_ht = (W_h_c + W_h_u) / 2
    # W_vt = (W_v_c + W_v_u) / 2
    W_vt = (0.5*W_v_c + 1.5*W_v_u) / 2
    # print(f'\n HT weight: \t Cessna: {W_h_c * LBS_TO_KG} \t USAF: {W_h_u * LBS_TO_KG}')
    # print(f'\n VT weight: \t Cessna: {W_v_c * LBS_TO_KG} \t USAF: {W_v_u * LBS_TO_KG}')

    return W_ht * LBS_TO_KG, W_vt * LBS_TO_KG

def W_fus(ac: Aircraft, update_ac: bool = False):
    fus = ac.fuselage
    l_f = fus.length / FT_TO_M  # Fuselage length in ft
    w_f = fus.width/ FT_TO_M  # Max fuselage width in ft
    h_f = fus.height / FT_TO_M  # Max fuselage height in ft
    n_ult = ac.requirements.general['n_ult']
    alt_c = ac.requirements.cruise['cr_altitude'] * FT_TO_M  # m
    rho_c = float(Atmosphere(height=alt_c).density) # kg/m3
    V_c = ac.requirements.cruise['cr_speed'] * np.sqrt(rho_c / 1.225) # KEAS
    l_f_n = l_f # Fuselage length - nacelle in ft
    P_max = ac.fuselage.max_perimeter / FT_TO_M  # Maximum fuselage perimeter in ft
    N_pax = ac.fuselage.n_pax  # Number of passengers including pilot(s)
    Wto = ac.weights.m_takeoff / LBS_TO_KG

    # Cessna
    W_fus_c = 14.86 * (Wto**0.144) * (l_f_n / P_max)**0.778 * (l_f_n**0.383) * (N_pax**0.455)
    
    # USAF
    W_fus_u = 200 * ((Wto * n_ult * 1e-5)**0.286 * (l_f / 10)**0.857 * ((w_f + h_f) / 10) * (V_c / 100)**0.338)**1.1
    # print(f'\n fus weight: \t Cessna: {W_fus_c * LBS_TO_KG} \t USAF: {W_fus_u * LBS_TO_KG}')
    return (W_fus_c + W_fus_u) / 2 * LBS_TO_KG

def W_nac(ac: Aircraft, update_ac: bool = False):
    alpha_p_id = ac.engine.alpha_p_id  # turboprop or hydrogen or piston
    engine_over_wing: bool = ac.engine.eng_above_wing
    P_to = ac.engine.power_to / HP_TO_W  # HP

    # Cessna NOTE: fill in methods
    if alpha_p_id == 'piston':
        engine_type = 'horizontally opposed'
        if engine_type == 'radial':
            K_n = 0.37   # lbs/hp 
        elif engine_type == 'horizontally opposed':
            K_n = 0.24   # lbs/hp 

        W_n_c = K_n * P_to
    # elif alpha_p_id == 'hydrogen':
    #     W_n_c = ...
    # elif alpha_p_id =='turboprop':
    #     W_n_c = ...

    # Torenbeek
    if alpha_p_id == 'piston':
        engine_type = 'horizontally opposed'
        if engine_type == 'radial':
            N_e = ac.engine.count
            W_n_t = 0.045 * P_to**1.25 * N_e**(-0.25)
        elif engine_type == 'horizontally opposed':
            W_n_t = 0.32 * P_to 
        if engine_over_wing:
            W_n_t += 0.11 * P_to
        # print(f'\n nac weight: \t Cessna: {W_n_c * LBS_TO_KG} \t Torenbeek: {W_n_t * LBS_TO_KG}')
        # return (W_n_c + W_n_t) / 2 * LBS_TO_KG
        return W_n_c * LBS_TO_KG

    # elif alpha_p_id == 'hydrogen':
    #     W_n_t = ...
    elif alpha_p_id =='turboprop':
        W_n_t = 0.14 * P_to 
        if engine_over_wing:
            W_n_t += 0.11 * P_to
        return W_n_t * LBS_TO_KG

def W_gear(ac: Aircraft, update_ac: bool = False):
    Wto = ac.weights.m_takeoff / LBS_TO_KG
    W_L = ac.requirements.landing['la_mass_frac'] * Wto
    shock_strut_frac_whole = 0.2
    d_tire = ac.landing_gear.selected_mlg_tire["Outside Diameter Max (In)"] * 2.54 / 100  # tire diameter
    d_tire_n = ac.landing_gear.selected_nlg_tire["Outside Diameter Max (In)"] * 2.54 / 100  # tire diameter
    l_s_m = shock_strut_frac_whole * (np.abs(ac.landing_gear.height_mlg) - d_tire / 2) / FT_TO_M # Shock strut length main gear [ft]
    l_s_n = shock_strut_frac_whole * (np.abs(ac.landing_gear.height_nlg) - d_tire_n / 2) / FT_TO_M # Shock strut length nose gear [ft]
    l_s_m = shock_strut_frac_whole * (np.abs(ac.landing_gear.height_mlg)) / FT_TO_M # Shock strut length main gear [ft]
    l_s_n = shock_strut_frac_whole * (np.abs(ac.landing_gear.height_nlg)) / FT_TO_M # Shock strut length nose gear [ft]
    n_ult = ac.requirements.general['n_ult']

    # Cessna
    W_mlg = (0.013 * Wto + 0.362 * (W_L**0.417) * (n_ult**0.950) * (l_s_m**0.183))
    # print(f'W_mlg components: {Wto, W_L, n_ult, l_s_m}')
    # print(W_mlg)
    W_nlg = (6.2 + 0.0013 * Wto + 0.007157 * (W_L**0.749) * (n_ult) * (l_s_n**0.788))
    # print(f'Landing gear: {Wto, W_L, n_ult, l_s_n, ac.landing_gear.height_nlg, d_tire_n / 2}')
    # print(f'\n LG weight: \t Cessna: mlg {W_mlg * LBS_TO_KG} \t nlg: {W_nlg * LBS_TO_KG}')
    
    # # USAF
    # W_g_u = 0.054 * (l_s_m**0.501) * (W_L * n_ult)**0.684

    return W_mlg * LBS_TO_KG, W_nlg * LBS_TO_KG

def x_cg_structural_from_nose(ac: Aircraft,
                            update_ac=False):
    w = ac.wing
    ht = ac.empennage.horizontal_tail
    vt = ac.empennage.vertical_tail

    x_le_w = w.x_le
    x_le_ht = ht['x_le']
    x_le_vt = vt['x_le']
    l_fus = ac.fuselage.length
    w_sweep_c_4_deg = ac.wing.sweep
    ht_sweep_c_4_deg = ac.empennage.horizontal_tail['sweep']
    vt_sweep_c_4_deg = ac.empennage.vertical_tail['sweep']
    w_c_r = w.c_root
    ht_c_r = ht['c_r_h']
    vt_c_r = vt['c_r_v']
    w_taper = ac.wing.taper_ratio
    ht_taper = ht['taper_ratio']
    vt_taper = vt['taper_ratio']
    b_w = ac.wing.span
    b_ht = ht['b_h']
    b_vt = vt['b_v']
    x_c_front_spar = ac.wing.x_c_front_spar
    x_c_rear_spar = ac.wing.x_c_rear_spar
    t_tail_condition: bool = ac.empennage.t_tail_condition  # is it a t-tail or not
    l_nacelle = ac.engine.length_nac
    if ac.engine.eng_x_pos == 'le':
        nac_mount_dist = 0  # distance nacelle is mounted behind le of wing (front of nacelle)
    x_nlg = ac.landing_gear.longitudinal_nlg
    x_mlg = ac.landing_gear.longitudinal_mlg
    
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
        y_pos_eng: float = ac.engine.eng_y_pos_fuselage + ac.fuselage.width / 2
        x_cg_nac = 0.4 * l_nacelle + x_pos_le_along_span_from_nose(LE_sweep_deg(w_sweep_c_4_deg, w_c_r, b_w, w_taper), y_pos_eng, x_le_w) + nac_mount_dist
    elif n_eng > 2:
        y_pos_eng: list = ac.engine.eng_y_pos_fuselage + ac.fuselage.width / 2  # NOTE: get a list/array or someting with engine y-positions for one side
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
    # print(f'Weights: {[W_wing(ac), W_ht, W_vt, W_fus(ac), W_mlg, W_nlg, W_nac(ac)]}')
    Weights = np.array([W_wing(ac), W_ht, W_vt, W_fus(ac), W_mlg, W_nlg, W_nac(ac)])
    cgs = np.array([x_cg_w, x_cg_ht, x_cg_vt, x_cg_fus, x_cg_mlg, x_cg_nlg, x_cg_nac])
    x_cg_struc = (Weights @ cgs) / np.sum(Weights)
    if update_ac:
        ac.weights.x_cg_structural = x_cg_data
    return x_cg_struc, x_cg_data, ac

def convert_x_cg_from_nose_to_lemac_frac_mac(x_le_w, x_cg_nose, ac: Aircraft):
    c_r = ac.wing.c_root
    taper = ac.wing.taper_ratio
    sweep_c_4 = ac.wing.sweep
    b = ac.wing.span
    mac = 2 * c_r / 3 * (1 + taper + taper**2) / (1 + taper)
    le_sweep_deg = LE_sweep_deg(sweep_c_4, c_r, b, taper)
    y_span_mac = (c_r - mac) / (2 / b * c_r * (1 - taper))
    x_lemac = x_pos_le_along_span_from_nose(le_sweep_deg, y_span_mac, x_le_w)
    return (x_cg_nose - x_lemac) / mac

def update_m_and_cg(m_current, cg_current, added_m, added_cg):
    new_cg = (m_current * cg_current + added_m * added_cg) / (m_current + added_m)
    new_m  =(m_current + added_m)
    return new_m, new_cg

def loading_diagram(x_le_w, ac: Aircraft, show_plot: bool=False, output_filepath=None, update_ac_cgs = False):
    c_r = ac.wing.c_root
    m_cargo = ac.weights.m_cargo
    ht_sweep_le_deg = ac.empennage.horizontal_tail['sweep_LE_deg']
    vt_sweep_le_deg = ac.empennage.vertical_tail['sweep_LE_deg']
    # x_le_ht = ac.empennage.horizontal_tail['x_h_frac_lf'] * ac.fuselage.length - ac.empennage.horizontal_tail['MAC_h']* np.tan(np.deg2rad(ht_sweep_le_deg)) - 0.4 * ac.empennage.horizontal_tail['y_MAC_h'] # HT root chord LE distance from nose
    # x_le_vt = ac.empennage.vertical_tail['x_v_frac_lf'] * ac.fuselage.length - ac.empennage.vertical_tail['MAC_v']* np.tan(np.deg2rad(vt_sweep_le_deg)) - 0.4 * ac.empennage.vertical_tail['y_MAC_v']  # ac.empennage.vertical_tail['x_v_frac_lf'] * ac.fuselage.length
    x_cg_seats = ac.fuselage.x_pos_seats  # Front to back cg positions of seats (length 3)
    assert len(x_cg_seats) == 3, f"There must be 3 seat cg positions, {len(x_cg_seats)} were given"
    x_cg_seats = convert_x_cg_from_nose_to_lemac_frac_mac(x_le_w, x_cg_seats, ac)
    n_pax = ac.fuselage.n_pax
    m_pax = ac.weights.m_pax
    n_rows = n_pax / len(x_cg_seats)
    n_window = ac.fuselage.n_window_seats
    n_aisle = ac.fuselage.n_aisle_seats
    n_middle = ac.fuselage.n_middle_seats
    m_fuel = ac.weights.m_fuel
    x_cg_fuel_tanks = ac.engine.x_cg_fuel_tanks_c_r * c_r + x_le_w
    x_cg_fuel_tanks = convert_x_cg_from_nose_to_lemac_frac_mac(x_le_w, x_cg_fuel_tanks, ac)
    if isinstance(x_cg_fuel_tanks, (list, np.ndarray, tuple)):
        # if len(x_cg_cargo_holds)>1:
        n_fuel_tanks = len(x_cg_fuel_tanks)
        m_fuel_tanks = fuel_tanks_mass_fracs * m_fuel
    else:
        fuel_tanks_mass_fracs = [1]
        m_fuel_tanks = [m_fuel]
        x_cg_fuel_tanks = [x_cg_fuel_tanks]
    m_fuel_tanks_reverse = m_fuel_tanks[::-1]
    x_cg_fuel_tanks_reverse = x_cg_fuel_tanks[::-1]

    x_cg_cargo_holds = ac.fuselage.x_cargo_holds
    x_cg_cargo_holds = convert_x_cg_from_nose_to_lemac_frac_mac(x_le_w, x_cg_cargo_holds, ac)
    cargo_mass_frac = 1.0
    if isinstance(x_cg_cargo_holds, (list, np.ndarray, tuple)):
        if len(x_cg_cargo_holds)>1:
            cargo_mass_frac = ac.fuselage.mass_frac_cargo_holds
        cargo_masses = cargo_mass_frac * m_cargo
        cargo_masses_reverse = cargo_masses[::-1]
        x_cg_cargo_reverse = x_cg_cargo_holds[::-1]
    else:
        x_cg_cargo_holds = [x_cg_cargo_holds]
        cargo_masses = [m_cargo]
        cargo_masses_reverse = cargo_masses
        x_cg_cargo_reverse = x_cg_cargo_holds

    # Initiate lists to plot:
    x_cg_lemac_plot_ftb = []
    W_plot_ftb = []
    x_cg_lemac_plot_btf = []
    W_plot_btf = []

    # OEW
    W_oe, x_cg_oe, ac = W_oe_and_cg_from_nose(ac)
    W_plot_ftb.append(W_oe)
    x_cg_lemac_plot_ftb.append(convert_x_cg_from_nose_to_lemac_frac_mac(x_le_w, x_cg_oe, ac))
    W_plot_btf.append(W_oe)
    x_cg_lemac_plot_btf.append(convert_x_cg_from_nose_to_lemac_frac_mac(x_le_w, x_cg_oe, ac))

    # Cargo/PL
    # cargo_masses = cargo_mass_frac * m_cargo
    # cargo_masses_reverse = cargo_masses[::-1]
    # x_cg_cargo_reverse = x_cg_cargo_holds[::-1]

    for i, x_cg in enumerate(x_cg_cargo_holds):
        new_m, new_cg = update_m_and_cg(W_plot_ftb[-1], x_cg_lemac_plot_ftb[-1], cargo_masses[i], x_cg)
        W_plot_ftb.append(new_m)
        x_cg_lemac_plot_ftb.append(new_cg)
        new_m, new_cg = update_m_and_cg(W_plot_btf[-1], x_cg_lemac_plot_btf[-1], cargo_masses_reverse[i], x_cg_cargo_reverse[i])
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
            new_m, new_cg = update_m_and_cg(W_plot_btf[-1], x_cg_lemac_plot_btf[-1], pax_row_mass, x_cg_seats_reverse[i])
            W_plot_btf.append(new_m)
            x_cg_lemac_plot_btf.append(new_cg)
        assert W_plot_btf[-1] == W_plot_ftb[-1], f"After loading all passengers total weights should match. \n ftb weight = {W_plot_ftb[-1]} \t btf weight = {W_plot_btf[-1]}"
        assert x_cg_lemac_plot_btf[-1] == x_cg_lemac_plot_ftb[-1], f"After loading all passengers total cg should match. \n ftb cg = {x_cg_lemac_plot_ftb[-1]} \t btf cg = {x_cg_lemac_plot_btf[-1]}"
    
    # Fuel
    for i, x_cg in enumerate(x_cg_fuel_tanks):
        new_m, new_cg = update_m_and_cg(W_plot_ftb[-1], x_cg_lemac_plot_ftb[-1], m_fuel_tanks[i], x_cg)
        W_plot_ftb.append(new_m)
        x_cg_lemac_plot_ftb.append(new_cg)
        new_m, new_cg = update_m_and_cg(W_plot_btf[-1], x_cg_lemac_plot_btf[-1], m_fuel_tanks_reverse[i], x_cg_fuel_tanks_reverse[i])
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
    t_tail_condition: bool = ac.empennage.t_tail_condition  # T-tail or not
    tail_type: str = ac.empennage.horizontal_tail['type']  # 'fully moving' or 'adjustable' or 'fixed'

    # Flight conditions
    V_cruise = ac.requirements.cruise['cr_speed'] * KTS_TO_MS
    alt_cr = ac.requirements.cruise['cr_altitude'] * FT_TO_M
    V_app = ac.requirements.approach['ap_speed'] * KTS_TO_MS
    alt_LD = ac.requirements.landing['la_altitude'] * FT_TO_M
    LD_temp_shift = ac.requirements.landing['la_temp_shift']

    # Fuselage params
    b_f = ac.fuselage.width
    l_f = ac.fuselage.length
    h_f = ac.fuselage.height

    # Wing params
    x_le = ac.wing.x_le  # wing le distance from nose
    A = ac.wing.aspect_ratio
    b = ac.wing.span
    S = ac.wing.area
    c_r = ac.wing.c_root
    taper = ac.wing.taper_ratio
    S_net = Snet(S, b_f, taper, b, c_r)
    sweep_c_4_deg = ac.wing.sweep
    sweep_c_2_deg = ac.wing.sweep_c_2_deg  # sweep_at_x_c_deg(LE_sweep_deg(sweep_c_4_deg, c_r, b, taper), c_r, b, taper, x_c=0.5)
    sweep_LE_deg = ac.wing.sweep_LE_deg  # LE_sweep_deg(sweep_c_4_deg, c_r, b, taper)
    C_m0_airfoil = ac.wing.cm_c4
    C_L_0 = np.abs(ac.hld_and_ailerons.landing_lift['alpha_zero_lift']) * ac.hld_and_ailerons.landing_lift['CL_alpha']  # CL of flapped wing at alpha = 0
    C_L_LD = ac.hld_and_ailerons.landing_lift['CL_max'] / 1.21
    mac = ac.wing.MAC
    mgc = S/b

    # HLD
    flap_type: str = ac.hld_and_ailerons.flaps['flap_type']  # containing fowler or slotter or none for Cm calculations
    flap_defl_ld_deg = ac.hld_and_ailerons.flaps['ld_deflection']
    ext_flap_chord_ratio_ld = ac.hld_and_ailerons.flaps['cdash_c_ld']  # c’/c = the ratio between the chord of the airfoil with extended flap and the chord in clean configuration
    cf_c = ac.hld_and_ailerons.flaps['cf_c']  # simple flaps: 0.25, highly efficient slotted flaps: 0.35-0.4
    cf_ext_flap_chord_ld = cf_c / ext_flap_chord_ratio_ld  # c_f/c'
    y_start_f = ac.hld_and_ailerons.flaps['y_flap_in']
    y_end_f = ac.hld_and_ailerons.flaps['y_flap_out']
    Delta_Cl_max_ld = D_Cl_max(flap_type, cdash_c=ext_flap_chord_ratio_ld)  # the airfoil lift coefficient increase due to flap extension at landing condition (estimated in the wing design module)
    flap_span_wing_span = np.abs(y_start_f - y_end_f) / b
    Swf = S_wf(y_start_f, y_end_f, taper, c_r, b)

    # Tail params
    height_ht = ac.empennage.vertical_tail['b_v']  # height of ht above wing (perpendicular to chord lines)
    l_h = ac.empennage.horizontal_tail['l_h']
    A_h = ac.empennage.horizontal_tail['aspect_ratio']
    r = 2 * l_h / b
    m_tv = height_ht * 2 / b

    # Nacelle params
    nac_y = b_f / 2 + ac.engine.eng_y_pos_fuselage
    if ac.engine.eng_x_pos == 'le':
        nac_mount_dist = x_pos_le_along_span_from_nose(sweep_LE_deg, nac_y, x_le)  # Distance front of nacelle from nose
    l_fn = x_pos_le_along_span_from_nose(sweep_LE_deg, b_f/2, x_le)
    b_n = ac.engine.nac_diameter  # nacelle diameter (max)
    x_mac_c_4 = x_pos_le_along_span_from_nose(sweep_LE_deg, y_pos_at_chord_length(c_r, taper, mac, b), x_le) + 0.25 * mac
    l_n = x_mac_c_4 - nac_mount_dist
    if any(word in flap_type.lower() for word in ['fowler', 'slotted']):
        fowler_data = np.array([[0.15, 0.2, 0.25, 0.3, 0.35, 0.4],
                                [0.235, 0.222, 0.213, 0.207, 0.202, 0.2]])
        cf_ext_flap_chord_ld = np.clip(cf_ext_flap_chord_ld, fowler_data[0,:].min(), fowler_data[0,:].max())
        mu1 = float(np.interp(cf_ext_flap_chord_ld, fowler_data[0,:], fowler_data[1,:]))
        # print(f' mu1: {mu1}, cf/cdash: {cf_ext_flap_chord_ld}')
    else: 
        mu1 = interpolate_csv_z(flap_defl_ld_deg, cf_ext_flap_chord_ld, "lookups/mu_1.csv")
    mu2 = interpolate_csv_z(taper, flap_span_wing_span, 'lookups/mu_2.csv')
    mu3 = interpolate_csv_z(taper, flap_span_wing_span, 'lookups/mu_3.csv')
    # print(f'mu2: {mu2}, flap_span_wing_span: {flap_span_wing_span}, taper: {taper}')
    # print(f'mu3: {mu3}')
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
    x_ac_f1_cr = -1.8 * b_f * h_f * l_fn / (C_L_alpha_A_less_h_cr) / (S * mac)
    x_ac_f1_ld = -1.8 * b_f * h_f * l_fn / (C_L_alpha_A_less_h_LD) / (S * mac)
    x_ac_f2 = 0.273 * b_f * mgc * (b - b_f) * np.tan(np.deg2rad(sweep_c_4_deg)) / ((1 + taper) * mac**2 * (b + 2.15 * b_f))
    x_ac_n_cr = n_eng * (-4) * b_n**2 * l_n / (S * mac * C_L_alpha_A_less_h_cr)
    x_ac_n_ld = n_eng * (-4) * b_n**2 * l_n / (S * mac * C_L_alpha_A_less_h_LD)
    x_ac_cr = x_ac_w_cr + x_ac_f1_cr + x_ac_f2 + x_ac_n_cr
    x_ac_ld = x_ac_w_ld + x_ac_f1_ld + x_ac_f2 + x_ac_n_ld
    # print(f'x_ac_cr: {x_ac_cr} \n x_ac_n_cr: {x_ac_n_cr} \n x_ac_f2: {x_ac_f2} \n x_ac_f1_cr: {x_ac_f1_cr} \n x_ac_w_cr: {x_ac_w_cr}')

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
    C_m_ac_flap_ld = mu2 * (-mu1 * Delta_Cl_max_ld * ext_flap_chord_ratio_ld - (C_L_LD + Delta_Cl_max_ld * (1 - Swf / S)) / 8 * ext_flap_chord_ratio_ld * (ext_flap_chord_ratio_ld - 1)) + 0.7 * A / (1 + 2 / A) * mu3 * Delta_Cl_max_ld * np.tan(np.deg2rad(sweep_c_4_deg))
    C_m_ac_flap_ld += C_L_A_less_h_ld * -0.15
    C_m_ac_ld = C_m_ac_w + C_m_ac_flap_ld + C_m_ac_fus_ld
    print(f' Cm aerodynamic cente stuff: \n C_m_ac_w: {C_m_ac_w} \n C_m_ac_fus_ld: {C_m_ac_fus_ld} \n C_m_ac_flap_ld: {C_m_ac_flap_ld} \n total: {C_m_ac_ld}')

    # Tail sizes
    Sh_S_cont = (x_cg_lemac_mac - (x_ac_ld - C_m_ac_ld / C_L_A_less_h_ld)) / (C_L_h / C_L_A_less_h_ld * l_h / mac * V_h_V2)
    Sh_S_n_stab = (x_cg_lemac_mac - x_ac_cr) / (C_L_alpha_h_cr / C_L_alpha_A_less_h_cr * (1 - depsilon_dalpha_cr) * l_h / mac * V_h_V2)
    Sh_S_stab = (x_cg_lemac_mac - (x_ac_cr - SM)) / (C_L_alpha_h_cr / C_L_alpha_A_less_h_cr * (1 - depsilon_dalpha_cr) * l_h / mac * V_h_V2)

    print(f'x_cg_lemac_mac: {x_cg_lemac_mac}, \nx_ac_cr: {x_ac_cr}, \nSM: {SM}, \nC_L_alpha_h_cr: {C_L_alpha_h_cr}, \nC_L_alpha_A_less_h_cr: {C_L_alpha_A_less_h_cr}, \ndepsilon_dalpha_cr:{depsilon_dalpha_cr}, \nl_h: {l_h}, \nmac: {mac}, \nV_h_V2: {V_h_V2}')
    print(f'\nx_ac_ld: {x_ac_ld}, \nC_m_ac_ld: {C_m_ac_ld}, \nC_L_A_less_h_ld: {C_L_A_less_h_ld}, \nC_L_h: {C_L_h}, \nC_L_A_less_h_ld: {C_L_A_less_h_ld}')
    # PLotting
    min_x = min(x_cg_lemac_mac)
    max_x = max(x_cg_lemac_mac)
    range_x = np.abs(max_x - min_x)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(min_x - 0.1 * range_x, max_x + 0.1 * range_x)
    # ax.set_ylim(0, 0.7)
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
                                      show_plot: bool = False,
                                      update_ac: bool = False):
    x_cg_lemac_mac_plot = convert_x_cg_from_nose_to_lemac_frac_mac(x_le_w_fus_length_arr * ac.fuselage.length, x_le_w_fus_length_arr * ac.fuselage.length, ac)
    x_cg_lemac_mac_plot = np.arange(-1.5, 1.5, 0.01)
    print(f'x_Cg_lemac_mac_plot: {x_cg_lemac_mac_plot}')
    fwd_cg = np.zeros_like(x_le_w_fus_length_arr)
    aft_cg = np.zeros_like(x_le_w_fus_length_arr)
    for i, x_le_w_fus_length in enumerate(x_le_w_fus_length_arr):
        fwd_cg[i], aft_cg[i], ac1 = loading_diagram(x_le_w_fus_length * ac.fuselage.length, ac, update_ac_cgs=False)
    x_le_w_l_fus = x_le_w_fus_length_arr
    Sh_S_cont, Sh_S_n_stab, Sh_S_stab = scissor_plot(ac, x_cg_lemac_mac_plot, output_filepath='outputs/Initial_scissor_plot.png', show_plot=True)
    print(f'fwd_cg: {fwd_cg}, aft cg: {aft_cg}')
    # Plotting

    # Shared x-axis
    x = np.linspace(0, 1, 100)

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

    print(f'x1={x_cg_lemac_mac_plot}, x2={fwd_cg}, y1={Sh_S_cont}, y2={x_le_w_l_fus}')
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
        fwd_cg = float(np.interp(Sh_S, Sh_S_cont, x_cg_lemac_mac_plot))
        aft_cg = float(np.interp(Sh_S, Sh_S_stab, x_cg_lemac_mac_plot))
        if update_ac:
            ac.weights.x_cg_aft = aft_cg
            ac.weights.x_cg_fwd = fwd_cg
            ac.empennage.horizontal_tail['area div S'] = Sh_S
            ac.wing.x_le = x_le_from_x_lemac(wing_pos * ac.wing.MAC, ac.wing.y_MAC, ac.wing.sweep_LE_deg)
        return Sh_S, wing_pos, aft_cg, fwd_cg, x_le_from_x_lemac(wing_pos * ac.wing.MAC, ac.wing.y_MAC, ac.wing.sweep_LE_deg)
    else: 
        return None


if __name__ == "__main__":
    ac = Aircraft('Boosted_turboprop_tricycle',
                loader.load('concepts/reqs_turb.yaml', Requirements),
                loader.load('yamls/mission.yaml', Mission),
                loader.load('yamls/weights.yaml', Weights),
                loader.load('concepts/wing_electra.yaml', Wing),
                loader.load('concepts/fus_tri.yaml', Fuselage),
                loader.load('concepts/engine_tprop_b.yaml', Engine),
                loader.load('concepts/tricycle_empennage.yaml', Empennage),
                loader.load('yamls/HLD_and_ailerons.yaml', HLD_and_AIL),
                loader.load('concepts/tricycle_gear.yaml', Landing_Gear))
    scissor_plot(ac, x_cg_lemac_mac=np.arange(-1.5, 1.5, 0.01), SM=0.05, output_filepath='outputs/test_scissor_plot.png', show_plot=True)