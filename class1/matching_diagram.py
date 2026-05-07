"""
Matching diagram and design point selection utilities.

Uses Roskam correlations to generate matching plot and select design point.
Uses prelim_drag.py and its inputs, maximum lift coefficients during landing, cruise, ..., landing and take-off field
length requirements, stall speed, CS-23 regulations etc.

Type of plot (P/W or T/W on y-axis), equations and selected optimal design point are selected
based on aircraft engine configuration.
"""

import sys
import os
from classes.aircraft_2 import Aircraft, loader
from class1.prelim_drag import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from classes.isa import Atmosphere
from lookups.consts import *
from typing import Any
from numpy import dtype, ndarray
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))


def stall_speed_matching(ac: Aircraft,  # Change units
                         W_P: np.ndarray = np.arange(5,10000)
                         ):
    V_s0 = ac.requirements.general['stall_speed'] * KTS_TO_MS

    landing_altitude = ac.requirements.landing['la_altitude'] * FT_TO_M
    landing_temperature_shift = ac.requirements.landing['la_temp_shift']
    atmos_model = Atmosphere(landing_altitude, landing_temperature_shift)
    rho = atmos_model.density

    CL_max_landing = ac.requirements.landing['as_CL_max_la']

    W_S = (CL_max_landing * 0.5 * rho * V_s0 ** 2) * np.ones_like(W_P)

    return W_P, W_S


def landing_dist_matching(ac: Aircraft,
                          W_P: np.ndarray = np.arange(5,10000)
                          ):
    W_L = ac.weights.m_takeoff * ac.requirements.landing['la_mass_frac']
    C_L_max_landing = ac.requirements.landing['as_CL_max_la']

    balance_fl_condition = ac.requirements.landing['dist_balance_fl_condition']

    landing_altitude = ac.requirements.landing['la_altitude'] * FT_TO_M
    landing_temperature_shift = ac.requirements.landing['la_temp_shift']
    atmos_model = Atmosphere(landing_altitude, landing_temperature_shift)
    rho = atmos_model.density

    if balance_fl_condition:
        # Roskam 3.14
        S_L = ac.requirements.landing['la_distance'] / FT_TO_M
        V_S_L = np.sqrt(S_L / 0.5136) * KTS_TO_MS

    else:
        # Roskam 3.12
        S_LG = ac.requirements.landing['la_distance'] / FT_TO_M
        V_S_L = np.sqrt(S_LG / 0.265) * KTS_TO_MS

    W_S_LD = (V_S_L ** 2 / (1 / 2 * rho * C_L_max_landing)) * np.ones_like(W_P)
    W_S = W_S_LD * ac.weights.m_takeoff / W_L

    return W_P, W_S


def cruise_speed_matching(ac: Aircraft,
                          type_to_use : str = "Single Engine Propeller Driven",
                          friction_source : str = 'lookups/skin_fric.csv',
                          s_wet_source : str = 'lookups/s_wets.csv',
                          W_S: np.ndarray = np.arange(5,10000)
                          ):

    V_cr = ac.requirements.cruise['cr_speed'] * KTS_TO_MS

    cruise_altitude = ac.requirements.cruise['cr_altitude'] * FT_TO_M
    atmos_model = Atmosphere(cruise_altitude, 0)
    rho = atmos_model.density
    sigma = atmos_model.density_ratio

    cruise_mass_frac = ac.requirements.cruise['cr_mass_frac']
    CD0 = prelim_drag.cd0(ac, type_to_use, friction_source, s_wet_source)
    A = ac.wing.aspect_ratio
    _, e = prelim_drag.k(ac)

    eta_p = ac.engine.eta_3

    # NOTE: first eqn for piston, second for turboprop, third for electric, check ADSEE I book for details (dep on propeller type)
    # engine.Phi: power split parameter
    alpha_p_piston = 1.132 * sigma - 0.132
    alpha_p_turboprop = sigma ** 0.75
    alpha_p_electric = 1

    alpha_p = alpha_p_electric * ac.engine.Phi + alpha_p_turboprop * (1 - ac.engine.Phi)

    W_P = eta_p * alpha_p / cruise_mass_frac * (CD0 * 0.5 * rho * V_cr ** 3 / (cruise_mass_frac * W_S) + cruise_mass_frac * W_S / (np.pi * A * e * 0.5 * rho * V_cr)) ** (-1)

    return W_P, W_S


def takeoff_dist_matching(ac: Aircraft,  # Change units
                          W_S: np.ndarray = np.arange(5,10000)
                          ):

    W_S = W_S*PA_TO_LBSpFT2
    to = ac.requirements.take_off
    balance_fl_condition = to['dist_balance_fl']  # True or False dep on yaml

    alt = to['to_altitude'] * FT_TO_M
    temp_shift = to['to_temperature']
    sigma = Atmosphere(alt, temp_shift).density_ratio
    C_Lmax_TO = to['as_CL_max_to']

    if balance_fl_condition:
        S_TO = to['to_distance'] / FT_TO_M
        # Eqn 3.6 Roskam
        x = (-8.134 + np.sqrt(8.134 ** 2 + S_TO * 4 * 0.0149)) / (2 * 0.0149)
    else:
        # Roskam eq 3.4
        S_TOG = to['to_distance'] / FT_TO_M
        x = (-4.9 + np.sqrt(4.9 ** 2 + 4 * 0.009 * S_TOG)) / (2 * 0.009)

    W_P = (x * sigma*C_Lmax_TO / W_S) * LBSpHP_TO_NpW

    return W_P, W_S


def C_L3_2_C_D_max(C_D0,
                   A,
                   e
                   ):
    return 1.345 * (A * e) ** (3 / 4) / (C_D0 ** (1 / 4))


def delta_e(ac: Aircraft, flap_deflection):
    if ac.requirements.climb['eng_placement'] == 'fuselage':
        delta = 0.0046 * flap_deflection
        return delta
    elif ac.requirements.climb['eng_placement'] == 'wing':
        delta = 0.0026 * flap_deflection
        return delta
    else:
        return None


def delta_cd0(ac: Aircraft, flap_deflection):
    delta_from_flaps = 0.0013 * flap_deflection
    delta_from_lg = 0.0250

    return delta_from_flaps + delta_from_lg


def W_P_for_ROC(RC,
                eta_p,
                W_S,
                C_L_C_D_param,
                sigma):
    W_P = eta_p / (RC / 33000 + np.sqrt(W_S) / (19 * C_L_C_D_param * np.sqrt(sigma)))  # lbs/hp
    W_P = W_P * LBSpHP_TO_NpW

    return W_P


def W_P_for_CGR(CGR,
                L_D,
                C_L_climb,
                eta_p,
                sigma,
                W_S):
    CGRP = (CGR + 1 / L_D) / (np.sqrt(C_L_climb))
    W_P = 18.97 * eta_p * np.sqrt(sigma) / (CGRP * np.sqrt(W_S))
    W_P = W_P * LBSpHP_TO_NpW

    return W_P


def all_engine_operative(ac : Aircraft,
                         type_to_use : str = "Single Engine Propeller Driven",
                         friction_source : str = 'lookups/skin_fric.csv',
                         s_wet_source : str = 'lookups/s_wets.csv',
                         W_S: np.ndarray = np.arange(5,10000)
                         ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:

    eta_p = ac.engine.eta_3

    take_off_altitude = ac.requirements.take_off['to_altitude'] * FT_TO_M
    take_off_temperature_shift = ac.requirements.take_off['to_temp_shift']
    atmos_model = Atmosphere(take_off_altitude, take_off_temperature_shift)
    sigma = atmos_model.density_ratio

    CD0 = prelim_drag.cd0(ac, type_to_use, friction_source, s_wet_source)
    A = ac.wing.aspect_ratio
    k, e = prelim_drag.k(ac)

    take_off_flap_deflection = ac.requirements.climb['take_off_flap_deflection']
    e = e + delta_e(ac, take_off_flap_deflection)
    CD0 = CD0 + delta_cd0(ac, take_off_flap_deflection)
    C_L_C_D_param = C_L3_2_C_D_max(CD0, A, e)

    C_L_max_take_off = ac.requirements.take_off['as_CL_max_to']
    C_L_climb = C_L_max_take_off - 0.2
    L_D = C_L_climb / (CD0 + k * C_L_climb ** 2)

    RC = 300  # ft/min
    W_S *= PA_TO_LBSpFT2
    W_P_AEO_turbine_and_non_turbine_ROC = W_P_for_ROC(RC, eta_p, W_S, C_L_C_D_param, sigma)

    CGR = 1 / 12
    W_P_AEO_turbine_and_non_turbine_CGR = W_P_for_CGR(CGR, L_D, C_L_climb, eta_p, sigma, W_S)

    # Additional conditions for turbined aircraft
    W_P_AEO_turbine_additional_condition_CGR = W_S * 0
    if ac.requirements.climb['turbine_condition']:
        take_off_altitude_turbine = ac.requirements.climb_gradient['to_altitude_turbine'] * FT_TO_M
        take_off_temperature_shift_turbine = ac.requirements.climb_gradient['to_temperature_shift_turbine']
        atmos_model = Atmosphere(take_off_altitude_turbine, take_off_temperature_shift_turbine)
        sigma = atmos_model.density_ratio

        CGR = 0.04
        W_P_AEO_turbine_additional_condition_CGR = W_P_for_CGR(CGR, L_D, C_L_climb, eta_p, sigma, W_S)

    return W_P_AEO_turbine_and_non_turbine_ROC, W_P_AEO_turbine_and_non_turbine_CGR, W_P_AEO_turbine_additional_condition_CGR, W_S


def one_engine_inoperative(ac : Aircraft,
                           type_to_use : str = "Single Engine Propeller Driven",
                           friction_source : str = 'lookups/skin_fric.csv',
                           s_wet_source : str = 'lookups/s_wets.csv',
                           W_S: np.ndarray = np.arange(5,10000)
                           ):

    W_P_OEI_turbine_condition_1_ROC = W_S * 0
    W_P_OEI_turbine_condition_1_CGR = W_S * 0
    W_P_OEI_turbine_condition_2_ROC = W_S * 0
    W_P_OEI_turbine_condition_2_CGR = W_S * 0

    if ac.requirements.climb['turbine_condition']:
        eta_p = ac.engine.eta_3  # Propeller efficiency

        # Check
        take_off_altitude_turbine = ac.requirements.climb['to_altitude_turbine'] * FT_TO_M
        atmos_model = Atmosphere(take_off_altitude_turbine, 0)
        alt_turbine = ac.requirements.climb_gradient['to_altitude_turbine'] * FT_TO_M
        temp_turbine = atmos_model.temp_isa
        sigma = Atmosphere(alt_turbine, temp_turbine).density_ratio

        CD0 = prelim_drag.cd0(type_to_use, friction_source, s_wet_source)
        A = ac.wing.aspect_ratio
        k, e = prelim_drag.k(ac)
        e = e + delta_e(ac)
        CD0 = CD0 + delta_cd0(ac)
        C_L_C_D_param = C_L3_2_C_D_max(CD0, A, e)

        C_L_max_TO = ac.requirements.take_off['as_CL_max']
        C_L_climb = C_L_max_TO - 0.2
        L_D = C_L_climb / (CD0 + k * C_L_climb ** 2)

        RC = 0.027 * ac.requirements.general['stall_speed'] ** 2  # ft/min
        W_P_OEI_turbine_condition_1_ROC = W_P_for_ROC(RC, eta_p, W_S, C_L_C_D_param, sigma)

        CGR = 0.012
        W_P_OEI_turbine_condition_1_CGR = W_P_for_CGR(CGR, L_D, C_L_climb, eta_p, sigma, W_S)

        alt_turbine = ac.requirements.climb_gradient['to_altitude_turbine'] * FT_TO_M
        temp_shift_turbine = ac.requirements.climb_gradient['to_temperature_shift_turbine']
        sigma = Atmosphere(alt_turbine, temp_shift_turbine).density_ratio

        RC = 0.014 * ac.requirements.general['stall_speed'] ** 2  # ft/min
        W_P_OEI_turbine_condition_2_ROC = W_P_for_ROC(RC, eta_p, W_S, C_L_C_D_param, sigma)

        CGR = 0.006
        W_P_OEI_turbine_condition_2_CGR = W_P_for_CGR(CGR, L_D, C_L_climb, eta_p, sigma, W_S)

    return W_P_OEI_turbine_condition_1_ROC, W_P_OEI_turbine_condition_1_CGR, W_P_OEI_turbine_condition_2_ROC, W_P_OEI_turbine_condition_2_CGR, W_S


def balked_landing(ac : Aircraft,
        type_to_use : str = "Single Engine Propeller Driven",
        friction_source : str = 'lookups/skin_fric.csv',
        s_wet_source : str = 'lookups/s_wets.csv',
        W_S: np.ndarray = np.arange(5,10000)):

    eta_p = ac.engine.eta_3  # Propeller efficiency

    alt = ac.requirements.take_off['to_altitude'] * FT_TO_M
    temp_shift = ac.requirements.take_off['to_temperature_shift']
    sigma = Atmosphere(alt, temp_shift).density_ratio

    CD0 = prelim_drag.cd0(type_to_use, friction_source, s_wet_source)
    A = ac.wing.aspect_ratio
    k, e = prelim_drag.k(ac)
    e = e + delta_e(ac)
    CD0 = CD0 + delta_cd0(ac)
    C_L_C_D_param = C_L3_2_C_D_max(CD0, A, e)

    C_L_max_landing = ac.requirements.landing['as_CL_max']
    L_D = C_L_max_landing / (CD0 + k * C_L_max_landing ** 2)

    CGR = 1 / 30
    W_P_BL_turbine_and_non_turbine_CGR = W_P_for_CGR(CGR, L_D, C_L_max_landing, eta_p, sigma, W_S)

    W_P_BL_turbine_additional_condition_ROC = W_S * 0
    if ac.requirements.climb_gradient['engine type'] == 'turbine':
        alt_turbine = ac.requirements.climb_gradient['to_altitude_turbine'] * FT_TO_M
        temp_shift_turbine = ac.requirements.climb_gradient['to_temperature_shift_turbine']
        sigma = Atmosphere(alt_turbine, temp_shift_turbine).density_ratio

        RC = 0  # ft/min
        W_P_BL_turbine_additional_condition_ROC = W_P_for_ROC(RC, eta_p, W_S, C_L_C_D_param, sigma)

    return W_P_BL_turbine_and_non_turbine_CGR, W_P_BL_turbine_additional_condition_ROC, W_S


# NOTE: add this function to logbook because Claude helped
def find_design_point(datasets, max_wingloading,
                      ws_margin_frac=0.05, wp_margin_frac=0.05):
    """
    Automatically selects the optimal (most upper-right) design point
    from a matching diagram, respecting all constraints with margin.

    Parameters
    ----------
    datasets       : list of dicts with keys 'x', 'y', 'label'
                     x = W/S array, y = W/P or T/W array
    max_wingloading: float — hard upper limit on W/S (e.g. structural limit)
    ws_margin_frac : float — fractional inward margin on W/S  (default 5%)
    wp_margin_frac : float — fractional inward margin on W/P  (default 5%)

    Returns
    -------
    dict with keys: W_S, W_P, limiting_ws_constraint, limiting_wp_constraint
    """

    # --- 1. Identify vertical-line constraints (stall, landing) -----------
    # These are datasets where x is a scalar or near-constant (vertical line)
    vertical_labels = ["Stall speed", "Landing field length"]
    vertical_limits = {"max wing loading": max_wingloading}

    for ds in datasets:
        lbl = ds["label"].lower()
        if any(v in lbl for v in vertical_labels):
            xs = np.atleast_1d(ds["x"])
            vertical_limits[ds["label"]] = float(np.min(xs))

    # Optimal W/S = min of all vertical constraints, minus margin
    ws_limit = min(vertical_limits.values())
    ws_opt = ws_limit * (1.0 - ws_margin_frac)

    limiting_ws = min(vertical_limits, key=vertical_limits.get)

    # --- 2. Evaluate all curves at ws_opt to find upper W/P bound ----------
    # Curves where feasible region is BELOW: optimal W/P = min of all curves
    curve_labels = ["Take-off field length", "Cruise speed", "AEO RoC", "AEO Climb gradient", "AEO Climb gradient (turbine)", "Balked landing", "Balked landing (turbine)", "OEI RoC/Climb gradient I (turbine)", "OEI RoC/Climb gradient II (turbine)"]
    curve_values = {}

    for ds in datasets:
        lbl = ds["label"].lower()
        if any(c in lbl for c in curve_labels):
            xs = np.asarray(ds["x"], dtype=float)
            ys = np.asarray(ds["y"], dtype=float)
            if xs.size < 2:
                continue
            # Interpolate curve value at ws_opt (clamp to data range)
            ws_clamped = np.clip(ws_opt, xs.min(), xs.max())
            wp_at_ws = float(np.interp(ws_clamped, xs, ys))
            curve_values[ds["label"]] = wp_at_ws

    if not curve_values:
        raise ValueError("No curve datasets found to constrain W/P or T/W.")

    # Feasible = below all curves → pick the minimum as the upper bound
    wp_limit = min(curve_values.values())
    wp_opt = wp_limit * (1.0 - wp_margin_frac)

    limiting_wp = min(curve_values, key=curve_values.get)

    return {
        "W_S": ws_opt,
        "W_P": wp_opt,
        "limiting_ws_constraint": limiting_ws,
        "limiting_wp_constraint": limiting_wp,
        "all_ws_limits": vertical_limits,
        "all_wp_at_design": curve_values,
    }


def plot_matching_and_select_design_point(ac : Aircraft,  # Change units
        type_to_use : str = "Single Engine Propeller Driven",
        friction_source : str = 'lookups/skin_fric.csv',
        s_wet_source : str = 'lookups/s_wets.csv',
        W_S_plot: np.ndarray = np.arange(5,10000),
        W_P_plot: np.ndarray = np.arange(5,10000),
        output_filepath: str = 'outputs/Matching_Diagram.png') -> list:

    # Compute plots
    stall_W_P, stall_W_S = stall_speed_matching(ac, W_P_plot)
    to_W_P, to_W_S = takeoff_dist_matching(ac, W_S_plot)
    ld_W_P, ld_W_S = landing_dist_matching(ac, W_P_plot)
    cr_W_P, cr_W_S = cruise_speed_matching(ac, type_to_use, friction_source, s_wet_source, W_S_plot)

    AEO1_W_P, AEO2_W_P, AEO3_turb_W_P, AEO_W_S = all_engine_operative(ac, type_to_use, friction_source, s_wet_source, W_S_plot)
    OEI1a_W_P, OEI1b_W_P, OEI2a_W_P, OEI2b_W_P, OEI_W_S = one_engine_inoperative(ac, type_to_use, friction_source, s_wet_source, W_S_plot)
    BL_W_P, BL_turb_W_P, BL_W_S = balked_landing(ac, type_to_use, friction_source, s_wet_source, W_S_plot) 

    OEI1_W_P = np.min(OEI1a_W_P, OEI1b_W_P)
    OEI2_W_P = np.min(OEI2a_W_P, OEI2b_W_P)

    # Start actual plotting stuff
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(min(W_S_plot), max(W_S_plot))
    ax.set_ylim(min(W_P_plot), max(W_P_plot))

    datasets = [
        {"x": stall_W_S, "y": stall_W_P, "label": "Stall speed"},
        {"x": to_W_S, "y": to_W_P, "label": "Take-off field length"},
        {"x": ld_W_S, "y": ld_W_P, "label": "Landing field length"},
        {"x": cr_W_S, "y":cr_W_P, "label": "Cruise speed"},
        {"x": AEO_W_S, "y": AEO1_W_P, "label": "AEO RoC"},
        {"x": AEO_W_S, "y": AEO2_W_P, "label": "AEO Climb gradient"},
        {"x": BL_W_S, "y": BL_W_P, "label": "Balked landing"}
    ]

    if ac.requirements.climb['turbine_condition']:
        datasets.append({"x": AEO_W_S, "y": AEO3_turb_W_P, "label": "AEO Climb gradient (turbine)"})
        datasets.append({"x": BL_turb_W_P, "y":BL_W_S, "label": "Balked landing (turbine)"})

    if ac.requirements.climb['turbine_condition'] and ac.requirements.climb_gradient['n_eng']:
        {"x": OEI_W_S, "y": OEI1_W_P, "label": "OEI RoC/Climb gradient I (turbine)"},
        {"x": OEI_W_S, "y": OEI2_W_P, "label": "OEI RoC/Climb gradient II (turbine)"},

    colors = cm.tab10(np.linspace(0, 1, len(datasets)))

    for ds, color in zip(datasets, colors):
        ax.plot(ds["x"], ds["y"], color=color, label=ds["label"], linewidth=2)

    result = find_design_point(datasets, max_wingloading=ac.requirements.general["max_wing_loading"]*g,
                                ws_margin_frac=0.05,
                                wp_margin_frac=0.05)

    print(f"Design point:  W/S = {result['W_S']:.1f}  |  W/P = {result['W_P']:.4f}")
    print(f"Limited in W/S by: {result['limiting_ws_constraint']}")
    print(f"Limited in W/P by: {result['limiting_wp_constraint']}")

    # Plotting design point
    ax.scatter(result["W_S"], result["W_P"],
            marker="*", s=250, color="red", zorder=10,
            label=f"Design point ({result['W_S']:.0f}, {result['W_P']:.4f})")

    ax.annotate(
        f"  \t DESIGN POINT: \n"
        f"  W/S = {result['W_S']:.1f}\n"
        f"  W/P = {result['W_P']:.4f}\n"
        f"  [{result['limiting_ws_constraint']}]\n"
        f"  [{result['limiting_wp_constraint']}]",
        xy=(result["W_S"], result["W_P"]),
        xytext=(result["W_S"] * 0.75, result["W_P"] * 1.15),
        fontsize=8,
        arrowprops=dict(arrowstyle="->", color="red"),
        color="red",
    )

    # Labels and shit
    ax.set_xlabel(f"Wing loading $\\frac{{W}}{{S}}_{{TO}}$ [N/m$^2$]")
    ax.set_ylabel(f"Power loading $\\frac{{W}}{{P}}_{{TO}}$ [N/W]")
    ax.set_title("Matching diagram")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.show()
    plt.savefig(output_filepath, dpi=300)

    data = {
        "W/P": result['W_P'],
        "W/S": result['W_S'],
        "limiting_ws_constraint": result['limiting_ws_constraint'],
        "limiting_wp_constraint": result['limiting_wp_constraint'],
    }

    return data


if __name__ == '__main__':
    file_path = "../yamls/aircraft.yaml"
    target_class = Aircraft
    ac = loader.load(file_path, target_class)

    from pathlib import Path

    BASE_DIR = Path(__file__).resolve().parent

    friction_source = BASE_DIR / "../lookups/skin_fric.csv"
    s_wet_source = BASE_DIR / "../lookups/s_wets.csv"

    type_to_use: str = "Single Engine Propeller Driven"
    # friction_source: str = 'lookups/skin_fric.csv'
    # s_wet_source: str = 'lookups/s_wets.csv'
    W_S: np.ndarray = np.arange(5, 10000)

    V_cr = ac.requirements.cruise['cr_speed'] * KTS_TO_MS

    cruise_altitude = ac.requirements.cruise['cr_altitude'] * FT_TO_M
    atmos_model = Atmosphere(cruise_altitude, 0)
    rho = atmos_model.density
    sigma = atmos_model.density_ratio

    cruise_mass_frac = ac.requirements.cruise['cr_mass_frac']
    CD0 = cd0(ac, type_to_use, friction_source, s_wet_source)
    A = ac.wing.aspect_ratio
    _, e = k(ac)

    eta_p = 0.7

    # NOTE: first eqn for piston, second for turboprop, third for electric, check ADSEE I book for details (dep on propeller type)
    # engine.Phi: power split parameter
    alpha_p_piston = 1.132 * sigma - 0.132
    alpha_p_turboprop = sigma ** 0.75
    alpha_p_electric = 1
    ac.engine.Phi = 0.4

    alpha_p = alpha_p_electric * ac.engine.Phi + alpha_p_turboprop * (1 - ac.engine.Phi)

    W_P = eta_p * alpha_p / cruise_mass_frac * (CD0 * 0.5 * rho * V_cr ** 3 / (cruise_mass_frac * W_S) + cruise_mass_frac * W_S / (np.pi * A * e * 0.5 * rho * V_cr)) ** (-1)


    # print(target.requirements.take_off['to_distance'])
    # climb_rate_AEO_matching_23_65(ac = target, type_to_use = target.requirements.general['standard_type'])


# def alpha_thrust(alt_m, B, T_K, V_mps, theta_break=1.07):
#     p = Atmosphere(alt_m, delta_T=20) # ISA pressure at alt
#     M = V_mps/np.sqrt(1.4*287*T_K)
#     delta_t = p/101325*(1+0.2*M**2)**(1.4/0.4)
#     theta_t = T_K/288.15*(1+0.2*M**2)
#     if 0<B<5:
#         if theta_t <= theta_break:
#             return delta_t
#         elif theta_t>theta_break:
#             return delta_t*(1-2.1*(theta_t-theta_break)/theta_t)
#     elif 5<B<15:
#         if theta_t <= theta_break:
#             return delta_t*(1-(0.43+0.014*B)*np.sqrt(M))
#         elif theta_t>theta_break:
#             return delta_t*(1-(0.43+0.014*B)*np.sqrt(M) -3*(theta_t-theta_break)/(1.5+M))