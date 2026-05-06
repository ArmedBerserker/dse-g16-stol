"""
Matching diagram and design point selection utilities.

Uses Roskam correlations to generate matching plot and select design point. 
Uses preliminary_drag.py and its inputs, maximum lift coefficients during landing, 
cruise, ..., landing and take-off field length requirements, stall speed, CS-23 regulations etc.
Type of plot (P/W or T/W on y-axis), equations and selected optimal design point are selected
based on aircraft engine configuration.
"""

from classes.aircraft_2 import Aircraft
from lookups.consts import *
from class1.prelim_drag import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from classes.isa import Atmosphere

# density = Atmosphere(3000, 0).density
# print(density)

''' TO DO:
    - Change matching plotter to include climb requirements lines properly'''

''' Key points to pay attention to:
    - Finish RoC and Climb gradient (after trade-off)
    - Are conversions done right for climb gradient and climb rate
    - Are aircraft engine types checked for right?
    - Check OEI requirement for non-turbines
    - Add hard limit W/S and design point selection and assessing AR/CLmax impact'''

def stall_speed_matching(ac : Aircraft,  # Change units
        W_P_or_T_W: np.ndarray = np.arange(0,10000,5)) -> np.ndarray:
    ''' Returns list of arrays (W/P or T/W and W_S)'''
    V_s = (ac.requirements.general.stall_speed)*KTS_TO_MS
    rho = Atmosphere(height=ac.requirements.landing['la_altitude']*FT_TO_M, delta_T=ac.requirements.landing['la_temperature_shift']) 
    CL_max_LD = ac.requirements.cruise['as_CL_max']
    W_S = (CL_max_LD*0.5*rho*V_s**2)*np.ones_like(W_P_or_T_W)
    return (W_P_or_T_W, W_S)

def takeoff_dist_matching(ac : Aircraft,  # Change units
        W_S: np.ndarray = np.arange(0,10000,5)) -> float:
    ''' Returns list of arrays (W/P or T/W and W_S)'''
    W_S = W_S*PA_TO_LBSpFT2

    balance_fl = ac.requirements.take_off['dist_balance_fl']  # True or False dep on yaml

    alt = ac.requirements.take_off.to_altitude*FT_TO_M
    temp_shift = ac.requirements.take_off.to_temperature
    sigma = Atmosphere(alt, temp_shift).density_ratio 
    C_Lmax_TO = ac.requirements.take_off['as_CL_max']

    if balance_fl == True:
        S_TO = ac.requirements.take_off.to_distance/FT_TO_M
        ''' Eqn 3.6 Roskam'''
        x = (-8.134+np.sqrt(8.134**2+S_TO*4*0.0149))
    else:
        S_TOG = ac.requirements.take_off.to_distance/FT_TO_M
        x = (-4.9+np.sqrt(4.9**2+4*0.009*S_TOG))
        # S_TO = 1.66*S_TOG

    W_P = (x*sigma*C_Lmax_TO/W_S)*LBSpHP_TO_NpW
    return (W_P, W_S)


def landing_dist_matching(ac : Aircraft,
        W_P_or_T_W: np.ndarray = np.arange(0,10000,5)) -> float:
    ''' Returns list of arrays (W/P or T/W and W_S)'''
    # W_S = np.ones_like(W_P_or_T_W)

    W_L = ac.weights.m_takeoff*ac.requirements.landing.mass_frac 
    C_L_max_LD = ac.requirements.landing['as_CL_max']

    balance_fl = ac.requirements.take_off['dist_balance_fl']  # True or False dep on yaml

    alt = ac.requirements.landing.la_altitude*FT_TO_M
    temp_shift = ac.requirements.landing.la_temperature
    # sigma =  # Density ratio
    rho = Atmosphere(alt, temp_shift)

    if balance_fl == True:
        S_L = ac.requirements.landing.la_distance/FT_TO_M
        # S_LG = S_L/1.938
        V_S_L = np.sqrt(S_L/0.5136)*KTS_TO_MS

    else: 
        S_LG = ac.requirements.landing.la_distance/FT_TO_M
        V_S_L = np.sqrt(S_LG/0.265)*KTS_TO_MS

    W_S_LD = (V_S_L**2/(1/2*rho*C_L_max_LD))*np.ones_like(W_P_or_T_W)
    W_S = W_S_LD*ac.weights.m_takeoff/W_L

    return (W_P_or_T_W, W_S)

def alpha_thrust(alt_m, B, T_K, V_mps, theta_break=1.07):
    p = Atmosphere(alt_m, delta_T=20) # ISA pressure at alt
    M = V_mps/np.sqrt(1.4*287*T_K)
    delta_t = p/101325*(1+0.2*M**2)**(1.4/0.4)
    theta_t = T_K/288.15*(1+0.2*M**2)
    if 0<B<5:
        if theta_t <= theta_break:
            return delta_t
        elif theta_t>theta_break:
            return delta_t*(1-2.1*(theta_t-theta_break)/theta_t)
    elif 5<B<15:
        if theta_t <= theta_break:
            return delta_t*(1-(0.43+0.014*B)*np.sqrt(M))
        elif theta_t>theta_break:
            return delta_t*(1-(0.43+0.014*B)*np.sqrt(M) -3*(theta_t-theta_break)/(1.5+M))
        

def cruise_speed_matching(ac : Aircraft,
        type_to_use : str = "Single Engine Propeller Driven",
        friction_source : str = 'lookups/skin_fric.csv',
        s_wet_source : str = 'lookups/s_wets.csv',
        W_S: np.ndarray = np.arange(0,10000,5)) -> float:
    ''' Returns list of arrays (W/P or T/W and W_S)'''


    V_cr = ac.requirements.cruise.cr_speed*KTS_TO_MS
    alt = ac.requirements.cruise.cr_altitude*FT_TO_M
    mass_frac = ac.requirements.cruise.cr_mass_frac
    rho = Atmosphere(alt, delta_T=0)
    CD0 = prelim_drag.cd0(type_to_use, friction_source, s_wet_source)
    A = ac.wing.aspect_ratio
    e = prelim_drag.k(type_to_use, friction_source, s_wet_source)[1]

    sigma = Atmosphere(alt, delta_T=0)
    eta_p = ac.engine.eta_prop
    alpha_p = 1.132*sigma-0.132
    alpha_p = sigma**0.75
    alpha_p = 1
    alpha_p = # NOTE: first eqn for piston, second for turboprop, third for electric, check ADSEE I book for details (dep on propeller type)
    W_P = eta_p*alpha_p/mass_frac*(CD0*0.5*rho*V_cr**3/(mass_frac*W_S)+mass_frac*W_S/(np.pi*A*e*0.5*rho*V_cr))
    return (W_P, W_S)


def C_L3_2_C_D_max(C_D0, A, e):
    return 1.345*(A*e)**(3/4)/(C_D0**(1/4))

def climb_rate_AEO_matching_23_65(ac : Aircraft,
        type_to_use : str = "Single Engine Propeller Driven",
        friction_source : str = 'lookups/skin_fric.csv',
        s_wet_source : str = 'lookups/s_wets.csv',
        W_S: np.ndarray = np.arange(0,10000,5)) -> float:
    ''' Returns list of arrays (W/P or T/W and W_S)
    
    min climb rate 300 fpm at SL, 1:12 angle:
    - Not more than max cont. power
    - Gear not retracted
    - Flaps in TO config
    - Cowl flaps? 
    
    for turbines: 
    - climb grad at least 4% at pressure alt 5000ft and 81 deg F under same circumstances'''

    CD0 = prelim_drag.cd0(type_to_use, friction_source, s_wet_source)

    eta_p =  # Propeller efficiency
    sigma = 
    A = 
    e =  # NOTE: Add changes to e and CD0 because of gear and flaps
    C_L_C_D_param = C_L3_2_C_D_max(CD0, A, e)
    RC =  300 # ft/min
    W_P = eta_p/(RC/33000+np.sqrt(W_S)/(19*C_L_C_D_param*np.sqrt(sigma)))  # hp/lbs
    return (W_P, W_S)

def climb_rate_OEI_multiengine_matching_turbine():
    ''' Configuration:
    - OEI (mosti critical one)
    - Prop in minimum drag pos
    - Remaining engines at max continuous poower
    - Gear retracted (if possible)
    - Flaps in most favorable position
    - Cowl flaps?'''
    alt = 5000  # ft

    RC = 


def climb_angle_AEO_matching_23_65(ac : Aircraft,
        type_to_use : str = "Single Engine Propeller Driven",
        friction_source : str = 'lookups/skin_fric.csv',
        s_wet_source : str = 'lookups/s_wets.csv',
        W_S: np.ndarray = np.arange(0,10000,5)) -> float:
    ''' Returns list of arrays (W/P or T/W and W_S)'''
    CGR = 1/12
    eta_p = 
    sigma = 
    C_L_max_TO = 
    C_L_climb = C_L_max_TO - 0.2
    CD0 = prelim_drag.cd0(type_to_use, friction_source, s_wet_source)
    k = prelim_drag.k  # NOTE: add to CD0 bc config and call e + change due to config
    L_D = C_L_climb/(CD0+k*C_L_climb**2)

    CGRP = (CGR+1/L_D)/(np.sqrt(C_L_climb))
    W_P = 18.97*eta_p*np.sqrt(sigma)/CGRP/np.sqrt(W_S)
    return (W_P, W_S)


def climb_angle_AEO_matching_23_77(ac : Aircraft,
        type_to_use : str = "Single Engine Propeller Driven",
        friction_source : str = 'lookups/skin_fric.csv',
        s_wet_source : str = 'lookups/s_wets.csv',
        W_S: np.ndarray = np.arange(0,10000,5)) -> float:
    ''' Returns list of arrays (W/P or T/W and W_S)'''

    ''' Climb angle at least 1:30:
    - TO power all engines
    - Gear down
    - Flaps in landing position'''
    CGR = 1/30
    eta_p = 
    sigma = 
    C_L_max_LD = 
    C_L_climb = C_L_max_LD - 0.2
    CD0 = prelim_drag.cd0(type_to_use, friction_source, s_wet_source)
    k = prelim_drag.k  # NOTE: add to CD0 bc config and call e + change due to config
    L_D = C_L_climb/(CD0+k*C_L_climb**2)

    CGRP = (CGR+1/L_D)/(np.sqrt(C_L_climb))
    W_P = 18.97*eta_p*np.sqrt(sigma)/CGRP/np.sqrt(W_S)
    return (W_P, W_S)


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
    curve_labels = ["Take-off", "Cruise", "Climb", "Series C", "Turbine"]
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
        W_S_plot: np.ndarray = np.arange(0,10000,5),
        W_P_or_T_W_plot: np.ndarray = np.arange(0,10000,5), 
        output_filepath: str = 'outputs/Matching_Diagram.png') -> list:
    
    # Is it a turbine (does it need those climb lines?)
    turbine = 
    # Compute plots
    stall_W_P_or_T_W, stall_W_S = stall_speed_matching(ac, W_P_or_T_W_plot)
    to_W_P_or_T_W, to_W_S = takeoff_dist_matching(ac, W_S_plot)
    ld_W_P_or_T_W, ld_W_S = landing_dist_matching(ac, W_P_or_T_W_plot)
    cr_W_P_or_T_W, cr_W_S = cruise_speed_matching(ac, type_to_use, friction_source, s_wet_source, W_S_plot)

    cl1_W_P_or_T_W, cl1_W_S = climb_angle_AEO_matching_23_65(ac, type_to_use, friction_source, s_wet_source, W_S_plot)
    cl2_W_P_or_T_W, cl2_W_S = climb_rate_AEO_matching_23_65(ac, type_to_use, friction_source, s_wet_source, W_S_plot)
    cl3_W_P_or_T_W, cl3_W_S = climb_angle_AEO_matching_23_77(ac, type_to_use, friction_source, s_wet_source, W_S_plot)

    if turbine:
        turb1_W_P_or_T_W, turb1_W_S = climb_rate_OEI_multiengine_matching_turbine(ac, type_to_use, friction_source, s_wet_source, W_S_plot)
        turb2_W_P_or_T_W, turb2_W_S = climb_rate_OEI_multiengine_matching_turbine(ac, type_to_use, friction_source, s_wet_source, W_S_plot)

    # Start actual plotting stuff
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(0, 100)
    ax.set_ylim(-10, 10)

    datasets = [
        {"x": stall_W_S, "y": stall_W_P_or_T_W, "label": "Stall speed"},
        {"x": to_W_S, "y": to_W_P_or_T_W, "label": "Take-off field length"},
        {"x": ld_W_S, "y": ld_W_P_or_T_W, "label": "Landing field length"},
        {"x": cr_W_S, "y":cr_W_P_or_T_W, "label": "Cruise speed"},
        {"x": cl1_W_S, "y": cl1_W_P_or_T_W, "label": "Series C"},
        {"x": cl2_W_S, "y": cl2_W_P_or_T_W, "label": "Series C"},
        {"x": cl3_W_S, "y": cl3_W_P_or_T_W, "label": "Series C"}
    ]

    if turbine:
        datasets.append({"x": turb1_W_S, "y": turb1_W_P_or_T_W, "label": "Turbine climb"})
        datasets.append({"x": turb2_W_S, "y": turb2_W_P_or_T_W, "label": "Turbine cruise"})

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