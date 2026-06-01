# Fix path FIRST, before any local imports
import sys
import os

from c2_m import D_Cl_max, S_wf, sweep_at_x_c_deg, chord_at_y_span, closest_value
from class1.c2_drag_new import interp_value

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from classes.aircraft_2 import loader, Aircraft, Requirements, Mission, Weights, Wing, Fuselage, Engine
from lookups.consts import *
from class1.prelim_drag import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines
from classes.isa import Atmosphere
from class1.matching_diagram import plot_matching_and_select_design_point
from class1 import c1_m
from pathlib import Path

def print_combo_table(combo_matrix, margin_matrix, passing_mask,
                    slat_types, flap_types,
                    best_slat_idx, best_flap_idx,
                    DCL_max_to_req, show=True):
    """
    Print a formatted ΔCLmax combination table.

    Parameters
    ----------
    combo_matrix    : (n_slats, n_flaps) ndarray  - combined ΔCLmax values
    margin_matrix   : (n_slats, n_flaps) ndarray  - combo - requirement
    passing_mask    : (n_slats, n_flaps) bool array
    slat_types      : list of slat labels
    flap_types      : list of flap labels
    best_slat_idx   : int or None
    best_flap_idx   : int or None
    DCL_max_to_req  : float - the raw requirement (without buffer)
    show            : bool  - set False to suppress output
    """
    if not show:
        return

    df = pd.DataFrame(
        combo_matrix,
        index=[f"Slat: {s}" for s in slat_types],
        columns=[f"Flap: {f}" for f in flap_types],
    )

    print("\n" + "=" * 60)
    print(f"  ΔCLmax combination table  |  Requirement: {DCL_max_to_req:.4f}")
    print("=" * 60)
    print(df.to_string(float_format=lambda x: f"{x:+.4f}"))

    print("\nMargin above requirement (positive = passes):")
    df_margin = pd.DataFrame(
        margin_matrix,
        index=df.index,
        columns=df.columns,
    )
    print(df_margin.to_string(float_format=lambda x: f"{x:+.4f}"))

    print("\nPassing combinations (margin ≥ 0.05):")
    df_pass = pd.DataFrame(
        passing_mask.astype(int),
        index=df.index,
        columns=df.columns,
    )
    print(df_pass.to_string())

    if best_slat_idx is not None:
        print(f"\n★  Selected: Slat [{best_slat_idx}] {slat_types[best_slat_idx]}"
            f" + Flap [{best_flap_idx}] {flap_types[best_flap_idx]}"
            f"  →  ΔCLmax = {combo_matrix[best_slat_idx, best_flap_idx]:.4f}")
    print("=" * 60 + "\n")

# def lift_slope1(ac: Aircraft, flight_condition): # This lift slope is also already defined in c2_m
#     V = ac.requirements.cruise['cr_speed'] * KTS_TO_MS
#     atm = Atmosphere(ac.requirements.cruise['cr_altitude'] * FT_TO_M, delta_T=0)
#     M = V / np.sqrt(1.4 * 287 * atm.temp)
#     beta = np.sqrt(1 - M**2)                                                        #We also have a beta function defined in c2_m, use that??
#     CLa = 2 * np.pi * ac.wing.aspect_ratio / (2 + np.sqrt(4 + (ac.wing.aspect_ratio * beta / 0.95)**2 * (1 + ac.wing.sweep_c_2_deg**2 / beta**2)))
#     return CLa

def lift_slope_rad(ac: Aircraft, M):
    # V = ac.requirements.cruise['cr_speed'] * KTS_TO_MS
    # atm = Atmosphere(ac.requirements.cruise['cr_altitude'] * FT_TO_M, delta_T=0)
    # M = V / np.sqrt(1.4 * 287 * atm.temp)
    beta = np.sqrt(1 - M**2)                                                        #We also have a beta function defined in c2_m, use that??
    CLa = 2 * np.pi * ac.wing.aspect_ratio / (2 + np.sqrt(4 + (ac.wing.aspect_ratio * beta / 0.95)**2 * (1 + ac.wing.sweep_c_2_deg**2 / beta**2)))
    return CLa

# def high_AR_method1(ac: Aircraft, flight_condition,):      #make this flexible to calculate cl at various airspeeds, thus landing and take-off?
#     V = ac.requirements.cruise['cr_speed'] * KTS_TO_MS
#     atm = Atmosphere(ac.requirements.cruise['cr_altitude'] * FT_TO_M, delta_T=0)
#     M = V / np.sqrt(1.4 * 287 * atm.temp)
#     CL_l_ratio = ...                                                    #Naomi - ratio from figure on slide 22
#     if M > 0.2:                                                         # needs to be local M, thus again dependt on what airspeed
#         CL = CL_l_ratio * ... + DCLmax               #Naomi - CL-l ratio is from figure on slide 22, ... needs to be airfoil clmax at M=0.2, DCLmax function from slide 23, but if we only do take-off and landing then ignore this
#     else:
#         CL = CL_l_ratio * ...                            #Naomi - ... needs to be airfoil clmax at M=0.2
#     return CL

def high_AR_method(ac: Aircraft, DY, M):
    ''' If we apply twist we take zero lift angle of attack at mac'''
    lift_slope = lift_slope_rad(ac, M) * np.pi / 180
    DY_p = round(closest_value(DY, values=[1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.5]), 1)
    DY_pl = round(closest_value(DY, values=[2.001, 2.251, 2.501, 3.001, 4.001, 4.501]), 2)
    CLmax_Clmax = interp_value(pd.read_csv('lookups/Lift/CLmax_Clmax_hihg_ar.csv'), ac.wing.sweep_LE_deg, 'SALE', f'CLL (DELTAY={DY_p:.1f})')
    Delta_CLmax = 0
    if M > 0.2:
        LE = closest_value(ac.wing.sweep_LE_deg, values=[0,20])
        Delta_CLmax = interp_value(pd.read_csv(f'lookups/Lift/dCLmax_high_ar_sweep{LE}.csv'), M, 'AMACH',f'DYA={DY_pl:.2f}')
    CLmax = ac.wing.cl_max_2d * CLmax_Clmax + Delta_CLmax
    
    DY_plt = round(closest_value(DY, values=[1.2, 2, 3, 4]), 1)
    Dalpha_CLmax = interp_value(pd.read_csv('lookups/Lift/dalpha_clmax_high_ar.csv'), ac.wing.sweep_LE_deg, 'SALE', f'DACLL (DY={DY_plt:.1f})')
    alpha_stall = CLmax / lift_slope + ac.wing.alpha_0L + Dalpha_CLmax
    return CLmax, alpha_stall, lift_slope

def low_AR_method(ac: Aircraft, DY, C1, M):
    lift_slope = lift_slope_rad(ac, M) * np.pi / 180
    if ac.wing.x_c_t_c_max<0.35:
        file = 'lookups/Lift/CLmax_base_0.35c.csv'
    else:
        file = 'lookups/Lift/CLmax_base_0.5c.csv'
    C2 = interp_value(pd.read_csv('lookups/Lift/C2.csv'), ac.wing.taper_ratio, 'TR', 'C2')
    M_p = closest_value(M, values=[0.2, 0.4, 0.6])
    x = (C1 + 1) * ac.wing.aspect_ratio / np.sqrt(1 - M**2) * np.cos(np.deg2rad(ac.wing.sweep_LE_deg))
    x1 = (C2 + 1) * ac.wing.aspect_ratio * np.tan(np.deg2rad(ac.wing.sweep_LE_deg))

    v = round(closest_value(DY, values=[0, 0.25, 0.5, 0.75, 1, 1.35]), 2)
    CL_maxbase = interp_value(pd.read_csv(file), x, 'C1ABC', 
                              f'DYAG={v:.2f}')
    DCLmax = interp_value(pd.read_csv('lookups/Lift/dCLmax_low_ar.csv'), x1, 'C2A', f'DE (AMN={M_p})')
    CLmax = CL_maxbase + DCLmax

    alpha_CLmax_b = interp_value(pd.read_csv('lookups/Lift/alpha_CLmax_base.csv'), x, 'T418A', 'D418A')
    if x1<4.5:
        file1 = 'lookups/Lift/dalpha_CLmax1.csv'
        val = ac.wing.aspect_ratio * np.cos(np.deg2rad(ac.wing.sweep_LE_deg)) * (1 + 4 * ac.wing.taper_ratio**2)
        val_p = round(closest_value(val, values=[0, 2, 3, 4, 5, 6, 7, 8, 9, 30]), 1)
        Dalpha_CLmax = interp_value(pd.read_csv(file1), x1, 'D18B2', f'C18B2 (T18B2={val_p:.1f})')
    else:
        file1 = 'lookups/Lift/dalpha_CLmax2.csv'
        Dalpha_CLmax = interp_value(pd.read_csv(file1), x1, 'D18B1', f'C18B1 (T18B1={M_p})')
    alpha_stall = alpha_CLmax_b + Dalpha_CLmax
    return CLmax, alpha_stall, lift_slope

def finit_lift_slope(ac: Aircraft, M):            #CL max of clean wing, needs to be made usable a multiple different speeds
    if ac.wing.airfoil_name == "NACA 23012":
        DY = 26 * ac.wing.t_c_max
        C1 = interp_value(pd.read_csv('lookups/Lift/C1.csv'), ac.wing.taper_ratio, 'taper', 'C1')
        if ac.wing.aspect_ratio > (4 / ((C1 + 1) * np.cos(np.deg2rad(ac.wing.sweep_LE_deg)))):
            CLmax, alpha_stall, lift_slope = high_AR_method(ac, DY, M)
        else:
            CLmax, alpha_stall, lift_slope = low_AR_method(ac, DY, C1, M)
    else:
        raise ValueError(
            'Airfoil changed, so redefine sharpness factor to proceed.'
        )
    return CLmax, alpha_stall, lift_slope


######## Flaps #########
# def D_Cl_max(flap_type: str, # 'plain' or 'split' or 'slotted' or 'fowler' or 'double slotted'or 'triple slotted' or 'fowler'
#              cdash_c: float) -> float:
#     flap_values = {
#         'plain': 0.9,
#         'split': 0.9,
#         'slotted': 1.3,
#         'fowler': 1.3 * cdash_c,
#         'double slotted': 1.6 * cdash_c,
#         'triple slotted': 1.9 * cdash_c,
#         'fixed slot': 0.2,
#         'leading edge flap': 0.3,
#         'kruger': 0.3,
#         'slat': 0.4 * cdash_c
#     }

#     return flap_values.get(flap_type)
def S_wf(ac, y_in, y_out):
    w = ac.wing
    c_in = chord_at_y_span(w.c_root, w.taper_ratio, y_in, w.span)
    c_out = chord_at_y_span(w.c_root, w.taper_ratio, y_out, w.span)
    return (c_in + c_out) * (y_out - y_in)

def flap_deflections_to_ld(flap_type: str,
                     to_or_ld: str = 'landing') -> float:
    flap_values_to = {
        'plain': 20,
        'slotted': 20,
        'fowler': 15,
        'double slotted': 20,
        'triple slotted': 20,
    }
    flap_values_ld = {
        'plain': 60,
        'slotted': 40,
        'fowler': 40,
        'double slotted': 50,
        'triple slotted': 40,
    }
    if to_or_ld == 'landing':
        return flap_values_ld.get(flap_type)
    else:
        return flap_values_to.get(flap_type)

def deltaC_Cf(flap_type: str, flap_deflection_deg: float) -> float:
    if flap_type == 'double slotted' or flap_type == 'triple slotted' or flap_type == 'fowler':
        line = '4b'
    elif flap_type == 'slotted':
        line = '2a'
    elif flap_type == 'plain':
        line = '1b'
    else:
        raise ValueError(f'Flap type given: {flap_type}, if not in [plain, slotted, fowler, double slotted, triple slotted]')
    return interp_value(pd.read_csv('lookups/HLD/deltaC_Cf.csv'), flap_deflection_deg, 'df', line)

def Delta_CLmax_and_alpha0L_deg_flaps(ac: Aircraft, flap_type: str, y_in, y_out, flight_condition: str = 'landing'):
    hinge_line_x_c = ac.hld_and_ailerons.flaps['hinge_line_x_c']
    flap_deflection = flap_deflections_to_ld(flap_type, flight_condition)
    if flap_type not in ['slat', 'fixed slot', 'kruger', 'leading edge flap']:
        cdash_c = (deltaC_Cf(flap_type, flap_deflection) + 1 / ac.hld_and_ailerons.flaps['cf_c']) * ac.hld_and_ailerons.flaps['cf_c']
    else:
        cdash_c = 1.05
    DClmax = D_Cl_max(flap_type, cdash_c)
    Swf = S_wf(ac, y_in, y_out)
    hinge_sweep = sweep_at_x_c_deg(ac.wing.sweep_LE_deg, ac.wing.c_root, ac.wing.span, ac.wing.taper_ratio, hinge_line_x_c)
    DCLmax = 0.9 * DClmax * Swf / ac.wing.area * np.cos(np.deg2rad(hinge_sweep))
    Dalpha0L = -10 * Swf / ac.wing.area * np.cos(np.deg2rad(hinge_sweep))
    if flight_condition == 'landing':
        return DCLmax, Dalpha0L*1.5
    else: 
        return DCLmax, Dalpha0L
    
def flapped_wing_slope(ac: Aircraft, wing_lift_slope, Swf, flap_type: str, flight_condition: str = 'landing'):
    flap_deflection = flap_deflections_to_ld(flap_type, flight_condition)
    if flap_type not in ['slat', 'fixed slot', 'kruger', 'leading edge flap']:
        cdash_c = (deltaC_Cf(flap_type, flap_deflection) + 1 / ac.hld_and_ailerons.flaps['cf_c']) * ac.hld_and_ailerons.flaps['cf_c']
    else:
        cdash_c = 1.05
    if flap_type == 'double slotted' or flap_type == 'triple slotted' or flap_type == 'fowler':
        return wing_lift_slope * (1 + Swf / ac.wing.area * (cdash_c - 1))
    else:
        return wing_lift_slope
    

# ailerons
def Aileron_control_dtive(ac: Aircraft, y_in, y_out):
    w = ac.wing
    cf_c = ac.hld_and_ailerons.flaps['cf_c']
    tau = interp_value(pd.read_csv('lookups/HLD/tau.csv'), cf_c, 'ratio', 'tau')
    integral = w.c_root / 2 * (y_out**2 - y_in**2) - (w.c_root - w.c_tip) * 2 / (3 * w.span) * (y_out**3 - y_in**3)
    return 2 * 2 * np.pi * tau * integral / (w.area * w.span)
def Roll_damping_coefficient(ac: Aircraft):
    w = ac.wing
    integral = w.c_root / 3 * ((w.span / 2)**3) - (w.c_root - w.c_tip) / (2 * w.span) * ((w.span / 2)**4)
    return -4 * (2 * np.pi + w.cd0) * integral / (w.area * w.span**2)
def Roll_rate_deg(ac: Aircraft, speed, aileron_defl_deg, Clda, Clp):
    return (-Clda / Clp * aileron_defl_deg / 180 * np.pi * 2 * speed / ac.wing.span) / np.pi * 180

def size_ailerons(ac: Aircraft, update_ac: bool = False):
    # Sets outer position of aileron along span and finds required inner position based on landing and take-off CS-23 requirements
    w = ac.wing
    y_out = 0.98 * w.span
    y_in = np.arange(0.5, 1, 0.01) * w.span
    Clp = Roll_damping_coefficient(ac)
    Clda = Aileron_control_dtive(ac, y_in, y_out)
    # Calc roll rates at TO and LD speeds required
    P_ld = Roll_rate_deg(ac, speed=ac.requirements.general['stall_speed'] * KTS_TO_MS * 1.3, 
                         aileron_defl_deg=ac.hld_and_ailerons.ailerons['deflection'], Clda=Clda, Clp=Clp)
    P_to = Roll_rate_deg(ac, speed=ac.requirements.general['stall_speed'] * KTS_TO_MS * 1.2, 
                         aileron_defl_deg=ac.hld_and_ailerons.ailerons['deflection'], Clda=Clda, Clp=Clp)
    # Locate y_in position required
    idx_ld = next(i for i, v in enumerate(P_ld) if v < 60/4*1.03)
    min_y_ld = y_in[idx_ld]
    idx_to = next(i for i, v in enumerate(P_to) if v < 60/5*1.03)
    min_y_to = y_in[idx_to]

    if update_ac:
        ac.hld_and_ailerons.ailerons['y_aileron_in'] = min(min_y_ld, min_y_to)
        ac.hld_and_ailerons.ailerons['y_aileron_out'] = y_out
    return min(min_y_ld, min_y_to), y_out

def size_HLD(ac: Aircraft, update_ac: bool = False):
    y_in = ac.fuselage.width / 2 + 0.3
    y_out = ac.hld_and_ailerons.ailerons['y_aileron_in'] - 0.05
    y_out_slat = 0.96 * ac.wing.span

    Swf = S_wf(ac, y_in, y_out)
    Swf_slat = S_wf(ac, y_in, y_out_slat)
    alpha_0L = ac.wing.alpha_0L

    flap_types = ['plain', 'slotted', 'fowler', 'double slotted', 'triple slotted']
    slat_types = ['fixed slot', 'leading edge flap', 'kruger', 'slat']

    # Cruise
    Atm_cr = Atmosphere(ac.requirements.cruise['cr_altitude'], delta_T=0)
    T_cr = float(Atm_cr.temp)
    M_cr = ac.requirements.cruise['cr_speed'] * KTS_TO_MS / np.sqrt(1.4 * 287 * T_cr)
    CLmax_cr, alpha_stall_cr, lift_slope_cr = finit_lift_slope(ac, M_cr)

    # TO
    Atm_to = Atmosphere(ac.requirements.take_off['to_altitude'], ac.requirements.take_off['to_temp_shift'])
    rho_to = float(Atm_to.density)
    T_to = float(Atm_to.temp)
    CLmax_to_as = ac.requirements.take_off['as_CL_max_to']
    V_to = np.sqrt(ac.weights.m_takeoff * 9.81 / (0.5 * rho_to * ac.wing.area * CLmax_to_as / 1.21))
    M_to = V_to / np.sqrt(1.4 * 287 * T_to)

    # Set up table
    CLmax_to, alpha_stall_to, lift_slope_to = finit_lift_slope(ac, M_to)
    flap_DCLmax_to = np.zeros(5)
    slat_DCLmax_to = np.zeros(4)
    for i, flap in enumerate(flap_types):
        flap_DCLmax_to[i], _ = Delta_CLmax_and_alpha0L_deg_flaps(ac, flap, y_in, y_out, flight_condition='take-off')
        wing_lift_slope = flapped_wing_slope(ac, lift_slope_to, Swf, flap, flight_condition='take-off')
    for i, slat in enumerate(slat_types):
        slat_DCLmax_to[i], _ = Delta_CLmax_and_alpha0L_deg_flaps(ac, slat, y_in, y_out_slat, flight_condition='take-off')
        wing_lift_slope = flapped_wing_slope(ac, lift_slope_to, Swf_slat, slat, flight_condition='take-off')
    DCL_max_to_req = CLmax_to_as - CLmax_to

    # 4×5 slat and flap combination matrix
    # combo_matrix[i, j] = slat_DCLmax[i] + flap_DCLmax[j]
    combo_matrix = slat_DCLmax_to[:, None] + flap_DCLmax_to[None, :]  # (4, 5)

    # Margin matrix: 
    margin_matrix = combo_matrix - DCL_max_to_req  # positive -> meets req.

    # Find the first combo that passes by at least 0.05
    MARGIN_BUFFER = 0.05
    passing_mask = margin_matrix >= MARGIN_BUFFER  # bool (4, 5)

    # Iterate in row-major order (slat 0 -> 3, flap 0 -> 4 within each slat row)
    best_slat_idx, best_flap_idx = None, None
    for i in range(len(slat_types)):
        for j in range(len(flap_types)):
            if passing_mask[i, j]:
                best_slat_idx, best_flap_idx = i, j
                break
        if best_slat_idx is not None:
            break

    if best_slat_idx is not None:
        print(f"First combination that meets ΔCLmax_req + {MARGIN_BUFFER:.2f}:")
        print(f"  Slat : {slat_types[best_slat_idx]}  (index {best_slat_idx})")
        print(f"  Flap : {flap_types[best_flap_idx]}  (index {best_flap_idx})")
        print(f"  ΔCLmax = {combo_matrix[best_slat_idx, best_flap_idx]:.4f}  "
            f"(req = {DCL_max_to_req:.4f}, margin = "
            f"{margin_matrix[best_slat_idx, best_flap_idx]:.4f})")
    else:
        print("No combination meets the requirement with the specified margin.")

    # Printint table
    print_combo_table(
        combo_matrix, margin_matrix, passing_mask,
        slat_types, flap_types,
        best_slat_idx, best_flap_idx,
        DCL_max_to_req,
        show=True,
    )

    # LD
    CLmax_ld_as = ac.requirements.landing['as_CL_max_la']
    Atm_ld = Atmosphere(ac.requirements.landing['la_altitude'], ac.requirements.landing['la_temp_shift'])
    T_ld = float(Atm_ld.temp)
    M_ld = ac.requirements.general['stall_speed'] * KTS_TO_MS * 1.3 / np.sqrt(1.4 * 287 * T_ld)
    # Set up table
    CLmax_ld, alpha_stall_ld, lift_slope_ld = finit_lift_slope(ac, M_ld)
    flap_DCLmax_ld = np.zeros(5)
    slat_DCLmax_ld = np.zeros(4)
    for i, flap in enumerate(flap_types):
        flap_DCLmax_ld[i], Dalpha0L = Delta_CLmax_and_alpha0L_deg_flaps(ac, flap, y_in, y_out, flight_condition='landing')
        wing_lift_slope = flapped_wing_slope(ac, lift_slope_ld, Swf, flap, flight_condition='landing')
    for i, slat in enumerate(slat_types):
        slat_DCLmax_ld[i], Dalpha0L = Delta_CLmax_and_alpha0L_deg_flaps(ac, slat, y_in, y_out_slat, flight_condition='landing')
        wing_lift_slope = flapped_wing_slope(ac, lift_slope_ld, Swf_slat, slat, flight_condition='landing')
    DCL_max_ld_req = CLmax_ld_as - CLmax_ld
    # 4×5 slat and flap combination matrix
    # combo_matrix[i, j] = slat_DCLmax[i] + flap_DCLmax[j]
    combo_matrix_ld = slat_DCLmax_ld[:, None] + flap_DCLmax_ld[None, :]  # (4, 5)

    # Margin matrix: 
    margin_matrix_ld = combo_matrix_ld - DCL_max_ld_req  # positive -> meets req.

    # Find the first combo that passes by at least 0.05
    MARGIN_BUFFER = 0.05
    passing_mask_ld = margin_matrix_ld >= MARGIN_BUFFER  # bool (4, 5)

    # Iterate in row-major order (slat 0 -> 3, flap 0 -> 4 within each slat row)
    best_slat_idx_ld, best_flap_idx_ld = None, None
    for i in range(len(slat_types)):
        for j in range(len(flap_types)):
            if passing_mask_ld[i, j]:
                best_slat_idx_ld, best_flap_idx_ld = i, j
                break
        if best_slat_idx_ld is not None:
            break

    if best_slat_idx_ld is not None:
        print(f"First combination that meets ΔCLmax_req + {MARGIN_BUFFER:.2f}:")
        print(f"  Slat : {slat_types[best_slat_idx_ld]}  (index {best_slat_idx_ld})")
        print(f"  Flap : {flap_types[best_flap_idx_ld]}  (index {best_flap_idx_ld})")
        print(f"  ΔCLmax = {combo_matrix_ld[best_slat_idx_ld, best_flap_idx_ld]:.4f}  "
            f"(req = {DCL_max_ld_req:.4f}, margin = "
            f"{margin_matrix_ld[best_slat_idx_ld, best_flap_idx_ld]:.4f})")
    else:
        print("No combination meets the requirement with the specified margin.")

    # Printint table
    print_combo_table(
        combo_matrix_ld, margin_matrix_ld, passing_mask_ld,
        slat_types, flap_types,
        best_slat_idx_ld, best_flap_idx_ld,
        DCL_max_ld_req,
        show=True,
    )

    flap_index = int(input('Select the index of the flaps required'))
    slat_index = int(input('Select the index of the slats required'))

    flap_DCLmax_to, Dalpha0L_fto = Delta_CLmax_and_alpha0L_deg_flaps(ac, flap, y_in, y_out, flight_condition='take-off')
    wing_lift_slope_flapped_to = flapped_wing_slope(ac, lift_slope_to, Swf, flap, flight_condition='take-off')
    slat_DCLmax_to, Dalpha0L_sto = Delta_CLmax_and_alpha0L_deg_flaps(ac, slat, y_in, y_out_slat, flight_condition='take-off')
    
    flap_DCLmax_ld, Dalpha0L_fld = Delta_CLmax_and_alpha0L_deg_flaps(ac, flap, y_in, y_out, flight_condition='landing')
    wing_lift_slope_flapped_ld = flapped_wing_slope(ac, lift_slope_to, Swf, flap, flight_condition='landing')
    slat_DCLmax_ld, Dalpha0L_sld = Delta_CLmax_and_alpha0L_deg_flaps(ac, slat, y_in, y_out_slat, flight_condition='landing')
    
    if update_ac:
        ac.hld_and_ailerons.flaps['flap_type'] = flap_types[flap_index]
        if flap_types != 'slat':
            ac.hld_and_ailerons.flaps['Delta_c_cf_to'] = deltaC_Cf(flap_types[flap_index], flap_deflections_to_ld(flap_types[flap_index], 'take-off'))
            ac.hld_and_ailerons.flaps['Delta_c_cf_ld'] = deltaC_Cf(flap_types[flap_index], flap_deflections_to_ld(flap_types[flap_index], 'landing'))
            ac.hld_and_ailerons.flaps['cdash_cf_to'] = ac.hld_and_ailerons.flaps['Delta_c_cf_to'] + 1 / ac.hld_and_ailerons.flaps['cf_c']
            ac.hld_and_ailerons.flaps['cdash_cf_ld'] = ac.hld_and_ailerons.flaps['Delta_c_cf_ld'] + 1 / ac.hld_and_ailerons.flaps['cf_c']
            ac.hld_and_ailerons.flaps['cdash_c_to'] = ac.hld_and_ailerons.flaps['cdash_cf_to'] * ac.hld_and_ailerons.flaps['cf_c']
            ac.hld_and_ailerons.flaps['cdash_c_ld'] = ac.hld_and_ailerons.flaps['cdash_cf_ld'] * ac.hld_and_ailerons.flaps['cf_c']
        else:
            ac.hld_and_ailerons.flaps['Delta_c_cf_to'] = 0
            ac.hld_and_ailerons.flaps['Delta_c_cf_ld'] = 0
            ac.hld_and_ailerons.flaps['cdash_c_to'] = 1.05
            ac.hld_and_ailerons.flaps['cdash_c_ld'] = 1.05
            ac.hld_and_ailerons.flaps['cdash_cf_to'] = 0
            ac.hld_and_ailerons.flaps['cdash_cf_ld'] = 0
        ac.hld_and_ailerons.flaps['S_wf'] = Swf
        ac.hld_and_ailerons.flaps['y_flap_in'] = y_in
        ac.hld_and_ailerons.flaps['y_flap_out'] = y_out_slat
        ac.hld_and_ailerons.flaps['ld_deflection'] = flap_deflections_to_ld(flap_types[flap_index], 'landing')
        ac.hld_and_ailerons.flaps['to_deflection'] = flap_deflections_to_ld(flap_types[flap_index], 'take-off')
        ac.hld_and_ailerons.slats['slat_type'] = slat_types[slat_index]
        ac.hld_and_ailerons.slats['y_slat_in'] = y_in
        ac.hld_and_ailerons.slats['y_slat_out'] = y_out_slat
        if slat_types[slat_index] == 'slat':
            ac.hld_and_ailerons.slats['cdash_c_to'] = 1.05
            ac.hld_and_ailerons.slats['cdash_c_ld'] = 1.05
        else:
            ac.hld_and_ailerons.slats['cdash_c_to'] = (deltaC_Cf(slat_types[slat_index], flap_deflections_to_ld(slat_types[slat_index], 'take-off')) + 1 / ac.hld_and_ailerons.flaps['cf_c']) * ac.hld_and_ailerons.flaps['cf_c']
            ac.hld_and_ailerons.slats['cdash_c_ld'] = (deltaC_Cf(slat_types[slat_index], flap_deflections_to_ld(slat_types[slat_index], 'landing')) + 1 / ac.hld_and_ailerons.flaps['cf_c']) * ac.hld_and_ailerons.flaps['cf_c']
        ac.hld_and_ailerons.slats['S_wf'] = Swf_slat
        ac.hld_and_ailerons.landing_lift['CL_alpha'] = wing_lift_slope_flapped_ld
        ac.hld_and_ailerons.landing_lift['alpha_zero_lift'] = alpha_stall_ld + Dalpha0L_fld + Dalpha0L_sld
        ac.hld_and_ailerons.landing_lift['CL_max'] = CLmax_ld + flap_DCLmax_ld + slat_DCLmax_ld
        ac.hld_and_ailerons.take_off_lift['CL_alpha'] = wing_lift_slope_flapped_to
        ac.hld_and_ailerons.take_off_lift['alpha_zero_lift'] = alpha_stall_to + Dalpha0L_fto + Dalpha0L_sto
        ac.hld_and_ailerons.take_off_lift['CL_max'] = CLmax_to + flap_DCLmax_to + slat_DCLmax_to
        ac.hld_and_ailerons.clean_lift['CL_alpha'] = lift_slope_cr
        ac.hld_and_ailerons.clean_lift['CL_max'] = CLmax_cr
        ac.hld_and_ailerons.clean_lift['alpha_stall'] = alpha_stall_cr
    return flap_types[flap_index], slat_types[slat_index], y_in, y_out, y_out_slat




# #### Take off
# def flap_dc_cf_TO(flap_type: str, # 'plain' or 'slotted' or 'fowler' or 'double slotted'or 'triple slotted' or 'fowler'
#              ):
#     flap_values = {
#         'plain': 0.15,
#         'slotted': 0.2,
#         'fowler': 0.5,
#         'double slotted': 0.6,
#         'triple slotted': 0.6,
#     }
#     return flap_values.get(flap_type)

# def cdash_c_TO(ac: Aircraft, flap_type: str):
#     if flap_type == "slat":
#         cdash_c = 1.05
#     else:
#         cdash_cf = flap_dc_cf_TO(flap_type) + (1 / ac.hld_and_ailerons.flaps['cf_c'])
#         cdash_c = cdash_cf * ac.hld_and_ailerons.flaps['cf_c']
#     return cdash_c


# def D_CL_max_TO(ac: Aircraft, flap_type: str, y_start_f, y_end_f, taper, c_r, b, LE_sweep, x_c_hinge):         #Naomi x_c_hinge needs to be defined somewhere
#     return 0.8 * 0.9 * D_Cl_max(flap_type) * (S_wf(y_start_f, y_end_f, taper, c_r, b)/ac.wing['area']) * np.cos(sweep_at_x_c_deg(LE_sweep, c_r, b, taper, x_c_hinge))             # Naomi fill in functions arguements, also this function is the same for slats...jsut output in yaml is under different heading


# ##### Landing



# def flap_dc_cf_LA(flap_type: str, # 'plain' or 'slotted' or 'fowler' or 'double slotted'or 'triple slotted' or 'fowler'
#              ):
#     flap_values = {
#         'plain': 0.45,
#         'slotted': 0.3,
#         'fowler': 0.6,
#         'double slotted': 0.85,
#         'triple slotted': 0.85,
#     }

#     return flap_values.get(flap_type)

# def cdash_c_LA(flap_type: str, cf_c):
#     if flap_type == "slat":
#         cdash_c = 1.05
#     else:
#         cdash_cf = flap_dc_cf_LA(flap_type) + (1/cf_c)
#         cdash_c = cdash_cf * cf_c
#     return cdash_c


# def D_CL_max_LA(ac: Aircraft, flap_type: str, y_start_f, y_end_f, taper, c_r, b, LE_sweep, x_c_hinge):         #Naomi x_c_hinge needs to be defined somewhere
#     return 0.9 * D_Cl_max(flap_type) * (S_wf(y_start_f, y_end_f, taper, c_r, b)/ac.wing['area']) * np.cos(sweep_at_x_c_deg(LE_sweep, c_r, b, taper, x_c_hinge))             # Naomi fill in functions arguements


# # Now how to check if it meets the CLmax demands, implement that
# #make a matrix to pick and choose combination

# ######## Aileron #############

# # define aileron geom
# # y_start_a
# # y_end_a
# # cf_c
# # alieron max deflect - check what value we want
# # does the roll rate make sense, which req
# # define x_c_hinge line


# def aileron_control_derivative(ac: Aircraft,
#     Cl_alpha,
#     tau,
#                                c_r,
#                                taper,
#                                b,
#     n_points=1000
# ):
#     """
#     Compute Cl_delta_a using a DATCOM-style integration.

#     Parameters
#     ----------
#     CL_alpha : float
#         Wing lift curve slope [1/rad]

#     y1 : float
#         Inboard aileron span station [m]

#     y2 : float
#         Outboard aileron span station [m]

#     chord_function : function
#         Function c(y) returning local chord length
#     """

#     y = np.linspace(ac.hld_and_ailerons.ailerons['y_aileron_in'], ac.hld_and_ailerons.ailerons['y_aileron_out'], n_points)

#     c = chord_at_y_span(c_r, taper, y, b)

#     integrand = c * y

#     integral = np.trapz(integrand, y)

#     Cl_da = (2 * Cl_alpha * tau / (ac.wing['area'] * ac.wing['span'])) * integral

#     return Cl_da

# def roll_damping_derivative(ac: Aircraft,
#     Cl_alpha,
#                             cd0,
#                                c_r,
#                                taper,
#                                b,
#     n_points=1000
# ):

#     y = np.linspace(ac.hld_and_ailerons.ailerons['y_aileron_in'], ac.hld_and_ailerons.ailerons['y_aileron_out'], n_points)

#     c = chord_at_y_span(c_r, taper, y, b)

#     integrand = c * y **2

#     integral = np.trapz(integrand, y)

#     Cl_p = - ((4 * (Cl_alpha * cd0)**(b/2)) / (ac.wing['area'] * ac.wing['span']**2)) * integral

#     return Cl_p


# def AC_roll_rate(ac: Aircraft,):
#     P = - (aileron_control_derivative()/roll_damping_derivative()) * ... * (2 * ... /ac.wing['span']) #input max aileron deflection and chosen speeds
#     return P

# #check is meets the roll rate requirements, for different speeds


# ##### CL-a curve for flapped wing


# def Sdash_S(ac: Aircraft, S, cdash_c):
#     return 1 + (ac.hld_and_ailerons.flaps['S_wf']/S) * (cdash_c - 1)            #check if inputs make sense to call from yamls or not


# def lift_slope_flapped(ac: Aircraft):
#     return Sdash_S() * lift_slope()                                 #maybe want to put Sdash_S as a yaml output, same wiht lfit slope?

# def delta_zerolift(ac: Aircraft, flight_condition, LE_sweep, c_r, b, taper, x_c_hinge):
#     if flight_condition == "landing":
#         d_0lift_airfoil = - 15
#     elif flight_condition == "take_off":
#         d_0lift_airfoil = - 10
#     else:
#         raise ValueError(
#             'Invalid flight condition for this function, please select another one.'
#         )
#     return d_0lift_airfoil * (ac.hld_and_ailerons.flaps['S_wf']/ac.wing['area']) * np.cos(sweep_at_x_c_deg(LE_sweep, c_r, b, taper, x_c_hinge))         #check inputs, last part is sweep at hinge line