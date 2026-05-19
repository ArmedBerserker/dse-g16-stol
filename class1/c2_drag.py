from classes.aircraft_2 import Aircraft
from c2_m import Snet, LE_sweep_deg, sweep_at_x_c_deg, lift_slope, S_wf, closest_value
from lookups.consts import *
import pandas as pd
import numpy as np
from scipy.interpolate import RegularGridInterpolator, interp1d
from scipy.optimize import brentq
from classes.isa import Atmosphere
import matplotlib.pyplot as plt

def mu_air(T):
    return 1.81 * 1e-5 * (T / 293.15)**1.5 * (293.15 + 110.4) / (T + 110.4)

def exposed_wing_mgc(S, b, taper, c_r, b_f):
    S_net = Snet(S, b_f, taper, b, c_r)
    return S_net / (b - b_f)

def CD0_w(R_wf, R_LS, C_f_w, l_dash, t_c_max, S_wet_w, S):
    return R_wf * R_LS * C_f_w * (1 + l_dash * t_c_max + 100 * t_c_max**4) * S_wet_w / S

def CD0_fuselage(density, speed, l_f, mu, S_fus, # max fuselage cross section area
                 S_b_fus, S, M, R_N_fus):
    R_wf = Rwf(M, R_N_fus)
    C_f_fus = C_f(R_N=density * speed * l_f / mu, M=M)
    d_b = np.sqrt(4 * S_b_fus / np.pi)
    d_f_eq = np.sqrt(4 * S_fus / np.pi)
    S_wet_fus =  # NOTE: add method
    C_D0_b_fus = R_wf * C_f_fus * (1 + 60 / (l_f / d_f_eq)**3 + 0.0025 * (l_f / d_f_eq)) * S_wet_fus / S
    C_D_b_fuse = (0.029 * (d_b / d_f_eq)**3 / (C_D0_b_fus * (S / S_fus))**0.5) * (S / S_fus)
    return C_D0_b_fus + C_D_b_fuse

def S_wet_wing(t_c_r, t_c_t, S, b_f, taper, b, c_r):
    tau = t_c_t / t_c_r
    return 2 * Snet(S, b_f, taper, b, c_r) * (1 + 0.25 * t_c_r * (1 + tau * taper) / (1 + taper))

def C_f(R_N, M):  # DATCOM fig 4.1.5.1-26
    x = np.log10(R_N)
    CF3 = 3.92725e-6 * x**5 -1.30370e-4 * x**4 -1.65388e-3 * x**3 -9.59519e-3 * x**2 +2.18366e-2 * x 
    CF0 = 4.12963e-6 * x**5 -1.36204e-4 * x**4 + 1.71620e-3 * x**3 - 9.88935e-3 * x**2 +2.23641e-2 * x
    if M > 0.3:
        return CF3
    else: 
        return CF0 - (CF0 - CF3) * (M / 0.3)

def flap_interference_factor(flap_type: str, # 'split' or 'plain' or 'slotted' 'fowler' or 'krueger'
                             ):
    if flap_type == 'split':
        return -0.15
    elif flap_type == 'plain':
        return 0
    elif flap_type == 'slotted':
        return 0.4
    elif flap_type == 'fowler':
        return 0.25
    else:
        return 0.1
    
def interp_value(df : pd.DataFrame,
                 x_query,
                 x_col : str,
                 y_col : str,
                 log_x = False) -> float:
    data = df[[x_col, y_col]].dropna().sort_values(x_col)
 
    x = data[x_col].to_numpy(dtype = float)
    y = data[y_col].to_numpy(dtype = float)
 
    if log_x:
        x = np.log10(x)
        x_query = np.log10(x_query)
       
 
    return float(np.interp(x_query, x, y))

def Rwf(M, R_N_fus):
    M = closest_value(M, values=[0.25, 0.4])
    if M == 0.25:
        return interp_value(pd.read_csv('lookups/roskam_p6_fig_4_1_rwf.csv'), R_N_fus, x_col='Fuselage Reynolds Number (M = 0.25)', y_col='Wing-Fuse Interference Factor (M = 0.25)', log_x=True)
    else:
        return interp_value(pd.read_csv('lookups/roskam_p6_fig_4_1_rwf_2.csv'), R_N_fus, x_col='Fuselage Reynolds Number (M = 0.4)', y_col='Wing-Fuse Interference Factor (M = 0.4)', log_x=True)
    

def C_D0(ac: Aircraft, 
         n_engine_operative: int,
         flap_type: str, # 'split' or 'plain' or 'slotted' 'fowler' or 'krueger'
         flight_condition: str = 'cruise', # 'cruise' or 'landing' or 'take-off'
         nacelle_on_top_of_wing: bool = True, 
         ):
    S = ac.wing.area
    if flight_condition == 'cruise':
        temp_shift = 0
        alt = ac.requirements.cruise['cr_altitude'] * FT_TO_M
        speed = ac.requirements.cruise['cr_speed'] * KTS_TO_MS
        Atm = Atmosphere(alt, temp_shift)
        temp = Atm.temp
        density = Atm.density
    if flight_condition == 'take-off':
        C_L = 
        temp_shift = ac.requirements.take_off['to_temp_shift']
        alt = ac.requirements.take_off['to_altitude'] * FT_TO_M
        Atm = Atmosphere(alt, temp_shift)
        density = Atm.density
        mass_frac = ac.requirements.cruise['to_mass_frac']
        speed = np.sqrt(mass_frac * ac.weights.m_takeoff / (0.5 * density * C_L * S))
    if flight_condition == 'take-off':
        temp_shift = ac.requirements.landing['la_temp_shift']
        alt = ac.requirements.landing['la_altitude'] * FT_TO_M
        Atm = Atmosphere(alt, temp_shift)
        density = Atm.density
        mass_frac = ac.requirements.landing['la_mass_frac']
        speed = 1.3 * ac.requirements.general['stall_speed'] * KTS_TO_MS
        
    # Wing
    mu = mu_air(temp)
    M = speed / np.sqrt(1.4 * 287 * temp)
    b = ac.wing.span
    sweep_c_4_deg = ac.wing.sweep
    taper = ac.wing.taper_ratio
    c_r = ac.wing.c_root
    b_f = ac.fuselage.width
    l_f = ac.fuselage.length
    c_w_e = exposed_wing_mgc(S, b, taper, c_r, b_f)
    # R_N_w = density * speed * c_w_e / mu
    R_N_fus = density * speed * l_f / mu
    R_wf = Rwf(M, R_N_fus)
    C_f_w = C_f(R_N=density * speed * c_w_e / mu, M=M)
    l_dash = 1.2
    t_c_max = ac.wing.t_c_max # at mean geometric chord
    x_c_t_c_max = 
    if x_c_t_c_max < 0.3:
        l_dash = 2.0
    t_c_r = t_c_max
    t_c_t = t_c_max
    sweep_t_c_max_deg = sweep_at_x_c_deg(ac.wing.sweep_LE_deg, c_r, b, taper, x_c=x_c_t_c_max)
    R_LS = interp_value(pd.read_csv('lookups/roskam_p6_fig_4_2_rls.csv'), np.cos(np.deg2rad(sweep_t_c_max_deg)), 'cos(quarter chord)', 'lifting surface correction', log_x=False)
    S = ac.wing.area

    wing = CD0_w(R_wf, R_LS, C_f_w, l_dash, t_c_max, S_wet_wing(t_c_r, t_c_t, S, b_f, taper, b, c_r), S)  # R_wf * R_LS * C_f_w * (1 + l_dash * t_c_max + 100 * t_c_max**4) * S_wet_w / S

    # Fuselage
    S_fus =  # max fuselage cross section area
    S_b_fus = 
    fuselage = CD0_fuselage(density, speed, l_f, mu, S_fus, S_b_fus, S, M, R_N_fus=R_N_fus)

    # HT
    t_tail_condition: bool =  # t-tail or not
    ht_sweep_c_4_deg = 
    c_r_ht = 
    taper_ht = 
    b_ht = 
    S_ht = 
    S_wet_ht = 
    ht_t_c_max = 
    ht_x_c_t_c_max = 
    ht_sweep_t_c_max_deg = sweep_at_x_c_deg(ac.wing.sweep_LE_deg, c_r_ht, b_ht, taper_ht, x_c=ht_x_c_t_c_max)
    R_LS_ht = interp_value(pd.read_csv('lookups/roskam_p6_fig_4_2_rls.csv'), np.cos(np.deg2rad(ht_sweep_t_c_max_deg)), 'cos(quarter chord)', 'lifting surface correction', log_x=False)
    b_f_ht =  # fuselage width at ht intersection position
    R_N_ht = density * speed * exposed_wing_mgc(S_ht, b_ht, taper_ht, c_r_ht, b_f_ht) / mu
    if not t_tail_condition:
        R_N_ht *= np.sqrt(0.85)
    C_f_ht = C_f(R_N_ht, M)
    l_dash_ht = 1.2
    if ht_x_c_t_c_max < 0.3:
        l_dash_ht = 2.0
    t_c_r_ht = 
    t_c_t_ht = 
    ht = CD0_w(R_wf=1.0, R_LS=R_LS_ht, C_f_w=C_f_ht, l_dash=l_dash_ht, t_c_max=ht_t_c_max, S_wet_w=S_wet_wing(t_c_r_ht, t_c_t_ht, S_ht, b_f_ht, taper_ht, b_ht, c_r_ht), S=S_ht)

    # VT
    vt_sweep_c_4_deg = 
    c_r_vt = 
    taper_vt = 
    b_vt = 
    S_vt = 
    vt_t_c_max = 
    vt_x_c_t_c_max = 
    vt_sweep_t_c_max_deg = sweep_at_x_c_deg(LE_sweep_deg(vt_sweep_c_4_deg, c_r_vt, b_vt, taper_vt), c_r_vt, b_vt, taper_vt, x_c=vt_x_c_t_c_max)
    R_LS_vt = interp_value(pd.read_csv('lookups/roskam_p6_fig_4_2_rls.csv'), np.cos(np.deg2rad(vt_sweep_t_c_max_deg)), 'cos(quarter chord)', 'lifting surface correction', log_x=False)
    b_f_vt =  # fuselage width at vt intersection position
    C_f_vt = C_f(R_N=density * speed * exposed_wing_mgc(S_vt, b_vt, taper_vt, c_r_vt, b_f_vt) / mu, M=M)
    l_dash_vt = 1.2
    if vt_x_c_t_c_max < 0.3:
        l_dash_vt = 2.0
    t_c_r_vt = 
    t_c_t_vt = 
    vt = CD0_w(R_wf=1.0, R_LS=R_LS_vt, C_f_w=C_f_vt, l_dash=l_dash_vt, t_c_max=vt_t_c_max, S_wet_w=S_wet_wing(t_c_r_vt, t_c_t_vt, S_vt, b_f_vt, taper_vt, b_vt, c_r_vt), S=S_vt)

    # Nacelle/pylon
    n_eng = ac.engine.count
    l_nac =  # nacelle length
    S_nac_max =  # max nacelle cross section area
    S_b_nac =  # nacelle base area
    isolated_nac = CD0_fuselage(density, speed, l_nac, mu, S_nac_max, S_b_nac, S, M, R_N_fus)

    c_r_nac = 
    b_nac = 
    nac_t_c_max = 
    nac_x_c_t_c_max = 
    nac_sweep_t_c_max_deg = sweep_at_x_c_deg(LE_sweep_deg(sweep_c_4=0, c_r=c_r_nac, b=b_nac, taper_ratio=0), c_r_nac, b_nac, taper_ratio=0, x_c=nac_x_c_t_c_max)
    R_LS_nac = interp_value(pd.read_csv('lookups/roskam_p6_fig_4_2_rls.csv'), np.cos(np.deg2rad(nac_sweep_t_c_max_deg)), 'cos(quarter chord)', 'lifting surface correction', log_x=False)
    C_f_nac = C_f(R_N=density * speed * c_r_nac / mu, M=M)
    l_dash_nac = 1.2
    if nac_x_c_t_c_max < 0.3:
        l_dash_nac = 2.0
    isolated_pylon = CD0_w(R_wf=1.0, R_LS=R_LS_nac, C_f_w=C_f_nac, l_dash=l_dash_nac, t_c_max=nac_t_c_max, S_wet_w=S_wet_wing(nac_t_c_max, nac_t_c_max, S=b_nac*c_r_nac, b_f=t_c_max*(c_r + taper*c_r)/2, taper=0, b=b_nac, c_r=c_r_nac))

    b_nac =  # nacelle width
    c_nac =  # chord at nacelle 
    i_n =  # nacelle incidence angle [deg]
    D_cl_1 = -0.3
    if nacelle_on_top_of_wing:
        D_cl_1 = 0.2
    D_cl_2 = 0.056 * i_n
    wing_nac_interference = 0.036 * (c_nac * b_nac / S) * (D_cl_1 + D_cl_2)**2

    SHP =  # shaft horse power
    D_prop =  # propeller diameter
    wind_milling = 33 / (0.5 * density * speed**2 * PA_TO_LBSpFT2 * S) * SHP / (speed * MpS_TO_FpS)
    if n_engine_operative != n_eng:
        wind_milling_inoperative = 0.00125 * ac.engine.eta_prop * D_prop**2 / S
    propulsion = n_eng * (wind_milling + wing_nac_interference + isolated_pylon + isolated_nac) + (n_eng - n_engine_operative) * wind_milling_inoperative
    
    # Flap:
    if flight_condition != 'cruise':
        y_start_flap = 
        y_end_flap = 
        cf_c =  # flap chord length / chord length
        flap_deflection =  # degrees
        if flap_type = 'split':  # 'split' or 'plain' or 'slotted' 'fowler' or 'krueger'
            fd = closest_value(t_c_max, values=[10, 20, 30])
            D_CD_flap_stuff = interp_value(pd.read_csv(f'lookups/t_c_0.{fd}.csv'), cf_c, f'cf/c ({flight_condition})', f'dCdp ({flight_condition})', log_x=False)
        elif flap_type == 'plain':
            fd = closest_value(flap_deflection, values=[15, 60])
            D_CD_flap_stuff = interp_value(pd.read_csv(f'lookups/d_f_{fd}.csv'), cf_c, 'cf/c', 'dCdp', log_x=False)
        elif flap_type == 'slotted':
            fd = closest_value(cf_c, values=[0.1, 0.2, 0.3])
            D_CD_flap_stuff = interp_value(pd.read_csv('lookups/cf_c_comb2.csv'), flap_deflection, 'df', f'dCdp (cf={fd}0)', log_x=False)
        elif flap_type == 'fowler':
            fd = closest_value(cf_c, values=[0.1, 0.2, 0.3, 0.4])
            D_CD_flap_stuff = interp_value(pd.read_csv('lookups/roskam_p6_fig_4_48.csv'), flap_deflection, f'df(cf/c={fd})', f'dCdp(cf/c={fd})', log_x=False)
        elif flap_type == 'krueger':
            D_CD_flap_stuff = wing * (cf_c * np.cos(np.deg2rad(flap_deflection)) + 1)
        D_CD_flap_stuff =   # See eqn 4.71 Roskam VI
        flap_profile = D_CD_flap_stuff * np.cos(np.deg2rad(sweep_c_4_deg)) * S_wf(y_start_flap, y_end_flap, taper, c_r, b) / S

        b_fi_b = y_start_flap * 2 / b
        b_fo_b = y_end_flap * 2 / b
        K = 
        Delta_CL_max_flapped = 
        induced_flap = K**2 * Delta_CL_max_flapped**2 * np.cos(np.deg2rad(sweep_c_4_deg))

        interference_flap = flap_profile * flap_interference_factor(flap_type)

        flaps = interference_flap + flap_profile + induced_flap
    else: 
        flaps = 0

    # Landing gear
    w_tire =  # tire width NOTE: check all lg dictionary names
    d_tire = ac.landing_gear.selected_mlg_tire['Tire Radius (In)'] * 2.54 / 100  # tire diameter
    w_strut = 
    l_strut = np.abs(ac.landing_gear.height_mlg) - ac.landing_gear.selected_mlg_tire['Tire Radius (In)'] * 2.54 / 100
    m = (w_tire * d_tire + l_strut * w_strut) / ((w_tire + w_strut) * (l_strut + 0.5 * d_tire))
    landing_gear = ((w_tire + w_strut) * (l_strut + 0.5 * d_tire)) * 0.04955 * np.exp(5.615 * m)

    # Miscelaneous (+5%)
    C_D0 = (wing + fuselage + ht + vt + propulsion + flaps + landing_gear) * 1.05

    return C_D0

def C_D_L(ac:Aircraft, 
          flap_deflection: float, # deg
          flight_condition: str = 'cruise' # 'cruise' or 'landing' or 'take-off'
          ):
    S = ac.wing.area
    if flight_condition == 'cruise':
        temp_shift = 0
        alt = ac.requirements.cruise['cr_altitude'] * FT_TO_M
        speed = ac.requirements.cruise['cr_speed'] * KTS_TO_MS
        Atm = Atmosphere(alt, temp_shift)
        density = Atm.density
        mass_frac = ac.requirements.cruise['cr_mass_frac']
        C_L = mass_frac * ac.weights.m_takeoff * g / (0.5 * density * speed**2 * S)
    if flight_condition == 'take-off':
        C_L = 
        temp_shift = ac.requirements.take_off['to_temp_shift']
        alt = ac.requirements.take_off['to_altitude'] * FT_TO_M
        Atm = Atmosphere(alt, temp_shift)
        density = Atm.density
        mass_frac = ac.requirements.cruise['to_mass_frac']
        speed = np.sqrt(mass_frac * ac.weights.m_takeoff / (0.5 * density * C_L * S))
    if flight_condition == 'take-off':
        temp_shift = ac.requirements.landing['la_temp_shift']
        alt = ac.requirements.landing['la_altitude'] * FT_TO_M
        Atm = Atmosphere(alt, temp_shift)
        density = Atm.density
        mass_frac = ac.requirements.landing['la_mass_frac']
        speed = 1.3 * ac.requirements.general['stall_speed'] * KTS_TO_MS
        C_L = mass_frac * ac.weights.m_takeoff * g / (0.5 * density * speed**2 * S)
    
    # General
    A = ac.wing.aspect_ratio
    A_eff = A + ... # NOTE: add wing tip effect here
    c_r = 
    b = ac.wing.span
    taper = ac.wing.taper_ratio
    if ac.wing.sweep == 0:
        e = 1.78 * (1 - 0.045 * A_eff**0.68) - 0.64
    else:
        e = 4.61 * (1 - 0.045 * A_eff**0.68) * (np.cos(np.deg2rad(LE_sweep_deg(ac.wing.sweep, c_r, b, taper))))**0.15 - 3.1

    # Wing

    # Fuselage

    # Empennage

    # Nacelle/pylon

    # Flap
    if flap_deflection != 0:
        e += 0.0046 * flap_deflection

    K = 1 / (np.pi * A_eff * e)
    CDi = C_L**2 * K
    tip_twist =  # degrees
    if tip_twist != 0:
        CDi += 0.00004 * 2 / 3 * tip_twist
    return CDi, e, K
