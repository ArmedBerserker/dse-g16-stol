from classes.aircraft_2 import Aircraft
from c2_m import Snet, LE_sweep_deg, sweep_at_x_c_deg, lift_slope, S_wf
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
                 S_b_fus, S):
    R_N_fus = density * speed * l_f / mu
    R_wf = 
    C_f_fus = 
    d_b = np.sqrt(4 * S_b_fus / np.pi)
    d_f_eq = np.sqrt(4 * S_fus / np.pi)
    S_wet_fus =  # NOTE: add method
    C_D0_b_fus = R_wf * C_f_fus * (1 + 60 / (l_f / d_f_eq)**3 + 0.0025 * (l_f / d_f_eq)) * S_wet_fus / S
    C_D_b_fuse = (0.029 * (d_b / d_f_eq)**3 / (C_D0_b_fus * (S / S_fus))**0.5) * (S / S_fus)
    return C_D0_b_fus + C_D_b_fuse

def S_wet_wing(t_c_r, t_c_t, S, b_f, taper, b, c_r):
    tau = t_c_t / t_c_r
    return 2 * Snet(S, b_f, taper, b, c_r) * (1 + 0.25 * t_c_r * (1 + tau * taper) / (1 + taper))

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

def C_D0(ac: Aircraft, 
         n_engine_operative: int,
         flap_type: str, # 'split' or 'plain' or 'slotted' 'fowler' or 'krueger'
         flight_condition: str = 'cruise', # 'cruise' or 'landing' or 'take-off'
         nacelle_on_top_of_wing: bool = True, 
         ):
    if flight_condition == 'cruise':
        temp_shift = 0
        alt = ac.requirements.cruise['cr_altitude'] * FT_TO_M
        speed = ac.requirements.cruise['cr_speed'] * KTS_TO_MS
        Atm = Atmosphere(alt, temp_shift)
        temp = Atm.temp
        density = Atm.density
        
    # Wing
    mu = mu_air(temp)
    M = speed / np.sqrt(1.4 * 287 * temp)
    S = ac.wing.area
    b = ac.wing.span
    sweep_c_4_deg = ac.wing.sweep
    taper = ac.wing.taper_ratio
    c_r = 
    b_f = ac.fuselage.width
    l_f = ac.fuselage.length
    c_w_e = exposed_wing_mgc(S, b, taper, c_r, b_f)
    R_N_w = density * speed * c_w_e / mu
    R_N_fus = density * speed * l_f / mu
    R_wf = 
    R_LS = 
    C_f_w = 
    l_dash = 1.2
    t_c_max =  # at mean geometric chord
    x_c_t_c_max = 
    if x_c_t_c_max < 0.3:
        l_dash = 2.0
    t_c_r = 
    t_c_t = 
    sweep_t_c_max_deg = sweep_at_x_c_deg(LE_sweep_deg(sweep_c_4_deg, c_r, b, taper), c_r, b, taper, x_c=x_c_t_c_max)
    S = ac.wing.area

    wing = CD0_w(R_wf, R_LS, C_f_w, l_dash, t_c_max, S_wet_wing(t_c_r, t_c_t, S, b_f, taper, b, c_r), S)  # R_wf * R_LS * C_f_w * (1 + l_dash * t_c_max + 100 * t_c_max**4) * S_wet_w / S

    # Fuselage
    S_fus =  # max fuselage cross section area
    S_b_fus = 
    fuselage = CD0_fuselage(density, speed, l_f, mu, S_fus, S_b_fus, S)

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
    ht_sweep_t_c_max_deg = sweep_at_x_c_deg(LE_sweep_deg(ht_sweep_c_4_deg, c_r_ht, b_ht, taper_ht), c_r_ht, b_ht, taper_ht, x_c=ht_x_c_t_c_max)
    R_LS_ht = 
    b_f_ht =  # fuselage width at ht intersection position
    R_N_ht = density * speed * exposed_wing_mgc(S_ht, b_ht, taper_ht, c_r_ht, b_f_ht) / mu
    if not t_tail_condition:
        R_N_ht *= np.sqrt(0.85)
    C_f_ht = 
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
    R_LS_vt = 
    b_f_vt =  # fuselage width at vt intersection position
    R_N_vt = density * speed * exposed_wing_mgc(S_vt, b_vt, taper_vt, c_r_vt, b_f_vt) / mu
    C_f_vt = 
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
    isolated_nac = CD0_fuselage(density, speed, l_nac, mu, S_nac_max, S_b_nac, S)

    c_r_nac = 
    b_nac = 
    nac_t_c_max = 
    nac_x_c_t_c_max = 
    nac_sweep_t_c_max_deg = sweep_at_x_c_deg(LE_sweep_deg(sweep_c_4=0, c_r=c_r_nac, b=b_nac, taper_ratio=0), c_r_nac, b_nac, taper_ratio=0, x_c=nac_x_c_t_c_max)
    R_LS_nac = 
    R_N_nac = density * speed * c_r_nac / mu
    C_f_nac = 
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
        D_CD_flap_stuff =  # See eqn 4.71 Roskam VI
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
    w_tire =  # tire width
    d_tire =  # tire diameter
    l_strut = 
    w_strut = 
    m = (w_tire * d_tire + l_strut * w_strut) / ((w_tire + w_strut) * (l_strut + 0.5 * d_tire))
    landing_gear = ((w_tire + w_strut) * (l_strut + 0.5 * d_tire)) * 0.04955 * np.exp(5.615 * m)

    # Miscelaneous (+5%)
    C_D0 = (wing + fuselage + ht + vt + propulsion + flaps + landing_gear) * 1.05

    return C_D0

def C_D_L(ac:Aircraft, 
         flight_condition: str = 'cruise' # 'cruise' or 'landing' or 'take-off'
         ):
    
    if flight_condition == 'cruise':
        temp_shift = 0
        alt = ac.requirements.cruise['cr_altitude'] * FT_TO_M
        speed = ac.requirements.cruise['cr_speed'] * KTS_TO_MS
        Atm = Atmosphere(alt, temp_shift)
        temp = Atm.temp
        density = Atm.density
        
    # Wing

    # Fuselage

    # Empennage

    # Nacelle/pylon

    # Flap

    # Landing gear

    # Miscelaneous
