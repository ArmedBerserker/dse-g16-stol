from classes.aircraft_2 import Aircraft
from c2_m import Snet, LE_sweep_deg, sweep_at_x_c_deg, lift_slope
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

def C_D0(ac: Aircraft, 
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
    l_dash = 
    t_c_max =  # at mean geometric chord
    x_c_t_c_max = 
    t_c_r = 
    t_c_t = 
    tau = t_c_t / t_c_r
    sweep_t_c_max_deg = sweep_at_x_c_deg(LE_sweep_deg(sweep_c_4_deg, c_r, b, taper), c_r, b, taper, x_c=x_c_t_c_max)
    S_wet_w = 2 * Snet(S, b_f, taper, b, c_r) * (1 + 0.25 * t_c_r * (1 + tau * taper) / (1 + taper))
    S = ac.wing.area

    wing = R_wf * R_LS * C_f_w * (1 + l_dash * t_c_max + 100 * t_c_max**4) * S_wet_w / S

    # Fuselage

    # Empennage

    # Nacelle/pylon

    # Flap

    # Landing gear

    # Canopy/windshield

    # Store(s)

    # Trim

    # Interference

    # Miscelaneous

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

    # Canopy/windshield

    # Store(s)

    # Trim

    # Interference

    # Miscelaneous
