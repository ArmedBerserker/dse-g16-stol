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
from prelim_drag import *
import pandas as pd
import numpy as np

''' Key points to pay attention to:
    - W_F (final weight): is this calculated correctly with the fuel used, is it set up and 
        loaded correctly to work for all fuel types?
    - Are conversions done right?
    - Are aircraft engine types checked for right?'''

def stall_speed_matching(ac : Aircraft,  # Change units
        type_to_use : str = "Single Engine Propeller Driven",
        friction_source : str = 'lookups/skin_fric.csv',
        s_wet_source : str = 'lookups/s_wets.csv',
        W_P_or_T_W: np.ndarray = np.arange(0,10000,5)) -> np.ndarray:
    ''' Returns list of arrays (W/P or T/W and W_S)'''
    V_s = (ac.requirements.general.stall_speed)*KTS_TO_MS
    rho =  # at SL ISA+20
    CL_max = 
    W_S = (CL_max*0.5*rho*V_s**2)*np.ones_like(W_P_or_T_W)
    return (W_P_or_T_W, W_S)

def takeoff_dist_matching(ac : Aircraft,  # Change units
        type_to_use : str = "Single Engine Propeller Driven",
        friction_source : str = 'lookups/skin_fric.csv',
        s_wet_source : str = 'lookups/s_wets.csv',
        W_S: np.ndarray = np.arange(0,10000,5)) -> float:
    ''' Returns list of arrays (W/P or T/W and W_S)'''

    balance_fl =  # True or False dep on yaml
    propeller =  # True or False dep on yaml

    alt = ac.requirements.take_off.to_altitude
    temp = ac.requirements.take_off.to_temperature
    sigma =  # Density ratio
    C_Lmax_TO = 

    if propeller == True:

        if balance_fl == True:
            S_TO = ac.requirements.take_off.to_distance
            # Add check for prop or jet
            ''' Eqn 3.6 Roskam'''
            x = (-8.134+np.sqrt(8.134**2+S_TO*4*0.0149))
        else:
            S_TOG = ac.requirements.take_off.to_distance
            x = (-4.9+np.sqrt(4.9**2+4*0.009*S_TOG))
            # S_TO = 1.66*S_TOG

        W_P = x*sigma*C_Lmax_TO/W_S
        return (W_P, W_S)

    else:
        S_TOFL = ac.requirements.take_off.to_distance
        T_W = 37.5*W_S/(S_TOFL*sigma*C_Lmax_TO)
        return (T_W, W_S)


def landing_dist_matching(ac : Aircraft,
        type_to_use : str = "Single Engine Propeller Driven",
        friction_source : str = 'lookups/skin_fric.csv',
        s_wet_source : str = 'lookups/s_wets.csv',
        W_P_or_T_W: np.ndarray = np.arange(0,10000,5)) -> float:
    ''' Returns list of arrays (W/P or T/W and W_S)'''
    # W_S = np.ones_like(W_P_or_T_W)

    W_L = ac.weights.m_takeoff-ac.weights.m_fuel  # NOTE: Change to fraction of fuel burned not all
    # W_S_LD = W_S/ac.weights.m_takeoff*W_L
    C_L_max_LD =  
    # V_s = (ac.requirements.general.stall_speed)*KTS_TO_MS
    # V_a = 1.3*V_s  # Check if this is V_s_L (stall speed in landing config)

    balance_fl =  # True or False dep on yaml

    alt = ac.requirements.landing.la_altitude
    temp = ac.requirements.landing.la_temperature
    # sigma =  # Density ratio
    rho =  
    C_Lmax_LD =  # 

    if balance_fl == True:
        S_L = ac.requirements.landing.la_distance
        # S_LG = S_L/1.938
        V_S_L = np.sqrt(S_L/0.5136)

    else: 
        S_LG = ac.requirements.landing.la_distance
        V_S_L = np.sqrt(S_LG/0.265)

    W_S_LD = (V_S_L**2/(1/2*rho*C_L_max_LD))*np.ones_like(W_P_or_T_W)
    W_S = W_S_LD*ac.weights.m_takeoff/W_L

    return (W_P_or_T_W, W_S)

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


