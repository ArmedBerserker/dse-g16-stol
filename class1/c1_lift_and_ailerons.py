# Fix path FIRST, before any local imports
import sys
import os

from c2_m import D_Cl_max, S_wf, sweep_at_x_c_deg, chord_at_y_span

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


def lift_slope(ac: Aircraft, flight_condition): # This lift slope is also already defined in c2_m
    V = ac.requirements.cruise['cr_speed'] * KTS_TO_MS
    atm = Atmosphere(ac.requirements.cruise['cr_altitude'] * FT_TO_M, delta_T=0)
    M = V / np.sqrt(1.4 * 287 * atm.temp)
    beta = np.sqrt(1 - M**2)                                                        #We also have a beta function defined in c2_m, use that??
    CLa = 2 * np.pi * ac.wing.aspect_ratio / (2 + np.sqrt(4 + (ac.wing.aspect_ratio * beta / 0.95)**2 * (1 + ac.wing.sweep_c_2_deg**2 / beta**2)))
    return CLa

def high_AR_method(ac: Aircraft, flight_condition,):      #make this flexible to calculate cl at various airspeeds, thus landing and take-off?
    V = ac.requirements.cruise['cr_speed'] * KTS_TO_MS
    atm = Atmosphere(ac.requirements.cruise['cr_altitude'] * FT_TO_M, delta_T=0)
    M = V / np.sqrt(1.4 * 287 * atm.temp)
    CL_l_ratio = ...                                                    #Naomi - ratio from figure on slide 22
    if M > 0.2:                                                         # needs to be local M, thus again dependt on what airspeed
        CL = CL_l_ratio * ... + DCLmax               #Naomi - CL-l ratio is from figure on slide 22, ... needs to be airfoil clmax at M=0.2, DCLmax function from slide 23, but if we only do take-off and landing then ignore this
    else:
        CL = CL_l_ratio * ...                            #Naomi - ... needs to be airfoil clmax at M=0.2
    return CL

def low_AR_method(ac: Aircraft, DY):
    CL_maxbase = ...                        #Naomi - from figure on slide 26
    DCLmax = ...                            #Naomi - from figure on slide 27
    CL = CL_maxbase * DCLmax
    return CL

def CLmax(ac: Aircraft):            #CL max of clean wing, needs to be made usable a multiple different speeds
    if ac.wing.airfoil_name == "NACA 23012":
        DY = 26 * ac.wing.t_c_max          #sharpness factor, #t_c_max not defined anywhere
        df = pd.read_csv('lookups/Lift/C1.csv')
        C1 = np.interp(ac.wing.taper_ratio, df['taper'], df['C1'])
        if ac.wing.aspect_ratio > (4 / ((C1 + 1) * np.cos(np.deg2rad(ac.wing.sweep_LE_deg)))):
            CL_max_w = high_AR_method()                 #functions inputs
        else:
            CL_max_w = low_AR_method()                  #function inputs
    else: 
        print('Airfoil changed, so redefine sharpness factor to proceed.')
        break
    return CL_max_w


######## Flaps #########

#### TO
def flap_dc_cf_TO(flap_type: str, # 'plain' or 'slotted' or 'fowler' or 'double slotted'or 'triple slotted' or 'fowler'
             ):
    flap_values = {
        'plain': 0.15,
        'slotted': 0.2,
        'fowler': 0.5,
        'double slotted': 0.6,
        'triple slotted': 0.6,
    }

    return flap_values.get(flap_type)

def cdash_c_TO(ac: Aircraft, flap_type: str):
    if flap_type == "slat":
        cdash_c = 1.05
    else:
        cdash_cf = flap_dc_cf_TO(flap_type) + (1 / ac.hld_and_ailerons.flaps['cf_c'])
        cdash_c = cdash_cf * ac.hld_and_ailerons.flaps['cf_c']
    return cdash_c


def D_CL_max_TO(ac: Aircraft, flap_type: str, y_start_f, y_end_f, taper, c_r, b, LE_sweep, x_c_hinge):         #Naomi x_c_hinge needs to be defined somewhere
    return 0.8 * 0.9 * D_Cl_max(flap_type) * (S_wf(y_start_f, y_end_f, taper, c_r, b)/ac.wing['area']) * np.cos(sweep_at_x_c_deg(LE_sweep, c_r, b, taper, x_c_hinge))             # Naomi fill in functions arguements


##### Landing



def flap_dc_cf_LA(flap_type: str, # 'plain' or 'slotted' or 'fowler' or 'double slotted'or 'triple slotted' or 'fowler'
             ):
    flap_values = {
        'plain': 0.45,
        'slotted': 0.3,
        'fowler': 0.6,
        'double slotted': 0.85,
        'triple slotted': 0.85,
    }

    return flap_values.get(flap_type)

def cdash_c_LA(flap_type: str, cf_c):
    if flap_type == "slat":
        cdash_c = 1.05
    else:
        cdash_cf = flap_dc_cf_LA(flap_type) + (1/cf_c)
        cdash_c = cdash_cf * cf_c
    return cdash_c


def D_CL_max_LA(ac: Aircraft, flap_type: str, y_start_f, y_end_f, taper, c_r, b, LE_sweep, x_c_hinge):         #Naomi x_c_hinge needs to be defined somewhere
    return 0.9 * D_Cl_max(flap_type) * (S_wf(y_start_f, y_end_f, taper, c_r, b)/ac.wing['area']) * np.cos(sweep_at_x_c_deg(LE_sweep, c_r, b, taper, x_c_hinge))             # Naomi fill in functions arguements


# Now how to check if it meets the CLmax depands
#make it iterate or so?

######## Aileron #############

# define aileron geom
# y_start_a
# y_end_a
# cf_c
# alieron max deflect


def aileron_control_derivative(ac: Aircraft,
    Cl_alpha,
    tau,
                               c_r,
                               taper,
                               b,
    n_points=1000
):
    """
    Compute Cl_delta_a using a DATCOM-style integration.

    Parameters
    ----------
    CL_alpha : float
        Wing lift curve slope [1/rad]

    y1 : float
        Inboard aileron span station [m]

    y2 : float
        Outboard aileron span station [m]

    chord_function : function
        Function c(y) returning local chord length
    """

    y = np.linspace(ac.hld_and_ailerons.ailerons['y_aileron_in'], ac.hld_and_ailerons.ailerons['y_aileron_out'], n_points)

    c = chord_at_y_span(c_r, taper, y, b)

    integrand = c * y

    integral = np.trapz(integrand, y)

    Cl_da = (2 * Cl_alpha * tau / (ac.wing['area'] * ac.wing['span'])) * integral

    return Cl_da

def roll_damping_derivative(ac: Aircraft,
    Cl_alpha,
                            cd0,
                               c_r,
                               taper,
                               b,
    n_points=1000
):

    y = np.linspace(ac.hld_and_ailerons.ailerons['y_aileron_in'], ac.hld_and_ailerons.ailerons['y_aileron_out'], n_points)

    c = chord_at_y_span(c_r, taper, y, b)

    integrand = c * y **2

    integral = np.trapz(integrand, y)

    Cl_p = - ((4 * (Cl_alpha * cd0)**(b/2)) / (ac.wing['area'] * ac.wing['span']**2)) * integral

    return Cl_p


def AC_roll_rate(ac: Aircraft,):
    P = - (aileron_control_derivative()/roll_damping_derivative()) * ... * (2 * ... /ac.wing['span'])
