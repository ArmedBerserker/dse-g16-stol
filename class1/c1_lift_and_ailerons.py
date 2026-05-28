# Fix path FIRST, before any local imports
import sys
import os
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

def lift_slope(ac: Aircraft, flight_condition):
    V = ac.requirements.cruise['cr_speed'] * KTS_TO_MS
    atm = Atmosphere(ac.requirements.cruise['cr_altitude'] * FT_TO_M, delta_T=0)
    M = V / np.sqrt(1.4 * 287 * atm.temp)
    beta = 1 - M**2
    CLa = 2 * np.pi * ac.wing.aspect_ratio / (2 + np.sqrt(4 + (ac.wing.aspect_ratio * beta / 0.95)**2 * (1 + ac.wing.sweep_c_2_deg**2 / beta**2)))
    return CLa

def high_AR_method(ac: Aircraft):
    CL = ...
    return CL

def low_AR_method(ac: Aircraft):
    CL = ...
    return CL

def CLmax(ac: Aircraft):
    if ac.wing.airfoil_name == "NACA 23012":
        DY = 26 * ac.wing.t_c_max
        df = pd.read_csv('lookups/Lift/C1.csv')
        C1 = np.interp(ac.wing.taper_ratio, df['taper'], df['C1'])
        if ac.wing.aspect_ratio > (4 / ((C1 + 1) * np.cos(np.deg2rad(ac.wing.sweep_LE_deg)))):
            ...
        else:
            ...
    else: 
        print('Airfoil changed, so redefine sharpness factor to proceed.')
        break