# Fix path FIRST, before any local imports
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from classes.aircraft_2 import loader, Aircraft, Requirements, Mission, Weights, Wing, Fuselage, Engine
from lookups.consts import *
from class1.prelim_drag import *
# import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.cm as cm
# import matplotlib.lines as mlines
# from classes.isa import Atmosphere
# from class1.matching_diagram import plot_matching_and_select_design_point
# from class1 import c1_m
# from pathlib import Path
from aerodynamics.airfoils import AIRFOIL_DB


def sweep_at_x_c_deg(x_c, sweep_c_4_deg, c_r, taper, b):
    le_sweep = np.arctan(np.tan(np.deg2rad(sweep_c_4_deg)) + 0.5 * c_r / b * (1 - taper))
    sweep_at_x_c = np.arctan(np.tan(le_sweep) - x_c * 2 * c_r / b * (1 + taper))
    return np.rad2deg(sweep_at_x_c)

def size_wing_planform(ac: Aircraft):
    ''' Wing planform sizing:
    Constants:
        - taper ratio
        - wing incidence angle 
        - dihedral
        - twist
        - t/c
        - spar locations
        - airfoil
    Inputs:
        - A
        - S
        - taper ratio
        - sweep c/4 
        - Airfoil properties
        - airfoil
    Outputs: (updating ac object)
        - span
        - root and tip chord
        - LE c/2 and TE sweep
        - mac and its spanwise position
    '''
    w = ac.wing
    A = w.aspect_ratio
    S = w.area
    taper = w.taper_ratio
    sweep_c_4_deg = w.sweep
    b = np.sqrt(A * S)
    c_r = 2 * S / w.span / (1 + taper)
    c_t = c_r * taper
    le_sweep = sweep_at_x_c_deg(0, sweep_c_4_deg, w.c_root, taper, w.span)
    c_2_sweep = sweep_at_x_c_deg(0.5, sweep_c_4_deg, w.c_root, taper, w.span)
    te_sweep = sweep_at_x_c_deg(1, sweep_c_4_deg, c_r, taper, b)
    mac = 2 / 3 * c_r * (1 + taper + taper**2) / (1 + taper)
    y_mac = (c_r - mac) / (c_r * (1 - taper)) * b / 2
    airfoil_name = ac.wing.airfoil_name
    airfoil = next((a for a in AIRFOIL_DB if a.name == airfoil_name), None)

    # Update Aircraft
    if airfoil:
        w.t_c_max = airfoil.thickness_pct / 100
        w.camber_c = airfoil.camber_pct / 100
        w.cl_max_2d = airfoil.cl_max_2d
        w.cm_c4 = airfoil.cm_c4
        w.drag_bucket_cl = airfoil.drag_bucket_cl
        w.stall_type = airfoil.stall_type
        w.alpha_0L = airfoil.alpha_0L
    else:
        print(f'Error: airfoil {airfoil_name} not found in database')

    w.span = b
    w.c_root = c_r
    w.c_tip = c_t
    w.sweep_LE_deg = le_sweep
    w.sweep_c_2_deg = c_2_sweep
    w.sweep_TE_deg = te_sweep
    w.MAC = mac
    w.y_MAC = y_mac