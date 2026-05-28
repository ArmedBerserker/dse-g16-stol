from classes.aircraft_2 import Aircraft
from c2_m import Snet, LE_sweep_deg, sweep_at_x_c_deg, lift_slope, S_wf, closest_value, chord_at_y_span
from lookups.consts import *
import pandas as pd
import numpy as np
from scipy.interpolate import RegularGridInterpolator, interp1d
from scipy.optimize import brentq
from classes.isa import Atmosphere
import matplotlib.pyplot as plt

# Data and interpolation functions
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

# Intermediate variable functions
def mu_air(T):
    return 1.81 * 1e-5 * (T / 293.15)**1.5 * (293.15 + 110.4) / (T + 110.4)

def R_wf(M, R_N_fus):
    M = closest_value(M, values=[0.25, 0.4])
    if M == 0.25:
        return interp_value(pd.read_csv('lookups/roskam_p6_fig_4_1_rwf.csv'), R_N_fus, x_col='Fuselage Reynolds Number (M = 0.25)', y_col='Wing-Fuse Interference Factor (M = 0.25)', log_x=True)
    else:
        return interp_value(pd.read_csv('lookups/roskam_p6_fig_4_1_rwf_2.csv'), R_N_fus, x_col='Fuselage Reynolds Number (M = 0.4)', y_col='Wing-Fuse Interference Factor (M = 0.4)', log_x=True)

def R_LS(sweep_t_c_max_deg):
    return interp_value(pd.read_csv('lookups/roskam_p6_fig_4_2_rls.csv'), np.cos(np.deg2rad(sweep_t_c_max_deg)), 'cos(quarter chord)', 'lifting surface correction', log_x=False)

def R_N(rho, V, l, mu):
    return rho * V * l / mu

def C_f(R_N, M):  # DATCOM fig 4.1.5.1-26
    x = np.log10(R_N)
    CF3 = 3.92725e-6 * x**5 -1.30370e-4 * x**4 -1.65388e-3 * x**3 -9.59519e-3 * x**2 +2.18366e-2 * x 
    CF0 = 4.12963e-6 * x**5 -1.36204e-4 * x**4 + 1.71620e-3 * x**3 -9.88935e-3 * x**2 +2.23641e-2 * x
    if M > 0.3:
        return CF3
    else: 
        return CF3 + (CF0 - CF3) * (1 - M / 0.3)

def L_dash(x_c_t_c_max):
    if x_c_t_c_max < 0.3:
        return 2.0
    else:
        return 1.2

def S_wet_w(S_exp, surface_type: str = 'wing' # 'wing' or other like 'empennage' or 'pylon'
            ):
    if surface_type == 'wing':
        return 2 * 1.07 * S_exp
    else:
        return 2 * 1.05 * S_exp
    
def S_wet_fus(l_nosecone, l_tot, l_tailcone, d_max):
    l2 = l_tot - l_nosecone - l_tailcone
    return np.pi * d_max / 4 * (1 / (3 * l_nosecone**2) * ((4 * l_nosecone**2 + d_max**2 / 4)**1.5 - d_max**3 / 8) - d_max + 4 * l2 + 2 * np.sqrt(l_tailcone**2 + d_max**2 / 4))

def Snet(S, b_f, taper, b, c_r):
        c_fus_int = chord_at_y_span(c_r, taper, b_f/2, b)
        return S - (c_r + c_fus_int) * b_f / 2

def exposed_wing_mgc(S, b, taper, c_r, b_f):
    c_fus_int = chord_at_y_span(c_r, taper, b_f/2, b)
    S_net = S - (c_r + c_fus_int) * b_f / 2
    return S_net / (b - b_f)

def exposed_vt_mgc(S_fus_base, l_f, l_tc, x_le_vt, S_vt, b_vt, taper_vt, c_r_vt, h_f):
    d_min = np.sqrt(4 * S_fus_base / np.pi)
    d_fus_local = h_f - (x_le_vt - (l_f - l_tc)) / l_tc * (h_f - d_min)
    c_fus_int = chord_at_y_span(c_r_vt, taper_vt, d_fus_local/2, b_vt)
    S_net = S_vt - (c_r_vt + c_fus_int) / 2 * b_vt
    return S_net / (b_vt - d_fus_local / 2)

# Component drag eqns

def CD0_wing(M, R_N_fus, sweep_t_c_max_deg, R_N_w, t_c_max, x_c_t_c_max, S_exp, surface_type, S):
    Rwf = R_wf(M, R_N_fus)
    RLS = R_LS(sweep_t_c_max_deg)
    Cfw = C_f(R_N_w, M)
    Ldash = L_dash(x_c_t_c_max)
    Swet = S_wet_w(S_exp, surface_type)
    return Rwf * RLS * Cfw * (1 + Ldash * t_c_max + 100 * t_c_max**4) * Swet / S

def CD0_fusselage()