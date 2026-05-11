from classes.aircraft_2 import Aircraft
from lookups.consts import *
import pandas as pd
import numpy as np

def W_to(w_oe, w_f, w_pl, w_tfo, w_crew):
    return sum(w_oe, w_f, w_pl, w_tfo, w_crew)

def W_oe(w_struc, w_pwr, w_feq):
    return sum(w_struc, w_pwr, w_feq)

def W_to_new(w_oe, # from class II
             w_pl, # known
             m_ff, # from class I
             m_res, # from class I
             m_tfo, # from class I
             w_crew = 0.0, # included in PL probably
             ):
    return (w_oe + w_pl + w_crew) / (m_ff * (1 + m_res) - m_res - m_tfo)

def w_struc(w_wing, w_emp, w_fus, w_nac, w_gear):
    return sum(w_wing, w_emp, w_fus, w_nac, w_gear)

def Cessna_method_general_av(ac: Aircraft,
                             wing_type: str = 'cantilever'  # 'cantilever' or 'strut braced'
                             ):
    ''' Wing weight:
        - Includes: wing tip fairing, control surfaces
        - Excludes: fuel tanks, wing/fuselage spar carry-through structure, effect of sweep
        - Max thickness = 0.18c

        Empennage:
        - HT: no sweep included
        '''
    W_to = ac.weights.m_takeoff / LBS_TO_KG
    S = ac.wing.area * M2_TO_F2
    n_ult = ...  # NOTE: add later
    A = ac.wing.aspect_ratio

    # Wing
    if wing_type == 'cantilever':
        W_wing = 0.04674 * (W_to ** 0.397) * (S ** 0.36) * (n_ult ** 0.397) * (A ** 1.712)

    else:
        W_wing = 0.002933 * (S ** 1.018) * (A ** 2.473) * (n_ult ** 0.611)

    # Empennage NOTE: fill in values
    S_h = x * M2_TO_F2
    S_v = x * M2_TO_F2
    A_h = x
    W_h = (3.184 * (W_to ** 0.887) * (S_h ** 0.101) * (A_h ** 0.138)) / (174.04 * (t_r_h **0.223))
    W_v = (1.68 * (W_to ** 0.567) * (S_v ** 1.249) * (A_v ** 0.482)) / (639.95 * (t_r_v ** 0.747) * (np.cos(np.rad2deg(vt_sweep_c_4_deg)) ** 0.882))


def USAF_method_general_av(ac: Aircraft):
    
    W_to = ac.weights.m_takeoff / LBS_TO_KG
    S = ac.wing.area * M2_TO_F2
    n_ult = ...  # NOTE: add later
    A = ac.wing.aspect_ratio
    t_c_max = ...
    sweep_c_4_deg = ...
    taper = ...
    V_h = ...  # NOTE: add later: maximum level speed at sealevel in kts

    W_wing = 96.948 * (((W_to * n_ult * 1e-5) ** 0.65) * ((A / np.cos(np.rad2deg(sweep_c_4_deg))) ** 0.57) * ((S / 100) ** 0.61) * (((1 + taper) / (2 * t_c_max)) ** 0.36) * ((1 + V_h / 500) ** 0.5)) ** 0.993

    W_emp = 
    
def Torenbeek_method_general_av(ac: Aircraft):
    
    W_to = ac.weights.m_takeoff / LBS_TO_KG
    S = ac.wing.area * M2_TO_F2
    n_ult = ...  # NOTE: add later
    A = ac.wing.aspect_ratio
    b = ac.wing.span / FT_TO_M
    c_r = ... / FT_TO_M
    t_c_max = ...
    t_r = t_c_max * c_r  # max thickness of wing root chord in ft
    sweep_c_4_deg = ...
    taper = ...
    sweep_le_rad = np.arctan(np.tan(np.rad2deg(sweep_c_4_deg)) + c_r / 2 / b * (1 + taper))
    sweep_c_2_rad = np.arctan(np.tan(sweep_le_rad) - c_r / b * (1 + taper))

    W_wing = 0.00125 * W_to * ((b / np.cos(sweep_c_2_rad)) ** 0.75) * (1 + (6.3 * np.cos(sweep_c_2_rad) / b) ** 0.5) * (n_ult ** 0.55) * (b * S / (t_r * W_to * np.cos(sweep_c_2_rad)))**0.3
