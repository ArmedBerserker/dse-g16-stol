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
        '''
    W_to = ac.weights.m_takeoff / LBS_TO_KG
    S = ac.wing.area * M2_TO_F2
    n_ult = ...  # NOTE: add later
    A = ac.wing.aspect_ratio
    if wing_type == 'cantilever':
        W_w = 0.04674 * (W_to ** 0.397) * (S ** 0.36) * (n_ult ** 0.397) * (A ** 1.712)

    else:
        W_w = 0.002933 * (S ** 1.018) * (A ** 2.473) * (n_ult ** 0.611)

def USAF_method_general_av(ac: Aircraft,
                             wing_type: str = 'cantilever'  # 'cantilever' or 'strut braced'
                             ):
    ''' Wing weight:
        - Includes: wing tip fairing, control surfaces
        - Excludes: fuel tanks, wing/fuselage spar carry-through structure, effect of sweep
        - Max thickness = 0.18c
        '''
    W_to = ac.weights.m_takeoff / LBS_TO_KG
    S = ac.wing.area * M2_TO_F2
    n_ult = ...  # NOTE: add later
    A = ac.wing.aspect_ratio
    if wing_type == 'cantilever':
        W_w = 0.04674 * (W_to ** 0.397) * (S ** 0.36) * (n_ult ** 0.397) * (A ** 1.712)

    else:
        W_w = 0.002933 * (S ** 1.018) * (A ** 2.473) * (n_ult ** 0.611)