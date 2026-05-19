import sys
import os
import numpy as np
import yaml

# Fix path FIRST, before any local imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from classes.aircraft_2 import loader, Aircraft, Requirements, Mission, Weights, Wing, Fuselage, Engine
from lookups.consts import *
from class1.prelim_drag import *

def size_empennage_planform(ac: Aircraft):
    """ 
    Empennage planform sizing:
    Inputs (from YAML config):
        - Positioning: x_cgaft, x_h, x_v
        - Fixed Constants: Volume coefficients, Aspect ratios, Taper ratios
    Inputs (from ac.wing):
        - Area (S_w), Span (b_w), MAC (MAC_w)
    Outputs: (updating ac.htail and ac.vtail objects)
        - Moment arms (l_h, l_v)
        - Surface Areas (S_h, S_v)
        - Spans (b_h, b_v)
        - Root and tip chords
        - MAC and its spanwise position
    """
        
    w = ac.wing
    ht = ac.empennage.horizontal_tail
    vt = ac.empennage.vertical_tail

    S = w.area
    b = w.span
    
    # Extract positions to calculate the moment arms
    l_f = ac.fuselage.length
    x_cgaft = ac.weights.x_cg_aft
    x_h = ht['x_h_frac_lf'] * l_f
    x_v = vt['x_v_frac_lf'] * l_f
    
    # Calculate moment arms dynamically
    ht_moment_arm = x_h - x_cgaft
    vt_moment_arm = x_v - x_cgaft
    
    # Extract fixed constants from YAML
    V_h = ht['volume_coefficient']
    A_h = ht['aspect_ratio']
    taper_h = ht['taper_ratio']
    V_v = vt['volume_coefficient']
    A_v = vt['aspect_ratio']
    taper_v = vt['taper_ratio']

    # ---------------------------------------------------------
    # HORIZONTAL STABILIZER SIZING
    # ---------------------------------------------------------
    # 1. Calculate new area based on volume coefficient and calculated moment arm
    S_h = (V_h * S * w.MAC) / ht_moment_arm

    # 2. Update geometry dependent on the new area
    b_h = np.sqrt(A_h * S_h)
    c_r_h = 2 * ht.area / (ht.span * (1 + taper_h))
    c_t_h = c_r_t * taper_h
    
    MAC_h = (2 / 3) * c_r_h * (1 + taper_h + taper_h**2) / (1 + taper_h)
    y_MAC_h = (c_r_h - MAC_h) / (c_r_h * (1 - taper_h)) * (b_h / 2)

    ht['area'] = S_h
    ht['b_h'] = b_h
    ht['c_r_h'] = c_r_h
    ht['c_t_h'] = c_t_h
    ht['MAC_h'] = MAC_h
    ht['y_MAC_h'] = y_MAC_h

    # ---------------------------------------------------------
    # VERTICAL STABILIZER SIZING
    # ---------------------------------------------------------
    # 1. Calculate new area based on volume coefficient and calculated moment arm
    S_v = (V_v * S * b) / vt_moment_arm
    
    # 2. Update geometry dependent on the new area
    b_v = np.sqrt(A_v * S_v)
    c_r_v = 2 * S_v / (b_v * (1 + taper_v))
    c_t_v = c_r_v * taper_v
    
    MAC_v = (2 / 3) * c_r_v * (1 + taper_v + taper_v**2) / (1 + taper_v)
    y_MAC_v = (c_r_v - MAC_v) / (c_r_v * (1 - taper_v)) * b_v

    vt['area'] = S_v
    vt['b_v'] = b_v
    vt['c_r_v'] = c_r_v
    vt['c_t_v'] = c_t_v
    vt['MAC_v'] = MAC_v
    vt['y_MAC_v'] = y_MAC_v