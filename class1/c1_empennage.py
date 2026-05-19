import sys
import os
import numpy as np

# Fix path FIRST, before any local imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from classes.aircraft_2 import loader, Aircraft, Requirements, Mission, Weights, Wing, Fuselage, Engine
from lookups.consts import *
from class1.prelim_drag import *

def size_empennage_planform(ac: Aircraft, l_h: float = None, l_v: float = None):
    """ 
    Empennage planform sizing:
    Constants:
        - Volume coefficients (V_h, V_v)
        - Aspect ratios (A_h, A_v)
        - Taper ratios (lambda_h, lambda_v)
        - Sweep angles 
    Inputs:
        - Moment arms (l_h, l_v) [Optional overrides via function arguments]
        - From ac.wing: Area (S_w), Span (b_w), MAC (MAC_w)
    Outputs: (updating ac.htail and ac.vtail objects)
        - Surface Areas (S_h, S_v)
        - Spans (b_h, b_v)
        - Root and tip chords
        - MAC and its spanwise position
    """
    w = ac.wing
    ht = ac.htail
    vt = ac.vtail
    
    # Update moment arms if they are provided as inputs
    if l_h is not None:
        ht.moment_arm = l_h
    if l_v is not None:
        vt.moment_arm = l_v
    
    # ---------------------------------------------------------
    # HORIZONTAL STABILIZER SIZING
    # ---------------------------------------------------------
    # 1. Calculate new area based on volume coefficient and input moment arm
    ht.area = (ht.volume_coefficient * w.area * w.MAC) / ht.moment_arm
    
    # 2. Update geometry dependent on the new area
    ht.span = np.sqrt(ht.aspect_ratio * ht.area)
    ht.c_root = 2 * ht.area / (ht.span * (1 + ht.taper_ratio))
    ht.c_tip = ht.c_root * ht.taper_ratio
    
    ht.MAC = (2 / 3) * ht.c_root * (1 + ht.taper_ratio + ht.taper_ratio**2) / (1 + ht.taper_ratio)
    ht.y_MAC = (ht.c_root - ht.MAC) / (ht.c_root * (1 - ht.taper_ratio)) * (ht.span / 2)

    # ---------------------------------------------------------
    # VERTICAL STABILIZER SIZING
    # ---------------------------------------------------------
    # 1. Calculate new area based on volume coefficient and input moment arm
    vt.area = (vt.volume_coefficient * w.area * w.span) / vt.moment_arm
    
    # 2. Update geometry dependent on the new area
    vt.span = np.sqrt(vt.aspect_ratio * vt.area)
    vt.c_root = 2 * vt.area / (vt.span * (1 + vt.taper_ratio))
    vt.c_tip = vt.c_root * vt.taper_ratio
    
    vt.MAC = (2 / 3) * vt.c_root * (1 + vt.taper_ratio + vt.taper_ratio**2) / (1 + vt.taper_ratio)
    # y_MAC for a vertical tail (half-span surface) spans fully from root to tip
    vt.y_MAC = (vt.c_root - vt.MAC) / (vt.c_root * (1 - vt.taper_ratio)) * vt.span