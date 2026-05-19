import sys
import os
import numpy as np
import yaml

# Fix path FIRST, before any local imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from classes.aircraft_2 import loader, Aircraft, Requirements, Mission, Weights, Wing, Fuselage, Engine
from lookups.consts import *
from class1.prelim_drag import *

def size_empennage_planform(ac: Aircraft, config_path: str = 'empennage_config.yaml'):
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
    
    # Load parameters from the YAML file
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
        
    w = ac.wing
    ht = ac.htail
    vt = ac.vtail
    
    # Extract positions to calculate the moment arms
    x_cgaft = config['cg_and_positioning']['x_cgaft']
    x_h = config['cg_and_positioning']['x_h']
    x_v = config['cg_and_positioning']['x_v']
    
    # Calculate moment arms dynamically
    ht.moment_arm = x_h - x_cgaft
    vt.moment_arm = x_v - x_cgaft
    
    # Extract fixed constants from YAML
    ht.volume_coefficient = config['horizontal_tail']['volume_coefficient']
    ht.aspect_ratio = config['horizontal_tail']['aspect_ratio']
    ht.taper_ratio = config['horizontal_tail']['taper_ratio']
    
    vt.volume_coefficient = config['vertical_tail']['volume_coefficient']
    vt.aspect_ratio = config['vertical_tail']['aspect_ratio']
    vt.taper_ratio = config['vertical_tail']['taper_ratio']

    # ---------------------------------------------------------
    # HORIZONTAL STABILIZER SIZING
    # ---------------------------------------------------------
    # 1. Calculate new area based on volume coefficient and calculated moment arm
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
    # 1. Calculate new area based on volume coefficient and calculated moment arm
    vt.area = (vt.volume_coefficient * w.area * w.span) / vt.moment_arm
    
    # 2. Update geometry dependent on the new area
    vt.span = np.sqrt(vt.aspect_ratio * vt.area)
    vt.c_root = 2 * vt.area / (vt.span * (1 + vt.taper_ratio))
    vt.c_tip = vt.c_root * vt.taper_ratio
    
    vt.MAC = (2 / 3) * vt.c_root * (1 + vt.taper_ratio + vt.taper_ratio**2) / (1 + vt.taper_ratio)
    vt.y_MAC = (vt.c_root - vt.MAC) / (vt.c_root * (1 - vt.taper_ratio)) * vt.span