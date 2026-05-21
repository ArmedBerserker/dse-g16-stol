import matplotlib.pyplot as plt
import sys
import os
import numpy as np
import yaml

# Fix path FIRST, before any local imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from classes.aircraft_2 import loader, Aircraft, Requirements, Mission, Weights, Wing, Fuselage, Engine
from lookups.consts import *
from class1.prelim_drag import *

def loading_diagram_data(ac: Aircraft, tricyle_condition: bool):
    concepts = {

        "Taildragger Configuration": {
            "x": [2.90, 3.23, 3.19, 2.88, 2.90],
            "m": [1046.0, 1750.0, 1870.0, 1166.0, 1046.0],
            "labels": ["OEW", "OEW+WP", "OEW+WP+WF", "OEW+WF", ""]
        },

        "Tricycle Configuration": {
            "x": [2.85, 3.20, 3.16, 2.83, 2.85],
            "m": [1046.0, 1750.0, 1870.0, 1166.0, 1046.0],
            "labels": ["OEW", "OEW+WP", "OEW+WP+WF", "OEW+WF", ""]
        }
    }

    if tricyle_condition:
        concepts = concepts["Tricycle Configuration"]
    else:
        concepts = concepts["Taildragger Configuration"]
    return concepts

def class_I_loading_cgs(ac: Aircraft, tricycle_condition: bool, update_ac=False):
    concept_data = loading_diagram_data(ac, tricycle_condition)
    x = concept_data["x"]
    m = concept_data["m"]
    labels = concept_data["labels"]
    fwd_cg = np.min(x)
    aft_cg = np.max(x)
    if update_ac:
        ac.weights.x_cg_aft = aft_cg
        ac.weights.x_cg_fwd = fwd_cg
    else:
        return fwd_cg, aft_cg

    
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
    ht['l_h'] = ht_moment_arm
    
    # Extract fixed constants from YAML
    V_h = ht['volume_coefficient']
    A_h = ht['aspect_ratio']
    taper_h = ht['taper_ratio']
    sweep_h = ht['sweep']
    V_v = vt['volume_coefficient']
    A_v = vt['aspect_ratio']
    taper_v = vt['taper_ratio']
    sweep_v = vt['sweep']

    # ---------------------------------------------------------
    # HORIZONTAL STABILIZER SIZING
    # ---------------------------------------------------------
    # 1. Calculate new area based on volume coefficient and calculated moment arm
    S_h = (V_h * S * w.MAC) / ht_moment_arm

    # 2. Update geometry dependent on the new area
    b_h = np.sqrt(A_h * S_h)
    c_r_h = 2 * ht.area / (ht.span * (1 + taper_h))
    c_t_h = c_r_h * taper_h
    
    MAC_h = (2 / 3) * c_r_h * (1 + taper_h + taper_h**2) / (1 + taper_h)
    y_MAC_h = (c_r_h - MAC_h) / (c_r_h * (1 - taper_h)) * (b_h / 2)

    ht['area'] = S_h
    ht['b_h'] = b_h
    ht['c_r_h'] = c_r_h
    ht['c_t_h'] = c_t_h
    ht['MAC_h'] = MAC_h
    ht['y_MAC_h'] = y_MAC_h
    ht['sweep_LE_deg'] = np.rad2deg(np.arctan(np.tan(np.deg2rad(sweep_h)) + 0.5 * c_r_h / b_h * (1 - taper_h)))
    ht['x_le'] = x_h - y_MAC_h * np.tan(np.deg2rad(ht['sweep_LE_deg'])) - 0.4 * MAC_h

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
    vt['sweep_LE_deg'] = np.rad2deg(np.arctan(np.tan(np.deg2rad(sweep_v)) + 0.5 * c_r_v / b_v * (1 - taper_v)))
    vt['x_le'] = x_v - y_MAC_v * np.tan(np.deg2rad(vt['sweep_LE_deg'])) - 0.4 * MAC_v



if __name__ == '__main__':

    # ============================================
    # DATA
    # ============================================

    MTOW = 1870.0  # kg

    concepts = {

        "Taildragger Configuration": {
            "x": [2.90, 3.23, 3.19, 2.88, 2.90],
            "m": [1046.0, 1750.0, 1870.0, 1166.0, 1046.0],
            "labels": ["OEW", "OEW+WP", "OEW+WP+WF", "OEW+WF", ""]
        },

        "Tricycle Configuration": {
            "x": [2.85, 3.20, 3.16, 2.83, 2.85],
            "m": [1046.0, 1750.0, 1870.0, 1166.0, 1046.0],
            "labels": ["OEW", "OEW+WP", "OEW+WP+WF", "OEW+WF", ""]
        }
    }

    # ============================================
    # PLOT
    # ============================================

    plt.figure(figsize=(11,7))

    markers = ['o', 's']

    for i, (name, data) in enumerate(concepts.items()):

        xcg = data["x"]
        mass_fraction = [m / MTOW for m in data["m"]]
        labels = data["labels"]

        plt.plot(
            xcg,
            mass_fraction,
            marker=markers[i],
            linewidth=2.5,
            markersize=8,
            label=name
        )

        # ============================================
        # LABEL POSITIONS
        # ============================================

        if name == "Taildragger Configuration":

            offsets = [
                (25, -10),   # OEW
                (38, -5),    # OEW+WP
                (15, -30),   # OEW+WP+WF
                (-70, -15),  # OEW+WF
                (0, 0)
            ]

        else:  # Tricycle

            offsets = [
                (-60, 15),    # OEW
                (55, 20),     # OEW+WP
                (-120, 20),   # OEW+WP+WF
                (-85, 5),     # OEW+WF
                (0, 0)
            ]

        # ============================================
        # ANNOTATIONS WITH ARROWS
        # ============================================

        for x, y, label, offset in zip(xcg, mass_fraction, labels, offsets):

            if label != "":

                dx, dy = offset

                plt.annotate(
                    label,
                    xy=(x, y),
                    xytext=(dx, dy),
                    textcoords='offset points',
                    fontsize=10,
                    arrowprops=dict(
                        arrowstyle='->',
                        lw=1
                    )
                )

        # ============================================
        # MOST AFT CG
        # ============================================

        most_aft_cg = max(xcg)
        most_fwd_cg = min(xcg)

        print(f"{name}:")
        print(f"  Most aft CG = {most_aft_cg:.2f} m\n")
        print(f"  Most fwd CG = {most_fwd_cg:.2f} m\n")

    # ============================================
    # FIGURE SETTINGS
    # ============================================

    plt.xlabel(r'$x_{cg}$ [m]', fontsize=16)
    plt.ylabel(r'Mass fraction, $M/M_{TO}$ [-]', fontsize=16)

    plt.xlim(2.75, 3.30)
    plt.ylim(0.5, 1.05)

    plt.grid(True, alpha=0.4)

    plt.legend(fontsize=12)

    plt.title('Class I Loading Diagram', fontsize=20)

    plt.tight_layout()

    plt.show()