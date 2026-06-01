"""
Gust and Structural Envelope Generation Utilities.

Generates gust envelopes based on CS-23 regulations.
Calculates limit load factors (n) across the airspeed range, accounting for
structural design speeds (Vc, Vd).
"""

import sys
import os

# Fix path FIRST, before any local imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from classes.aircraft_2 import Aircraft, loader
from classes.isa import Atmosphere
from lookups.consts import *
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

BASE_DIR = Path(__file__).resolve().parent


def calculate_characteristic_speeds(ac: Aircraft, rho: float, weight: float):
    """Calculates V-n characteristic speeds."""
    m = ac.weights.m_takeoff * (1/LBS_TO_KG)
    S = ac.wing.area * M2_TO_F2

    #Stall Speeds
    V_s_clean = np.sqrt((2 * weight) / (rho * ac.wing.area * ac.requirements.climb['as_CL_max'])) * np.sqrt(rho / 1.225)

    # Design Speeds
    # Ensures V_c is at least 33 * sqrt(W/S)
    V_c_min = (33 * np.sqrt(m / S)) * KTS_TO_MS
    #V_c_a = ac.requirements.cruise['cr_speed'] * KTS_TO_MS * np.sqrt(rho/1.225)
    #V_c = max(V_c_a, V_c_min)
    V_c = V_c_min
    #if V_c_a < V_c_min:
        #print(f"Your cruise speed is too low to adhere to CS23, it needs to at least be: {V_c_min:.2f} m/s")

    V_d_min1 = 1.25 * V_c
    V_d_min2 = 1.4 * V_c_min
    V_d = max(V_d_min1, V_d_min2) #From CS 23.335 reqs


    return  V_c, V_d, V_s_clean


def compute_gust_lines(ac: Aircraft, V, rho, weight, Ude):
    """
    Returns gust envelopes for VB, VC, VD.
    """

    # Gust alleviation factor (CS-23 approximation)
    mu_g = (2 * (weight / ac.wing.area)) / (
            rho * np.rad2deg(ac.requirements.climb['lift_slope']) * g * ac.requirements.general['mac'])  #lift slope given in 1/deg

    Kg = (0.88 * mu_g) / (5.3 + mu_g)

    dn = (Kg * 1.225 * V * np.rad2deg(ac.requirements.climb['lift_slope']) * Ude) / (2 * (weight / ac.wing.area))

    n_upper = 1 + dn
    n_lower = 1 - dn

    return n_upper, n_lower

def generate_gust_envelope(ac: Aircraft, flight: str = 'cruise', condition: str = 'MTOW'):
    """
    Computes the load factor limits for the maneuvering envelope.
    """
    # Atmosphere Setup
    altitude = ac.requirements.cruise['cr_altitude'] * FT_TO_M if flight == 'cruise' else 0
    atmos = Atmosphere(altitude)
    rho = atmos.density
    temp = atmos.temp

    weight = (
        ac.weights.m_takeoff * g
        if condition == 'MTOW'
        else ac.weights.m_empty * g
    )

    V_c, V_d, V_s_clean =  calculate_characteristic_speeds(ac, rho, weight)

    # Airspeed Vector
    V_vec = np.linspace(0, V_d, 500)


    # Maneuvering limit loads
    n_max = ac.requirements.general['n_max']


    # -----------------------------
    # GUST ENVELOPES (CS-23.341)
    # -----------------------------

    # Design gust velocities (simplified CS-23 values)
    Ude_Vc = 50 * FT_TO_M    # 50 ft/s at VC
    Ude_Vd = 25 * FT_TO_M    # reduced at VD (25 ft/s)
    Ude_Vb = 66 * FT_TO_M    # intermediate (66 ft/s)

    # Compute gust curves
    n_g_vc_up, n_g_vc_low = compute_gust_lines(ac, V_vec, rho, weight, Ude_Vc)
    n_g_vd_up, n_g_vd_low = compute_gust_lines(ac, V_vec, rho, weight, Ude_Vd)
    n_g_vb_up, n_g_vb_low = compute_gust_lines(ac, V_vec, rho, weight, Ude_Vb)

    # ------------------------------------------------
    # MANEUVER ENVELOPE
    # ------------------------------------------------

    n_stall_pos = (V_vec / V_s_clean) ** 2

    # Positive maneuvering speed
    V_a = V_s_clean * np.sqrt(n_max)

    # ------------------------------------------------
    # VB = intersection of positive stall boundary
    #      and 66 ft/s gust line
    # ------------------------------------------------

    diff = n_stall_pos - n_g_vb_up

    crossings = np.where(np.diff(np.sign(diff)))[0]

    if len(crossings) > 0:
        idx = crossings[0]

        # Linear interpolation for better accuracy
        V1 = V_vec[idx]
        V2 = V_vec[idx + 1]

        d1 = diff[idx]
        d2 = diff[idx + 1]

        V_b = V1 - d1 * (V2 - V1) / (d2 - d1)

    else:
        # Fallback if no crossing is found
        V_b = min(V_c, V_a)

    # Regulatory bound
    V_b = np.clip(V_b, V_s_clean, V_c)

    vb_mask = V_vec <= V_b

    V_vb = V_vec[vb_mask]

    n_g_vb_up_plot = n_g_vb_up[vb_mask]
    n_g_vb_low_plot = n_g_vb_low[vb_mask]


    return {
        "V_vb": V_vb,

        "n_stall_pos": n_stall_pos,

        "n_g_vb_up_plot": n_g_vb_up_plot,
        "n_g_vb_low_plot": n_g_vb_low_plot,


        "V": V_vec,

        #Gusts
        "n_g_vc_up": n_g_vc_up,
        "n_g_vc_low": n_g_vc_low,
        "n_g_vd_up": n_g_vd_up,
        "n_g_vd_low": n_g_vd_low,
        "n_g_vb_up": n_g_vb_up,
        "n_g_vb_low": n_g_vb_low,

         "speeds": {"Vs": V_s_clean, "Va": V_a, "Vb": V_b, "Vc": V_c, "Vd": V_d}
    }


def plot_gust_diagram(ac: Aircraft, output_filepath: str = 'outputs/Gust_Diagram.png', show_plot: bool = False):
    """
    Generates and saves the V-n Diagram plot.
    """

    results = generate_gust_envelope(ac)
    V = results["V"]
    speeds = results["speeds"]

    fig, ax = plt.subplots(figsize=(10, 6))

    #Plot Envelopes
    n_max = ac.requirements.general['n_max']
    n_min = ac.requirements.general['n_min']

    # -----------------------------
    # MANEUVER ENVELOPE
    # -----------------------------

    Va_idx = np.argmin(np.abs(V - speeds["Va"]))

    ax.plot(
        V[:Va_idx + 1],
        results["n_stall_pos"][:Va_idx + 1],
        'b',
        linewidth=2,
        label='Positive stall boundary'
    )

    ax.plot(
        [speeds["Va"], speeds["Vd"]],
        [n_max, n_max],
        'b',
        linewidth=2
    )

    # Plot gust envelopes
    ax.plot(V, results["n_g_vc_up"], 'c--', label='Gust VC')
    ax.plot(V, results["n_g_vc_low"], 'c--')

    ax.plot(V, results["n_g_vd_up"], 'm--', label='Gust VD')
    ax.plot(V, results["n_g_vd_low"], 'm--')

    ax.plot(results["V_vb"], results["n_g_vb_up_plot"], 'g--', linewidth=2, label='VB Gust (66 ft/s)')
    ax.plot(results["V_vb"], results["n_g_vb_low_plot"], 'g--', linewidth=2)


    # Vertical line at Vd
    ax.plot([speeds["Vd"], speeds["Vd"]], [0, n_max],'k-', linewidth=2)

    # Reference Lines
    ax.axhline(0, color='black', lw=1)
    ax.axhline(1, color='gray', ls=':', alpha=0.7)

    # Annotate Speeds
    #ax.axvline(speeds["Vsclean"], color='r', ls='--', alpha=0.5, label=f'Vs ({speeds["Vsclean"]:.1f} m/s)')
    #ax.axvline(speeds["Vc"], color='g', ls='--', alpha=0.5, label=f'Vc ({speeds["Vc"]:.1f} m/s)')
    #ax.axvline(speeds["Vd"], color='m', ls='--', alpha=0.5, label=f'Vd ({speeds["Vd"]:.1f} m/s)')

    ax.axvline(speeds["Vb"], color='orange', linestyle='--', alpha=0.5, label=f'$V_B$ ({speeds["Vb"]:.1f} m/s)')

    # Vertical speed markers only inside envelope
    ax.plot([speeds["Vc"], speeds["Vc"]],[n_max, n_min],'g--', alpha=0.5, label=f'$V_C$ ({speeds["Vc"]:.1f} m/s)')
    ax.plot([speeds["Vd"], speeds["Vd"]], [n_max, 0], color='black', linestyle='--', alpha=0.5, label=f'$V_D$ ({speeds["Vd"]:.1f} m/s)')


    # Labels below axis
    y_text = -0.2
    ax.text(speeds["Vc"] + 1.5, y_text, r"$V_C$", ha='center', va='top', fontsize=11, color='g')
    ax.text(speeds["Vd"], y_text,r"$V_D$", ha='center', va='top', fontsize=11, color='black')
    ax.text(speeds["Vb"],-0.2,r"$V_B$", ha='center', va='top', fontsize=11, color='orange')


    ax.set_xlabel('Equivalent Airspeed (EAS) [m/s]')
    ax.set_ylabel('Load Factor (n)')
    #ax.set_title(f'V-n Diagram: MTOW={ac.weights.m_takeoff:.0f} kg')
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(output_filepath, dpi=300)
    if show_plot:
        plt.show()


if __name__ == '__main__':
    # Load aircraft parameters from the centralized YAML configuration
    file_path = BASE_DIR.parent / "yamls" / "aircraft.yaml"
    target_class = Aircraft
    ac = loader.load(file_path, Aircraft)

    # Generate the diagram
    plot_gust_diagram(ac)