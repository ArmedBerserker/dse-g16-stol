"""
Gust and Structural Envelope Generation

Generates one Gust diagram per (altitude × weight condition)
based on CS-23 maneuvering envelope rules.
"""

import sys
import os
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Ensure project root is in path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

from classes.aircraft_2 import Aircraft, loader
from classes.isa import Atmosphere
from lookups.consts import *

BASE_DIR = Path(__file__).resolve().parent


# ============================================================
# CHARACTERISTIC SPEEDS
# ============================================================

def calculate_characteristic_speeds(ac: Aircraft, rho: float, weight: float):
    """Calculates V-n characteristic speeds."""
    m = ac.weights.m_takeoff * (1/LBS_TO_KG)
    S = ac.wing.area * M2_TO_F2


    # Design Speeds
    # Ensures V_c is at least 33 * sqrt(W/S)
    V_c_min = (33 * np.sqrt(m / S)) * KTS_TO_MS
    #V_c_a = ac.requirements.cruise['cr_speed'] * KTS_TO_MS * np.sqrt(rho/1.225)
    #V_c = max(V_c_a, V_c_min)
    #V_c = V_c_min
    V_c = 132 * KTS_TO_MS
    #if V_c_a < V_c_min:
        #print(f"Your cruise speed is too low to adhere to CS23, it needs to at least be: {V_c_min:.2f} m/s")

    V_d_min1 = 1.25 * V_c
    V_d_min2 = 1.4 * V_c_min
    V_d = max(V_d_min1, V_d_min2) #From CS 23.335 reqs


    return  V_c, V_d


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

def gust_load_at_speed(ac, V, rho, weight, Ude):
    """
    Returns upper and lower gust load factor at a single speed.
    """
    mu_g = (2 * (weight / ac.wing.area)) / (
            rho * np.rad2deg(ac.requirements.climb['lift_slope']) * g * ac.requirements.general['mac'])  #lift slope given in 1/deg

    Kg = (0.88 * mu_g) / (5.3 + mu_g)

    dn = (Kg * 1.225 * V * np.rad2deg(ac.requirements.climb['lift_slope']) * Ude) / (2 * (weight / ac.wing.area))

    return 1 + dn, 1 - dn


# ============================================================
# ENVELOPE GENERATION
# ============================================================

def generate_gust_envelope(ac: Aircraft, altitude_m: float, weight: float):
    """
    Computes the load factor limits for the maneuvering envelope.
    """
    # Atmosphere Setup
    atmos = Atmosphere(altitude_m)
    rho = atmos.density

    V_c, V_d =  calculate_characteristic_speeds(ac, rho, weight)

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
    #Ude_Vb = 66 * FT_TO_M    # intermediate (66 ft/s)

    # Compute gust curves
    n_g_vc_up, n_g_vc_low = compute_gust_lines(ac, V_vec, rho, weight, Ude_Vc)
    n_g_vd_up, n_g_vd_low = compute_gust_lines(ac, V_vec, rho, weight, Ude_Vd)
    #n_g_vb_up, n_g_vb_low = compute_gust_lines(ac, V_vec, rho, weight, Ude_Vb)

    return {
        "V": V_vec,

        #Gusts
        "n_g_vc_up": n_g_vc_up,
        "n_g_vc_low": n_g_vc_low,
        "n_g_vd_up": n_g_vd_up,
        "n_g_vd_low": n_g_vd_low,

         "speeds": {"Vc": V_c, "Vd": V_d}
    }


# ============================================================
# PLOTTING FUNCTION (ONE CASE)
# ============================================================

def plot_gust_diagram(ac: Aircraft, altitude_ft: float, condition: str, show_plot: bool = False):

    altitude_m = altitude_ft * FT_TO_M
    atmos = Atmosphere(altitude_m)
    rho = atmos.density

    weight = (
        ac.weights.m_takeoff * g if condition == "MTOW"
        else ac.weights.m_empty * g
    )

    results = generate_gust_envelope(ac, altitude_m, weight)
    V = results["V"]
    speeds = results["speeds"]

    Vc = speeds["Vc"]
    Vd = speeds["Vd"]

    Ude_Vc = 50 * FT_TO_M
    Ude_Vd = 25 * FT_TO_M

    nCp, nCm = gust_load_at_speed(
        ac, Vc, rho, weight, Ude_Vc)

    nDp, nDm = gust_load_at_speed(
        ac, Vd, rho, weight, Ude_Vd)

    fig, ax = plt.subplots(figsize=(10, 6))

    #Plot Envelopes
    n_max = ac.requirements.general['n_max']
    n_min = ac.requirements.general['n_min']



    # A-C'-D'-E'-F'-A envelope

    poly_x = [0, Vc, Vd, Vd, Vc, 0]
    poly_y = [1, nCp, nDp, nDm, nCm, 1]
    n_max_env = max(poly_y)

    ax.plot(poly_x, poly_y, color='black', linewidth=2, label='Gust Envelope', zorder=10)

    # Plot gust envelopes
    ax.plot(V, results["n_g_vc_up"], 'c--', label=f'Gust line ({Ude_Vc} m/s ; 50 ft/s)')
    ax.plot(V, results["n_g_vc_low"], 'c--')

    ax.plot(V, results["n_g_vd_up"], 'm--', label=f'Gust line ({Ude_Vd} m/s ; 25 ft/s)')
    ax.plot(V, results["n_g_vd_low"], 'm--')

    # Positive gust rays
    #ax.plot([0, Vc],[1, nCp],'k--',alpha=0.7)
    #ax.plot([0, Vd],[1, nDp],'k--',alpha=0.7)

    # Negative gust rays
    #ax.plot([0, Vc],[1, nCm],'k--',alpha=0.7)
    #ax.plot([0, Vd],[1, nDm],'k--',alpha=0.7)

    ax.axhline(y=n_max_env, color='black', linestyle=':', linewidth=1.5)

    ax.text(0.98 * Vd,  # x-position
            n_max_env + 0.05,  # y-position
            fr'$n_{{max}} = {n_max_env:.2f}$', color='black', ha='right', va='bottom')


    # Vertical line at Vd
    ax.plot([Vd, Vd],[nDm, nDp],'k',linewidth=2)

    # Reference Lines
    ax.axhline(0, color='black', lw=1)
    ax.axhline(1, color='gray', ls=':', alpha=0.7)

    # Annotate Speeds
    #ax.axvline(speeds["Vsclean"], color='r', ls='--', alpha=0.5, label=f'Vs ({speeds["Vsclean"]:.1f} m/s)')
    #ax.axvline(speeds["Vc"], color='g', ls='--', alpha=0.5, label=f'Vc ({speeds["Vc"]:.1f} m/s)')
    #ax.axvline(speeds["Vd"], color='m', ls='--', alpha=0.5, label=f'Vd ({speeds["Vd"]:.1f} m/s)')

    #ax.axvline(speeds["Vb"], color='orange', linestyle='--', alpha=0.5, label=f'$V_B$ ({speeds["Vb"]:.1f} m/s)')

    # Vertical speed markers only inside envelope
    #ax.plot([speeds["Vc"], speeds["Vc"]],[n_max, n_min],'g--', alpha=0.5, label=f'$V_C$ ({speeds["Vc"]:.1f} m/s)')
    #ax.plot([speeds["Vd"], speeds["Vd"]], [n_max, 0], color='black', linestyle='--', alpha=0.5, label=f'$V_D$ ({speeds["Vd"]:.1f} m/s)')
    ax.plot([Vc, Vc],[nCm, nCp],'g--',alpha=0.6, label=f'$V_C$ ({speeds["Vc"]:.1f} m/s)')
    ax.plot([Vd, Vd],[nDm, nDp],'k--',alpha=0.6, label=f'$V_D$ ({speeds["Vd"]:.1f} m/s)')

    # Labels for gust velocities

    #ax.annotate(r"$50$ ft/s",xy=(Vc, nCp),xytext=(10, 10),textcoords="offset points")
    #ax.annotate(r"$25$ ft/s",xy=(Vd, nDp),xytext=(10, 10),textcoords="offset points")
    # Position labels at 90% of Vd
    x_label = 0.85 * Vd

    # Interpolate y-values on the gust lines
    y_vc_up = np.interp(x_label, V, results["n_g_vc_up"])
    y_vc_low = np.interp(x_label, V, results["n_g_vc_low"])

    y_vd_up = np.interp(x_label, V, results["n_g_vd_up"])
    y_vd_low = np.interp(x_label, V, results["n_g_vd_low"])

    angle_vc = np.degrees(np.arctan2(results["n_g_vc_up"][-1] - 1,V[-1]))
    angle_vd = np.degrees(np.arctan2(results["n_g_vd_up"][-1] - 1,V[-1]))

    ax.text(x_label - 1.5,y_vc_up + 0.25,fr"${Ude_Vc}$ m/s",rotation=angle_vc,fontsize=10,ha='left')
    ax.text(x_label - 1,y_vd_up + 0.25,fr"${Ude_Vd}$ m/s",rotation=angle_vd,fontsize=10,ha='left')



    # Labels below axis
    y_text = -0.2
    ax.text(speeds["Vc"] + 1.5, y_text, r"$V_C$", ha='center', va='top', fontsize=11, color='g')
    ax.text(speeds["Vd"] - 1.5, y_text,r"$V_D$", ha='center', va='top', fontsize=11, color='black')
    #ax.text(speeds["Vb"],-0.2,r"$V_B$", ha='center', va='top', fontsize=11, color='orange')


    ax.set_xlabel('Equivalent Airspeed (EAS) [m/s]')
    ax.set_ylabel('Load Factor (n)')
    #ax.set_title(f'V-n Diagram: MTOW={ac.weights.m_takeoff:.0f} kg')
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.5)


    plt.tight_layout()

    # Save output
    output_dir = BASE_DIR / "outputs" / "vn_diagrams"
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = f"Gust_{condition}_{int(altitude_ft)}ft.png"
    plt.savefig(output_dir / filename, dpi=300)

    if show_plot:
        plt.show()

    plt.close()


# ============================================================
# DRIVER (ALL CASES)
# ============================================================

def generate_all_gust_cases(ac: Aircraft):

    altitudes_ft = [0, ac.requirements.take_off['to_altitude'], ac.requirements.cruise['cr_altitude']]
    conditions = ["MTOW", "OEW"]

    for alt in altitudes_ft:
        for cond in conditions:
            plot_gust_diagram(ac, alt, cond)


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    file_path = BASE_DIR.parent / "yamls" / "aircraft.yaml"
    ac = loader.load(file_path, Aircraft)

    generate_all_gust_cases(ac)