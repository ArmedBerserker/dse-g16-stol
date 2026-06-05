"""
V-n Diagram and Structural Envelope Generation

Generates one V-n diagram per (altitude × weight condition)
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
    n_max = ac.requirements.general['n_max']

    # Stall Speeds
    #V_s_la = ac.requirements.general['stall_speed'] * KTS_TO_MS * np.sqrt(rho/1.225)
    V_s_la = np.sqrt((2 * weight)/(rho * ac.wing.area * ac.requirements.landing['as_CL_max_la'])) * np.sqrt(rho/1.225)
    V_s_to = np.sqrt((2 * weight)/(rho * ac.wing.area * ac.requirements.take_off['as_CL_max_to'])) * np.sqrt(rho/1.225)
    V_s_clean = np.sqrt((2 * weight)/(rho * ac.wing.area * ac.requirements.climb['as_CL_max'])) * np.sqrt(rho/1.225)

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

    V_ne = 0.9 * V_d

    V_a = V_s_clean * np.sqrt(n_max)

    return V_s_clean, V_c, V_d, V_a, V_s_la, V_s_to, V_ne


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

def generate_vn_envelope(ac: Aircraft, altitude_m: float, weight: float):

    atmos = Atmosphere(altitude_m)
    rho = atmos.density

    V_s_clean, V_c, V_d, V_a, V_s_la, V_s_to, V_ne = calculate_characteristic_speeds(ac, rho, weight)

    # Airspeed Vector
    V_vec = np.linspace(0, V_d, 500)

    # Maneuvering limit loads
    n_max = ac.requirements.general['n_max']
    n_min = ac.requirements.general['n_min']
    n_max_flaps = 2.0     #CS23 req

    n_pos = np.where(V_vec <= V_a, (V_vec / V_s_clean) ** 2, n_max)

    # Flap limit (n_flap)
    n_flap_la = np.minimum((V_vec / V_s_la) ** 2, n_max_flaps)

    n_flap_to = np.minimum((V_vec / V_s_to) ** 2, n_max_flaps)

    # Negative maneuvering limit
    n_neg_stall = -(V_vec / V_s_clean) ** 2
    n_neg = np.where(V_vec <= V_s_clean * np.sqrt(abs(n_min)), n_neg_stall, n_min)

    # Assumes linear reduction of negative maneuver loads
    # between Vc and Vd for conceptual envelope closure
    # Closure to V_D (Ramp from V_C to V_D)
    #mask_dive = V_vec > V_c
    #n_neg[mask_dive] = (n_min / (V_d - V_c)) * (V_vec[mask_dive] - V_d)

    mask = V_vec > V_c
    n_neg[mask] = n_min * ( 1 - (V_vec[mask] - V_c) / (V_d - V_c))

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
        "n_pos": n_pos,
        "n_neg": n_neg,
        # Gusts
        "n_g_vc_up": n_g_vc_up,
        "n_g_vc_low": n_g_vc_low,
        "n_g_vd_up": n_g_vd_up,
        "n_g_vd_low": n_g_vd_low,
        "n_flap_la": n_flap_la,
        "n_flap_to": n_flap_to,
        "speeds": {"Vsla": V_s_la, "Vsto": V_s_to, "Vsclean": V_s_clean, "Vc": V_c, "Vd": V_d, "Va": V_a, "Vne": V_ne}
    }


# ============================================================
# PLOTTING FUNCTION (ONE CASE)
# ============================================================

def plot_vn_diagram(ac: Aircraft, altitude_ft: float, condition: str, show_plot: bool = False):

    altitude_m = altitude_ft * FT_TO_M
    atmos = Atmosphere(altitude_m)
    rho = atmos.density

    weight = (
        ac.weights.m_takeoff * g if condition == "MTOW"
        else ac.weights.m_empty * g
    )

    results = generate_vn_envelope(ac, altitude_m, weight)

    V = results["V"]
    speeds = results["speeds"]

    n_max = ac.requirements.general['n_max']
    n_min = ac.requirements.general['n_min']
    n_max_flaps = 2.0

    Vc = speeds["Vc"]
    Vd = speeds["Vd"]

    Ude_Vc = 50 * FT_TO_M
    Ude_Vd = 25 * FT_TO_M

    nCp, nCm = gust_load_at_speed(
        ac, Vc, rho, weight, Ude_Vc)

    nDp, nDm = gust_load_at_speed(
        ac, Vd, rho, weight, Ude_Vd)


    fig, ax = plt.subplots(figsize=(14, 6))


    #Gusts
    # A-C'-D'-E'-F'-A envelope

    poly_x = [0, Vc, Vd, Vd, Vc, 0]
    poly_y = [1, nCp, nDp, nDm, nCm, 1]
    n_max_env = max(poly_y)



    # Plot gust envelopes
    #ax.plot(V, results["n_g_vc_up"], 'k--', label=f'Gust line ({Ude_Vc} m/s ; 50 ft/s)')
    #ax.plot(V, results["n_g_vc_low"], 'k--')

    ax.plot(V, results["n_g_vd_up"], 'gray', linestyle='--')
    ax.plot(V, results["n_g_vd_low"], 'gray', linestyle='--')

    # Positive gust rays
    #ax.plot([0, Vc],[1, nCp],'k--',alpha=0.7)
    #ax.plot([0, Vd],[1, nDp],'k--',alpha=0.7)

    # Negative gust rays
    #ax.plot([0, Vc],[1, nCm],'k--',alpha=0.7)
    #ax.plot([0, Vd],[1, nDm],'k--',alpha=0.7)

    #Plot Envelopes
    ax.plot(V, results["n_pos"], 'k-', linewidth=2.5, label='Maneuvering Envelope', zorder=100)
    ax.plot(V, results["n_neg"], 'k-', linewidth=2.5, zorder=10)

    ax.plot(poly_x, poly_y, 'gray', linestyle='--' , linewidth=2, label='Gust Envelope', zorder=10)

    # Plot Flaps (Limited to Vf)
    # Landing flaps
    V_f_la = max(1.8 * speeds["Vsla"], 1.4 * speeds["Vsclean"])
    mask_la = V <= V_f_la
    V_la_plot = np.append(V[mask_la], V_f_la)
    n_la_plot = np.append(results["n_flap_la"][mask_la], 0.0)

    # ax.plot(V_la_plot, n_la_plot, 'b', label='Flap Envelope Landing')
    ax.plot(V_la_plot, n_la_plot, 'b', linewidth=2, label=fr'$Flap_{{LA}}$ Envelope ($V_{{F,LA}}$ = {V_f_la:.1f} m/s)')

    # Takeoff flaps
    V_f_to = max(1.8 * speeds["Vsto"], 1.4 * speeds["Vsclean"])
    mask_to = V <= V_f_to
    V_to_plot = np.append(V[mask_to], V_f_to)
    n_to_plot = np.append(results["n_flap_to"][mask_to], 0.0)

    # ax.plot(V_to_plot, n_to_plot, 'r', label='Flap Envelope Take-Off')
    ax.plot(V_to_plot, n_to_plot, 'r', linewidth=2, label=fr'$Flap_{{TO}}$ Envelope ($V_{{F,TO}}$ = {V_f_to:.1f} m/s)')

    # Vertical line at Vd
    ax.plot([speeds["Vd"], speeds["Vd"]], [0, n_max],'k-', linewidth=2, zorder=20)

    # Reference Lines
    ax.axhline(0, color='black', lw=1)
    ax.axhline(1, color='gray', ls=':', alpha=0.7)

    # Annotate Speeds
    #ax.axvline(speeds["Vsclean"], color='r', ls='--', alpha=0.5, label=f'Vs ({speeds["Vsclean"]:.1f} m/s)')
    #ax.axvline(speeds["Vc"], color='g', ls='--', alpha=0.5, label=f'Vc ({speeds["Vc"]:.1f} m/s)')
    #ax.axvline(speeds["Vd"], color='m', ls='--', alpha=0.5, label=f'Vd ({speeds["Vd"]:.1f} m/s)')

    # Vertical speed markers only inside envelope
    ax.plot([speeds["Vsclean"], speeds["Vsclean"]],[0, 1], color='orange', linestyle='--', alpha=0.5, label=f'$V_S$ ({speeds["Vsclean"]:.1f} m/s)')
    ax.plot([speeds["Va"], speeds["Va"]], [n_max, 0], 'm--', alpha=0.5, label=f'$V_A$ ({speeds["Va"]:.1f} m/s)')
    ax.plot([speeds["Vc"], speeds["Vc"]],[n_max, n_min],'g--', alpha=0.5, label=f'$V_C$ ({speeds["Vc"]:.1f} m/s)')
    ax.plot([speeds["Vd"], speeds["Vd"]], [n_max, 0], color='black', linestyle='--', alpha=0.5, label=f'$V_D$ ({speeds["Vd"]:.1f} m/s)')
    #ax.plot([speeds["Vne"], speeds["Vne"]], [n_max, 0], color='black', linestyle='--', alpha=0.5, label=f'$V_{{NE}}$ ({speeds["Vne"]:.1f} m/s)')


    # Labels below axis
    y_text = 0.25
    ax.text(speeds["Vsclean"] + 1.5, y_text,r"$V_S$", ha='center', va='top', fontsize=11, color='orange')
    ax.text(speeds["Va"] + 1.5, y_text, r"$V_A$", ha='center', va='top', fontsize=11, color='m')
    ax.text(speeds["Vc"] + 1.5, y_text, r"$V_C$", ha='center', va='top', fontsize=11, color='g')
    ax.text(speeds["Vd"] + 1.5, y_text,r"$V_D$", ha='center', va='top', fontsize=11, color='black')
    #ax.text(speeds["Vne"] + 1.8, y_text, r"$V_{NE}$", ha='center', va='top', fontsize=11, color='black')

    # Flap speed labels
    y_text1 = -0.1
    ax.text(V_f_la, y_text1, r"$V_{F,LA}$", ha='center', va='top', fontsize=11, color='b')
    y_text2 = -0.3
    ax.text(V_f_to, y_text2, r"$V_{F,TO}$", ha='center', va='top', fontsize=11, color='r')

    # Labels
    ax.set_xlabel("Equivalent Airspeed (EAS) [m/s]")
    ax.set_ylabel("Load Factor (n)")
    #ax.set_title(f"V-n Diagram: {condition} @ {altitude_ft} ft")
    ax.legend(loc='upper left', framealpha=0.9, fontsize=10)
    ax.grid(True, linestyle="--", alpha=0.5)


    plt.tight_layout()

    # Save output
    output_dir = BASE_DIR / "outputs" / "vn_diagrams"
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = f"VN_G_{condition}_{int(altitude_ft)}ft.png"
    plt.savefig(output_dir / filename, dpi=300)

    if show_plot:
        plt.show()

    plt.close()


# ============================================================
# DRIVER (ALL CASES)
# ============================================================

def generate_all_vn_cases(ac: Aircraft):

    altitudes_ft = [0, ac.requirements.take_off['to_altitude'], ac.requirements.cruise['cr_altitude']]
    conditions = ["MTOW", "OEW"]

    for alt in altitudes_ft:
        for cond in conditions:
            plot_vn_diagram(ac, alt, cond)


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    file_path = BASE_DIR.parent / "yamls" / "aircraft.yaml"
    ac = loader.load(file_path, Aircraft)

    generate_all_vn_cases(ac)