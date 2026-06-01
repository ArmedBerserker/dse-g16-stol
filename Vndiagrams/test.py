"""
V-n Diagram and Structural Envelope Generation Utilities.

Generates maneuvering envelopes based on CS-23 regulations.
Creates individual plots per weight and altitude condition.
"""

import sys
import os
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

from classes.aircraft_2 import Aircraft, loader
from classes.isa import Atmosphere
from lookups.consts import *
import numpy as np
import matplotlib.pyplot as plt

BASE_DIR = Path(__file__).resolve().parent


# ============================================================
# CHARACTERISTIC SPEEDS
# ============================================================

def calculate_characteristic_speeds(ac: Aircraft, rho: float, weight: float):

    m = ac.weights.m_takeoff * (1 / LBS_TO_KG)
    S = ac.wing.area * M2_TO_F2
    n_max = ac.requirements.general['n_max']

    V_s_la = np.sqrt((2 * weight) / (rho * ac.wing.area * ac.requirements.landing['as_CL_max_la']))
    V_s_to = np.sqrt((2 * weight) / (rho * ac.wing.area * ac.requirements.take_off['as_CL_max_to']))
    V_s_clean = np.sqrt((2 * weight) / (rho * ac.wing.area * ac.requirements.climb['as_CL_max']))

    V_c_min = (33 * np.sqrt(m / S)) * KTS_TO_MS
    V_c = V_c_min

    V_d_min1 = 1.25 * V_c
    V_d_min2 = 1.4 * V_c_min
    V_d = max(V_d_min1, V_d_min2)

    V_a = V_s_clean * np.sqrt(n_max)

    return V_s_clean, V_c, V_d, V_a, V_s_la, V_s_to


# ============================================================
# ENVELOPE
# ============================================================

def generate_vn_envelope(ac: Aircraft, flight: str, condition: str,
                         altitude_m: float, weight: float):

    atmos = Atmosphere(altitude_m)
    rho = atmos.density

    V_s_clean, V_c, V_d, V_a, V_s_la, V_s_to = calculate_characteristic_speeds(ac, rho, weight)

    V_vec = np.linspace(0, V_d, 500)

    n_max = ac.requirements.general['n_max']
    n_min = ac.requirements.general['n_min']
    n_max_flaps = 2.0

    # Positive envelope
    n_pos = np.where(V_vec <= V_a, (V_vec / V_s_clean) ** 2, n_max)

    # Flaps
    n_flap_la = np.minimum((V_vec / V_s_la) ** 2, n_max_flaps)
    n_flap_to = np.minimum((V_vec / V_s_to) ** 2, n_max_flaps)

    # Negative envelope
    n_neg_stall = -(V_vec / V_s_clean) ** 2
    n_neg = np.where(V_vec <= V_s_clean * np.sqrt(abs(n_min)), n_neg_stall, n_min)

    mask = V_vec > V_c
    n_neg[mask] = n_min + (0 - n_min) * (V_vec[mask] - V_c) / (V_d - V_c)

    return {
        "V": V_vec,
        "n_pos": n_pos,
        "n_neg": n_neg,
        "n_flap_la": n_flap_la,
        "n_flap_to": n_flap_to,
        "speeds": {
            "Vsla": V_s_la,
            "Vsto": V_s_to,
            "Vsclean": V_s_clean,
            "Vc": V_c,
            "Vd": V_d,
            "Va": V_a
        }
    }


# ============================================================
# PLOTTING (SINGLE CASE)
# ============================================================

def plot_vn_diagram(ac: Aircraft, flight: str, condition: str,
                    altitude_ft: float, show_plot: bool = False):

    altitude_m = altitude_ft * FT_TO_M

    weight = (
        ac.weights.m_takeoff * g if condition == "MTOW"
        else ac.weights.m_empty * g
    )

    env = generate_vn_envelope(ac, flight, condition, altitude_m, weight)
    V = env["V"]
    speeds = env["speeds"]

    fig, ax = plt.subplots(figsize=(10, 6))

    # Maneuver envelope
    ax.plot(V, env["n_pos"], 'k-', linewidth=2.5, label="Maneuver Envelope")
    ax.plot(V, env["n_neg"], 'k-', linewidth=2.5)

    # Flaps (landing)
    V_f_la = max(1.8 * speeds["Vsla"], 1.4 * speeds["Vsclean"])
    mask_la = V <= V_f_la

    ax.plot(
        np.append(V[mask_la], V_f_la),
        np.append(env["n_flap_la"][mask_la], 0.0),
        'b', linewidth=2,
        label=f"Landing Flaps (Vf={V_f_la:.1f} m/s)"
    )

    # Flaps (takeoff)
    V_f_to = max(1.8 * speeds["Vsto"], 1.4 * speeds["Vsclean"])
    mask_to = V <= V_f_to

    ax.plot(
        np.append(V[mask_to], V_f_to),
        np.append(env["n_flap_to"][mask_to], 0.0),
        'r', linewidth=2,
        label=f"Takeoff Flaps (Vf={V_f_to:.1f} m/s)"
    )

    # Reference lines
    ax.axhline(0, color='black', lw=1)
    ax.axhline(1, color='gray', ls=':', alpha=0.7)

    ax.set_xlabel("Equivalent Airspeed (EAS) [m/s]")
    ax.set_ylabel("Load Factor (n)")
    ax.set_title(f"V-n Diagram: {condition} @ {altitude_ft} ft")
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend()

    plt.tight_layout()

    # Save per case
    output_dir = BASE_DIR / "outputs" / "vn_diagrams"
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = f"VN_{condition}_{int(altitude_ft)}ft_{flight}.png"
    plt.savefig(output_dir / filename, dpi=300)

    if show_plot:
        plt.show()

    plt.close()


# ============================================================
# DRIVER (ALL CASES)
# ============================================================

def generate_all_vn_cases(ac: Aircraft):

    altitudes_ft = [0, 5000, 10000]
    conditions = ["MTOW", "OEW"]
    flight = "cruise"

    for alt in altitudes_ft:
        for cond in conditions:
            plot_vn_diagram(ac, flight, cond, alt)


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    file_path = BASE_DIR.parent / "yamls" / "aircraft.yaml"
    ac = loader.load(file_path, Aircraft)

    generate_all_vn_cases(ac)