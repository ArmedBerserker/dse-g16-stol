"""
Gust and Structural Envelope Generation Utilities.

Generates gust envelopes based on CS-23 regulations.
Calculates limit load factors (n) across the airspeed range, accounting for
structural design speeds (Vc, Vd).
"""

import sys
import os

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
        else ac.weights.m_empty * g)

    V_c, V_d, V_s_clean =  calculate_characteristic_speeds(ac, rho, weight)

    V_vec = np.linspace(0, V_d, 500)

    n_max = ac.requirements.general['n_max']
    n_min = ac.requirements.general['n_min']

    # -----------------------------
    # GUST ENVELOPES (CS-23)
    # -----------------------------

    # Design gust velocities (simplified CS-23 values)
    Ude_Vc = 50 * FT_TO_M
    Ude_Vd = 25 * FT_TO_M
    Ude_Vb = 66 * FT_TO_M

    # Compute gust curves
    n_g_vc_up, n_g_vc_low = compute_gust_lines(ac, V_vec, rho, weight, Ude_Vc)
    n_g_vd_up, n_g_vd_low = compute_gust_lines(ac, V_vec, rho, weight, Ude_Vd)
    n_g_vb_up, n_g_vb_low = compute_gust_lines(ac, V_vec, rho, weight, Ude_Vb)

    # ------------------------------------------------
    # MANEUVER ENVELOPE
    # ------------------------------------------------

    # Positive maneuvering speed
    V_a = V_s_clean * np.sqrt(n_max)

    # ---- VB ----
    def n_stall(V):
        return (V / V_s_clean) ** 2

    W_S = weight / ac.wing.area
    a = ac.requirements.climb['lift_slope']

    Ude_Vb = 66 * FT_TO_M

    mu_g = (2 * W_S) / (rho * np.rad2deg(a) * g * ac.requirements.general['mac'])
    Kg = (0.88 * mu_g) / (5.3 + mu_g)

    C_vb = (Kg * rho * np.rad2deg(a) * Ude_Vb) / (2 * W_S)

    def n_gust_vb(V):
        return 1 + C_vb * V

    A = 1 / V_s_clean ** 2
    B = -C_vb
    Cq = -1

    roots = np.roots([A, B, Cq])

    V_b_candidates = [r.real for r in roots if np.isreal(r) and r > 0]

    V_b = min(V_b_candidates) if V_b_candidates else V_a

    vb_mask = V_vec <= V_b

    n_stall_pos = n_stall(V_vec)

    # ---- CERTIFICATION POINTS ----
    n_B_plus  = np.interp(V_b, V_vec, n_g_vb_up)
    n_B_minus = np.interp(V_b, V_vec, n_g_vb_low)

    n_C_plus  = np.interp(V_c, V_vec, n_g_vc_up)
    n_C_minus = np.interp(V_c, V_vec, n_g_vc_low)

    n_D_plus  = np.interp(V_d, V_vec, n_g_vd_up)
    n_D_minus = np.interp(V_d, V_vec, n_g_vd_low)

    gust_points = {
        "B+": (V_b, n_B_plus),
        "B-": (V_b, n_B_minus),
        "C+": (V_c, n_C_plus),
        "C-": (V_c, n_C_minus),
        "D+": (V_d, n_D_plus),
        "D-": (V_d, n_D_minus),
    }

    return {
        "V": V_vec,
        "n_stall_pos": n_stall_pos,

        "n_g_vc_up": n_g_vc_up,
        "n_g_vc_low": n_g_vc_low,
        "n_g_vd_up": n_g_vd_up,
        "n_g_vd_low": n_g_vd_low,
        "n_g_vb_up": n_g_vb_up,
        "n_g_vb_low": n_g_vb_low,

        "V_vb": V_vec[vb_mask],
        "n_g_vb_up_plot": n_g_vb_up[vb_mask],
        "n_g_vb_low_plot": n_g_vb_low[vb_mask],

        "gust_points": gust_points,

        "speeds": {
            "Vs": V_s_clean,
            "Va": V_a,
            "Vb": V_b,
            "Vc": V_c,
            "Vd": V_d
        }
    }


def plot_gust_diagram(ac, output_filepath='outputs/Gust_Diagram.png', show_plot=False):

    results = generate_gust_envelope(ac)

    V = results["V"]
    speeds = results["speeds"]
    gp = results["gust_points"]

    Vs = speeds["Vs"]
    Va = speeds["Va"]
    Vb = speeds["Vb"]
    Vc = speeds["Vc"]
    Vd = speeds["Vd"]

    n_max = ac.requirements.general['n_max']
    n_min = ac.requirements.general['n_min']

    fig, ax = plt.subplots(figsize=(10, 6))

    # -------------------------------------------------
    # Reference gust lines
    # -------------------------------------------------

    ax.plot(V, results["n_g_vc_up"], 'c--')
    ax.plot(V, results["n_g_vc_low"], 'c--')

    ax.plot(V, results["n_g_vd_up"], 'm--')
    ax.plot(V, results["n_g_vd_low"], 'm--')

    ax.plot(V, results["n_g_vb_up"], 'g--')
    ax.plot(V, results["n_g_vb_low"], 'g--')

    # -------------------------------------------------
    # POSITIVE STALL CURVE
    # -------------------------------------------------

    V_stall_pos = np.linspace(Vs, Va, 300)
    n_stall_pos = (V_stall_pos / Vs) ** 2

    # -------------------------------------------------
    # NEGATIVE STALL CURVE
    # -------------------------------------------------

    V_neg_limit = Vs * np.sqrt(abs(n_min))

    V_stall_neg = np.linspace(Vs, V_neg_limit, 300)
    n_stall_neg = -(V_stall_neg / Vs) ** 2

    # -------------------------------------------------
    # GUST POINTS
    # -------------------------------------------------

    Bp = gp["B+"]
    Bm = gp["B-"]

    Cp = gp["C+"]
    Cm = gp["C-"]

    Dp = gp["D+"]
    Dm = gp["D-"]

    # -------------------------------------------------
    # UPPER MANOEUVRE + GUST ENVELOPE
    # -------------------------------------------------

    ax.plot(
        V_stall_pos,
        n_stall_pos,
        'k',
        lw=2,
        label='Upper Envelope'
    )

    ax.plot(
        [Va, Bp[0]],
        [n_max, Bp[1]],
        'k',
        lw=2
    )

    ax.plot(
        [Bp[0], Cp[0], Dp[0]],
        [Bp[1], Cp[1], Dp[1]],
        'k',
        lw=2
    )

    # -------------------------------------------------
    # LOWER MANOEUVRE + GUST ENVELOPE
    # -------------------------------------------------

    ax.plot(
        V_stall_neg,
        n_stall_neg,
        'k',
        lw=2,
        label='Lower Envelope'
    )

    ax.plot(
        [V_neg_limit, Bm[0]],
        [n_min, Bm[1]],
        'k',
        lw=2
    )

    ax.plot(
        [Bm[0], Cm[0], Dm[0]],
        [Bm[1], Cm[1], Dm[1]],
        'k',
        lw=2
    )

    # -------------------------------------------------
    # CLOSE ENVELOPE AT VD
    # -------------------------------------------------

    ax.plot(
        [Vd, Vd],
        [Dm[1], Dp[1]],
        'k',
        lw=2
    )

    # -------------------------------------------------
    # FILL ENVELOPE
    # -------------------------------------------------

    poly_x = np.concatenate([
        V_stall_pos,
        [Bp[0], Cp[0], Dp[0]],
        [Vd],
        [Dm[0], Cm[0], Bm[0]],
        [V_neg_limit],
        V_stall_neg[::-1]
    ])

    poly_y = np.concatenate([
        n_stall_pos,
        [Bp[1], Cp[1], Dp[1]],
        [Dm[1]],
        [Dm[1], Cm[1], Bm[1]],
        [n_min],
        n_stall_neg[::-1]
    ])

    ax.fill(
        poly_x,
        poly_y,
        color='lightgrey',
        alpha=0.3
    )

    # -------------------------------------------------
    # CERTIFICATION POINTS
    # -------------------------------------------------

    ax.scatter(
        [
            Bp[0], Cp[0], Dp[0],
            Bm[0], Cm[0], Dm[0]
        ],
        [
            Bp[1], Cp[1], Dp[1],
            Bm[1], Cm[1], Dm[1]
        ],
        color='black',
        zorder=5
    )

    # -------------------------------------------------
    # SPEED MARKERS
    # -------------------------------------------------

    for speed, label in [
        (Va, "VA"),
        (Vb, "VB"),
        (Vc, "VC"),
        (Vd, "VD")
    ]:
        ax.axvline(
            speed,
            color='black',
            ls='--',
            alpha=0.5
        )

        ax.text(
            speed,
            ax.get_ylim()[1] * 0.95,
            label,
            ha='center'
        )

    ax.set_xlabel("EAS [m/s]")
    ax.set_ylabel("Load factor n")

    ax.grid(True, ls='--', alpha=0.4)

    ax.legend()

    plt.tight_layout()
    plt.savefig(output_filepath, dpi=300)

    if show_plot:
        plt.show()


if __name__ == '__main__':
    file_path = BASE_DIR.parent / "yamls" / "aircraft.yaml"
    ac = loader.load(file_path, Aircraft)

    plot_gust_diagram(ac)