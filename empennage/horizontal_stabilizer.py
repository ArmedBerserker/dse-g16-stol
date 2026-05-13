import numpy as np

# ======================================================
# INPUT PARAMETERS
# ======================================================

# Main
X_aftcg = 1
l_fuselage = 8.8
MAC = 1.5
S_w = 23
c = 1.5
AR = 11

# Horizontal stabilizer
X_h = 3
V_h = 0.786

AR_h = 4
Sweep_h_LE = 25
Taper_h = 0.4

t_c_h = 0.12


# ======================================================
# HORIZONTAL STABILIZER SIZING
# ======================================================

def horizontal_stabilizer(

    X_aftcg,
    l_fuselage,
    MAC,
    S_w,
    c,
    AR,
    X_h,
    V_h,
    AR_h,
    Sweep_h_LE,
    Taper_h,
    t_c_h

):

    # Area
    S_h = (
        V_h * S_w * c
    ) / (
        X_h - X_aftcg
    )

    # Span
    b_h = np.sqrt(AR_h * S_h)

    # Chords
    c_root = (
        2 * S_h
    ) / (
        b_h * (1 + Taper_h)
    )

    c_tip = Taper_h * c_root

    # Quarter-chord sweep
    Sweep_h_qc = (
        Sweep_h_LE
        - np.arctan(
            (c_root - c_tip) / b_h
        ) * (180 / np.pi)
    )

    return (
        S_h,
        b_h,
        c_root,
        c_tip,
        Sweep_h_qc
    )