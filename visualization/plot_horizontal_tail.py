import numpy as np
import matplotlib.pyplot as plt


import numpy as np
import matplotlib.pyplot as plt

from empennage.tail_sizing import (
    horizontal_stabilizer
)


def plot_horizontal_tail(

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

    # --------------------------------------------------
    # GET GEOMETRY FROM SIZING FUNCTION
    # --------------------------------------------------

    S_h, b_h, c_root, c_tip, Sweep_h_qc = horizontal_stabilizer(

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
    )

    # --------------------------------------------------
    # PLANFORM GEOMETRY
    # --------------------------------------------------

    half_span = b_h / 2

    x_le_root = 0

    x_le_tip = (
        half_span
        * np.tan(np.radians(Sweep_h_LE))
    )

    # Planform coordinates
    x = [
        x_le_root,
        x_le_tip,
        x_le_tip + c_tip,
        x_le_root + c_root,
        x_le_root
    ]

    y = [
        0,
        half_span,
        half_span,
        0,
        0
    ]

    # --------------------------------------------------
    # PLOT
    # --------------------------------------------------

    plt.figure(figsize=(8, 4))

    plt.plot(x, y)
    plt.fill(x, y, alpha=0.3)

    plt.axis("equal")

    plt.xlabel("x [m]")
    plt.ylabel("y [m]")

    plt.title("Horizontal Stabilizer Planform")

    plt.grid(True)

    plt.show()