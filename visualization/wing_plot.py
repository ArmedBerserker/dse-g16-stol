import matplotlib.pyplot as plt
import numpy as np


# ==========================================================
# WING PLANFORM PLOT
# ==========================================================
#
# PURPOSE:
# - Visualize trapezoidal wing planform
# - Plot:
#   * Leading edge
#   * Trailing edge
#   * Quarter-chord line
#   * MAC
#   * Front spar
#   * Rear spar
#
# INPUT:
# - Wing object from wing_planform.py
#
# ==========================================================


def plot_wing_planform(wing):

    # ======================================================
    # BASIC GEOMETRY
    # ======================================================

    half_span = wing.span / 2

    # Convert sweep angles to radians
    sweep_LE_rad = np.radians(
        wing.sweep_LE_deg
    )

    sweep_TE_rad = np.radians(
        wing.sweep_TE_deg
    )

    # ======================================================
    # ROOT GEOMETRY
    # ======================================================

    x_root_LE = 0.0
    x_root_TE = wing.c_root

    # ======================================================
    # TIP GEOMETRY
    # ======================================================

    x_tip_LE = (
        half_span
        * np.tan(sweep_LE_rad)
    )

    x_tip_TE = (
        x_tip_LE
        + wing.c_tip
    )

    # ======================================================
    # PLANFORM OUTLINE
    # ======================================================

    # Right wing
    x_outline = [
        x_root_LE,
        x_tip_LE,
        x_tip_TE,
        x_root_TE,
        x_root_LE
    ]

    y_outline = [
        0,
        half_span,
        half_span,
        0,
        0
    ]

    # Left wing (mirrored)
    y_outline_mirror = [
        -y for y in y_outline
    ]

    # ======================================================
    # CREATE FIGURE
    # ======================================================

    plt.figure(figsize=(12, 6))

    # ======================================================
    # PLOT WING OUTLINE
    # ======================================================

    plt.plot(
        x_outline,
        y_outline,
        linewidth=2,
        label="Wing outline"
    )

    plt.plot(
        x_outline,
        y_outline_mirror,
        linewidth=2
    )

    # ======================================================
    # QUARTER-CHORD LINE
    # ======================================================

    x_qc_root = 0.25 * wing.c_root

    x_qc_tip = (
        x_tip_LE
        + 0.25 * wing.c_tip
    )

    plt.plot(
        [x_qc_root, x_qc_tip],
        [0, half_span],
        linestyle="--",
        linewidth=2,
        label="Quarter-chord"
    )

    plt.plot(
        [x_qc_root, x_qc_tip],
        [0, -half_span],
        linestyle="--",
        linewidth=2
    )

    # ======================================================
    # FRONT SPAR
    # ======================================================

    x_fs_root = (
        wing.front_spar_frac
        * wing.c_root
    )

    x_fs_tip = (
        x_tip_LE
        + wing.front_spar_frac
        * wing.c_tip
    )

    plt.plot(
        [x_fs_root, x_fs_tip],
        [0, half_span],
        linewidth=2,
        label="Front spar"
    )

    plt.plot(
        [x_fs_root, x_fs_tip],
        [0, -half_span],
        linewidth=2
    )

    # ======================================================
    # REAR SPAR
    # ======================================================

    x_rs_root = (
        wing.rear_spar_frac
        * wing.c_root
    )

    x_rs_tip = (
        x_tip_LE
        + wing.rear_spar_frac
        * wing.c_tip
    )

    plt.plot(
        [x_rs_root, x_rs_tip],
        [0, half_span],
        linewidth=2,
        label="Rear spar"
    )

    plt.plot(
        [x_rs_root, x_rs_tip],
        [0, -half_span],
        linewidth=2
    )

    # ======================================================
    # MAC
    # ======================================================

    plt.plot(
        [wing.x_LE_MAC,
         wing.x_LE_MAC + wing.MAC],

        [wing.y_MAC,
         wing.y_MAC],

        linewidth=4,
        label="MAC"
    )

    plt.plot(
        [wing.x_LE_MAC,
         wing.x_LE_MAC + wing.MAC],

        [-wing.y_MAC,
         -wing.y_MAC],

        linewidth=4
    )

    # ======================================================
    # AERODYNAMIC CENTER
    #
    # Approximate at 25% MAC
    # ======================================================

    x_ac = (
        wing.x_LE_MAC
        + 0.25 * wing.MAC
    )

    plt.scatter(
        x_ac,
        wing.y_MAC,
        s=80,
        label="Aerodynamic center"
    )

    plt.scatter(
        x_ac,
        -wing.y_MAC,
        s=80
    )

    # ======================================================
    # PLOT SETTINGS
    # ======================================================

    plt.axis("equal")

    plt.grid(True)

    plt.xlabel("x-position [m]")

    plt.ylabel("Spanwise position y [m]")

    plt.title("Wing Planform")

    plt.legend()

    plt.tight_layout()

    plt.show()