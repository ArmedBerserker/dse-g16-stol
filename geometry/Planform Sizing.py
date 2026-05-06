import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass


# ==========================================================
# ROSKAM CLASS I WING PLANFORM SIZING
# ==========================================================
#
# INPUTS FROM CLASS I:
# - Wing area S
# - Aspect ratio AR
#
# INPUTS FROM REFERENCE AIRCRAFT / DESIGN CHOICES:
# - Quarter-chord sweep
# - Taper ratio
#
# THIS SCRIPT COMPUTES:
# - span
# - root chord
# - tip chord
# - MAC
# - MAC location
# - LE and TE sweep
#
# ADDITIONAL FEATURES:
# - wing planform plotting
# - quarter-chord visualization
# - MAC visualization
# - front and rear spar visualization
#
# NEXT DESIGN STEPS:
# - Select t/c ratio
# - Select airfoil
# - Compute fuel volume
# - High-lift device sizing
# - Dihedral angle
# - Incidence angle
#
# ==========================================================










# ==========================================================
# DATA STRUCTURE
# ==========================================================

@dataclass
class WingPlanform:

    # ------------------------------------------------------
    # INPUT PARAMETERS
    # ------------------------------------------------------

    S: float                         # Wing area [m^2]
    AR: float                        # Aspect ratio [-]

    taper_ratio: float               # lambda = ct/cr

    sweep_c4_deg: float              # Quarter-chord sweep [deg]

    dihedral_deg: float = 0.0        # Optional
    incidence_deg: float = 0.0       # Optional

    # Spar locations as chord fractions
    front_spar_frac: float = 0.15
    rear_spar_frac: float = 0.65

    # ------------------------------------------------------
    # COMPUTED PARAMETERS
    # ------------------------------------------------------

    span: float = None

    c_root: float = None
    c_tip: float = None

    MAC: float = None

    y_MAC: float = None
    x_LE_MAC: float = None

    sweep_LE_deg: float = None
    sweep_TE_deg: float = None


# ==========================================================
# MAIN SIZING FUNCTION
# ==========================================================

def size_wing_planform(wing: WingPlanform):

    # ------------------------------------------------------
    # 1. WING SPAN
    # ------------------------------------------------------

    wing.span = np.sqrt(wing.AR * wing.S)

    # ------------------------------------------------------
    # 2. ROOT + TIP CHORD
    # ------------------------------------------------------

    wing.c_root = (
        2 * wing.S
        /
        (wing.span * (1 + wing.taper_ratio))
    )

    wing.c_tip = wing.taper_ratio * wing.c_root

    # ------------------------------------------------------
    # 3. LEADING EDGE SWEEP
    #
    # Roskam trapezoidal relation:
    #
    # tan(LE) =
    # tan(c/4)
    # + (1/AR)*(1-lambda)/(1+lambda)
    #
    # ------------------------------------------------------

    sweep_c4_rad = np.radians(wing.sweep_c4_deg)

    sweep_LE_rad = np.arctan(
        np.tan(sweep_c4_rad)
        +
        (
            (1 / wing.AR)
            *
            (
                (1 - wing.taper_ratio)
                /
                (1 + wing.taper_ratio)
            )
        )
    )

    wing.sweep_LE_deg = np.degrees(sweep_LE_rad)

    # ------------------------------------------------------
    # 4. TRAILING EDGE SWEEP
    # ------------------------------------------------------

    wing.sweep_TE_deg = np.degrees(
        np.arctan(
            np.tan(sweep_LE_rad)
            -
            (
                (2 / wing.AR)
                *
                (
                    (1 - wing.taper_ratio)
                    /
                    (1 + wing.taper_ratio)
                )
            )
        )
    )

    # ------------------------------------------------------
    # 5. MEAN AERODYNAMIC CHORD
    # ------------------------------------------------------

    wing.MAC = (
        (2 / 3)
        * wing.c_root
        * (
            (1 + wing.taper_ratio + wing.taper_ratio**2)
            /
            (1 + wing.taper_ratio)
        )
    )

    # ------------------------------------------------------
    # 6. SPANWISE MAC LOCATION
    # ------------------------------------------------------

    wing.y_MAC = (
        wing.span / 6
        *
        (
            (1 + 2 * wing.taper_ratio)
            /
            (1 + wing.taper_ratio)
        )
    )

    # ------------------------------------------------------
    # 7. x-LOCATION OF MAC
    # ------------------------------------------------------

    wing.x_LE_MAC = (
        wing.y_MAC
        *
        np.tan(sweep_LE_rad)
    )

    return wing


# ==========================================================
# PLOT FUNCTION
# ==========================================================

def plot_wing_planform(wing):

    half_span = wing.span / 2

    sweep_LE_rad = np.radians(wing.sweep_LE_deg)
    sweep_TE_rad = np.radians(wing.sweep_TE_deg)

    # ------------------------------------------------------
    # ROOT GEOMETRY
    # ------------------------------------------------------

    x_root_LE = 0
    x_root_TE = wing.c_root

    # ------------------------------------------------------
    # TIP GEOMETRY
    # ------------------------------------------------------

    x_tip_LE = half_span * np.tan(sweep_LE_rad)
    x_tip_TE = x_tip_LE + wing.c_tip

    # ------------------------------------------------------
    # PLANFORM POLYGON
    # ------------------------------------------------------

    x_coords = [
        x_root_LE,
        x_tip_LE,
        x_tip_TE,
        x_root_TE,
        x_root_LE
    ]

    y_coords = [
        0,
        half_span,
        half_span,
        0,
        0
    ]

    # Mirror wing
    y_coords_mirror = [-y for y in y_coords]

    # ------------------------------------------------------
    # PLOTTING
    # ------------------------------------------------------

    plt.figure(figsize=(12, 6))

    # Right wing
    plt.plot(x_coords, y_coords, linewidth=2)

    # Left wing
    plt.plot(x_coords, y_coords_mirror, linewidth=2)

    # ------------------------------------------------------
    # QUARTER CHORD LINE
    # ------------------------------------------------------

    x_qc_root = 0.25 * wing.c_root
    x_qc_tip = x_tip_LE + 0.25 * wing.c_tip

    plt.plot(
        [x_qc_root, x_qc_tip],
        [0, half_span],
        linestyle='--',
        label='Quarter-chord'
    )

    plt.plot(
        [x_qc_root, x_qc_tip],
        [0, -half_span],
        linestyle='--'
    )

    # ------------------------------------------------------
    # FRONT SPAR
    # ------------------------------------------------------

    x_fs_root = wing.front_spar_frac * wing.c_root
    x_fs_tip = x_tip_LE + wing.front_spar_frac * wing.c_tip

    plt.plot(
        [x_fs_root, x_fs_tip],
        [0, half_span],
        label='Front spar'
    )

    plt.plot(
        [x_fs_root, x_fs_tip],
        [0, -half_span]
    )

    # ------------------------------------------------------
    # REAR SPAR
    # ------------------------------------------------------

    x_rs_root = wing.rear_spar_frac * wing.c_root
    x_rs_tip = x_tip_LE + wing.rear_spar_frac * wing.c_tip

    plt.plot(
        [x_rs_root, x_rs_tip],
        [0, half_span],
        label='Rear spar'
    )

    plt.plot(
        [x_rs_root, x_rs_tip],
        [0, -half_span]
    )

    # ------------------------------------------------------
    # MAC VISUALIZATION
    # ------------------------------------------------------

    plt.plot(
        [wing.x_LE_MAC, wing.x_LE_MAC + wing.MAC],
        [wing.y_MAC, wing.y_MAC],
        linewidth=4,
        label='MAC'
    )

    plt.plot(
        [wing.x_LE_MAC, wing.x_LE_MAC + wing.MAC],
        [-wing.y_MAC, -wing.y_MAC],
        linewidth=4
    )

    # ------------------------------------------------------
    # FIGURE SETTINGS
    # ------------------------------------------------------

    plt.axis('equal')

    plt.xlabel('x [m]')
    plt.ylabel('y [m]')

    plt.title('Wing Planform')

    plt.grid(True)

    plt.legend()

    plt.show()


# ==========================================================
# EXAMPLE
# ==========================================================

if __name__ == "__main__":

    wing = WingPlanform(

        # From Class I sizing
        S=21.0,
        AR=9.0,

        # From reference aircraft / design choice
        taper_ratio=0.5,
        sweep_c4_deg=2.0,

        # Structural assumptions
        front_spar_frac=0.15,
        rear_spar_frac=0.65
    )

    wing = size_wing_planform(wing)

    # ------------------------------------------------------
    # PRINT RESULTS
    # ------------------------------------------------------

    print("\n========== WING PLANFORM ==========\n")

    print(f"Wing area S               : {wing.S:.2f} m²")
    print(f"Aspect ratio AR           : {wing.AR:.2f}")

    print(f"\nSpan b                    : {wing.span:.2f} m")

    print(f"\nRoot chord cr             : {wing.c_root:.2f} m")
    print(f"Tip chord ct              : {wing.c_tip:.2f} m")

    print(f"\nTaper ratio               : {wing.taper_ratio:.2f}")

    print(f"\nQuarter-chord sweep       : {wing.sweep_c4_deg:.2f} deg")
    print(f"Leading-edge sweep        : {wing.sweep_LE_deg:.2f} deg")
    print(f"Trailing-edge sweep       : {wing.sweep_TE_deg:.2f} deg")

    print(f"\nMAC                       : {wing.MAC:.2f} m")
    print(f"y_MAC                     : {wing.y_MAC:.2f} m")
    print(f"x_LE_MAC                  : {wing.x_LE_MAC:.2f} m")

    # ------------------------------------------------------
    # DRAW PLANFORM
    # ------------------------------------------------------

    plot_wing_planform(wing)