import numpy as np
from dataclasses import dataclass




## Working
# Class 1 yields us the wing area and aspect ratio. We can use these to compute the span and root chord. We can then compute the tip chord from the taper ratio. With the root and tip chords, we can compute the mean aerodynamic chord (MAC) and its location. Finally, we can compute the quarter-chord and trailing edge sweep angles from the leading edge sweep.
# From reference AC in Roskam we obtain the following parameters: sweep & taper ratio
# The following parameters are then calculated:
# c_root, c_tip, span, MAC, y_MAC, x_LE_MAC, sweep_c4, sweep_TE
# 
# Select t/c ratio --> choose airfoil (Cd, Cldes, Mcrit, Cmc/4)
# Decide on type, size and location of control devices
# Draw the front and rear spar lines in the wing planform
# Compute wing fuel volume
# Decide on the wing dihedral angle
# Decide on the wing incidence angle
# ==========================================================
# DATA STRUCTURE
# ==========================================================

@dataclass
class WingPlanform:
    # Inputs
    S: float                 # Wing area [m^2]
    AR: float                # Aspect ratio [-]
    taper_ratio: float       # lambda = ct/cr
    sweep_LE_deg: float      # Leading edge sweep [deg]

    # Computed geometry
    span: float = None
    c_root: float = None
    c_tip: float = None
    MAC: float = None
    y_MAC: float = None
    x_LE_MAC: float = None

    # Sweep angles
    sweep_c4_deg: float = None
    sweep_TE_deg: float = None


# ==========================================================
# MAIN SIZING FUNCTION
# ==========================================================

def size_wing_planform(wing: WingPlanform):

    # ------------------------------------------------------
    # 1. Wing span from AR = b² / S
    # ------------------------------------------------------
    wing.span = np.sqrt(wing.AR * wing.S)

    # ------------------------------------------------------
    # 2. Root chord from trapezoidal wing relation
    #
    # S = (b/2) * (cr + ct)
    # ct = lambda * cr
    # ------------------------------------------------------
    wing.c_root = (
        2 * wing.S /
        (wing.span * (1 + wing.taper_ratio))
    )

    # Tip chord
    wing.c_tip = wing.taper_ratio * wing.c_root

    # ------------------------------------------------------
    # 3. Mean Aerodynamic Chord (MAC)
    # 2435
    # Roskam standard trapezoidal formula
    # ------------------------------------------------------
    wing.MAC = (
        (2 / 3)
        * wing.c_root
        * (
            (1 + wing.taper_ratio + wing.taper_ratio**2)
            / (1 + wing.taper_ratio)
        )
    )

    # ------------------------------------------------------
    # 4. Spanwise MAC location
    # ------------------------------------------------------
    wing.y_MAC = (
        wing.span / 6
        * (
            (1 + 2 * wing.taper_ratio)
            / (1 + wing.taper_ratio)
        )
    )

    # ------------------------------------------------------
    # 5. Leading edge x-location of MAC
    # ------------------------------------------------------
    sweep_LE_rad = np.radians(wing.sweep_LE_deg)

    wing.x_LE_MAC = (
        wing.y_MAC * np.tan(sweep_LE_rad)
    )

    # ------------------------------------------------------
    # 6. Quarter-chord sweep
    # Roskam trapezoidal relation
    # ------------------------------------------------------
    wing.sweep_c4_deg = np.degrees(
        np.arctan(
            np.tan(sweep_LE_rad)
            - (
                4 / wing.AR
                * (
                    (0.25 * (1 - wing.taper_ratio))
                    / (1 + wing.taper_ratio)
                )
            )
        )
    )

    # ------------------------------------------------------
    # 7. Trailing edge sweep
    # ------------------------------------------------------
    wing.sweep_TE_deg = np.degrees(
        np.arctan(
            np.tan(sweep_LE_rad)
            - (
                4 / wing.AR
                * (
                    (1 - wing.taper_ratio)
                    / (1 + wing.taper_ratio)
                )
            )
        )
    )

    return wing


# ==========================================================
# EXAMPLE
# ==========================================================

if __name__ == "__main__":

    wing = WingPlanform(
        S=21.0,  ## From Class 1             # m²
        AR=9.0, ## From Class 1
        taper_ratio=0.5,
        sweep_LE_deg=0
    )

    wing = size_wing_planform(wing)

    print("========== WING PLANFORM ==========")
    print(f"Wing area S           : {wing.S:.2f} m²")
    print(f"Aspect ratio AR       : {wing.AR:.2f}")
    print(f"Wingspan b            : {wing.span:.2f} m")
    print(f"Root chord cr         : {wing.c_root:.2f} m")
    print(f"Tip chord ct          : {wing.c_tip:.2f} m")
    print(f"MAC                   : {wing.MAC:.2f} m")
    print(f"y_MAC                 : {wing.y_MAC:.2f} m")
    print(f"x_LE_MAC              : {wing.x_LE_MAC:.2f} m")
    print(f"LE sweep              : {wing.sweep_LE_deg:.2f} deg")
    print(f"Quarter-chord sweep   : {wing.sweep_c4_deg:.2f} deg")
    print(f"Trailing edge sweep   : {wing.sweep_TE_deg:.2f} deg")