import yaml
import numpy as np
from dataclasses import dataclass


# ==========================================================
# WING DATA CLASS
# ==========================================================

@dataclass
class Wing:

    # ------------------------------------------------------
    # INPUT PARAMETERS
    # ------------------------------------------------------

    S: float
    AR: float

    taper_ratio: float

    sweep_c4_deg: float

    tc_ratio: float

    dihedral_deg: float

    incidence_deg: float

    twist_deg: float

    front_spar_frac: float
    rear_spar_frac: float

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
# YAML IMPORT
# ==========================================================

def load_wing_from_yaml(filepath):

    # ------------------------------------------------------
    # LOAD YAML FILE
    # ------------------------------------------------------

    with open(filepath, "r") as file:
        data = yaml.safe_load(file)

    # ------------------------------------------------------
    # ACCESS MAIN WING DICTIONARY
    # ------------------------------------------------------

    wing_data = data["wing"]

    # ------------------------------------------------------
    # ACCESS SUB-DICTIONARIES
    # ------------------------------------------------------

    geometry = wing_data["geometry"]

    airfoil = wing_data["airfoil"]

    settings = wing_data["settings"]

    structure = wing_data["structure"]

    aerodynamic = wing_data["aerodynamic"]

    correction = wing_data["correction_factors"]

    # ------------------------------------------------------
    # CREATE WING OBJECT
    # ------------------------------------------------------

    wing = Wing(

        # --------------------------------------------------
        # GEOMETRY
        # --------------------------------------------------

        S=geometry["area"],

        AR=geometry["aspect_ratio"],

        taper_ratio=geometry["taper_ratio"],

        sweep_c4_deg=geometry["sweep_c4_deg"],

        # --------------------------------------------------
        # AIRFOIL
        # --------------------------------------------------

        tc_ratio=airfoil["tc_ratio"],

        # --------------------------------------------------
        # SETTINGS
        # --------------------------------------------------

        dihedral_deg=settings["dihedral_deg"],

        incidence_deg=settings["incidence_deg"],

        twist_deg=settings["twist_deg"],

        # --------------------------------------------------
        # STRUCTURE
        # --------------------------------------------------

        front_spar_frac=structure["front_spar"],

        rear_spar_frac=structure["rear_spar"]

        # --------------------------------------------------
        # OPTIONAL FUTURE PARAMETERS
        # --------------------------------------------------
        #
        # c_f=aerodynamic["c_f"],
        #
        # ld=aerodynamic["ld"],
        #
        # phi=correction["phi"],
        #
        # psi=correction["psi"]
        #
        # Add these later if included in Wing dataclass
        # --------------------------------------------------
    )

    return wing


# ==========================================================
# PLANFORM SIZING
# ==========================================================

def size_wing_planform(wing):

    # ------------------------------------------------------
    # 1. SPAN
    #
    # AR = b^2 / S
    # ------------------------------------------------------

    wing.span = np.sqrt(wing.AR * wing.S)

    # ------------------------------------------------------
    # 2. ROOT CHORD
    #
    # S = (b/2)(cr + ct)
    # ct = lambda * cr
    # ------------------------------------------------------

    wing.c_root = (
        2 * wing.S
        /
        (
            wing.span
            * (1 + wing.taper_ratio)
        )
    )

    # ------------------------------------------------------
    # 3. TIP CHORD
    # ------------------------------------------------------

    wing.c_tip = (
        wing.taper_ratio
        * wing.c_root
    )

    # ------------------------------------------------------
    # 4. LEADING EDGE SWEEP
    #
    # Roskam relation
    # ------------------------------------------------------

    sweep_c4_rad = np.radians(
        wing.sweep_c4_deg
    )

    sweep_LE_rad = np.arctan(

        np.tan(sweep_c4_rad)

        +

        (
            1 / wing.AR
        )

        *

        (
            (1 - wing.taper_ratio)
            /
            (1 + wing.taper_ratio)
        )
    )

    wing.sweep_LE_deg = np.degrees(
        sweep_LE_rad
    )

    # ------------------------------------------------------
    # 5. TRAILING EDGE SWEEP
    # ------------------------------------------------------

    wing.sweep_TE_deg = np.degrees(

        np.arctan(

            np.tan(sweep_LE_rad)

            -

            (
                2 / wing.AR
            )

            *

            (
                (1 - wing.taper_ratio)
                /
                (1 + wing.taper_ratio)
            )
        )
    )

    # ------------------------------------------------------
    # 6. MEAN AERODYNAMIC CHORD
    # ------------------------------------------------------

    wing.MAC = (

        (2 / 3)

        * wing.c_root

        *

        (
            (
                1
                + wing.taper_ratio
                + wing.taper_ratio**2
            )

            /

            (
                1 + wing.taper_ratio
            )
        )
    )

    # ------------------------------------------------------
    # 7. SPANWISE MAC LOCATION
    # ------------------------------------------------------

    wing.y_MAC = (

        wing.span / 6

        *

        (
            (
                1
                + 2 * wing.taper_ratio
            )

            /

            (
                1 + wing.taper_ratio
            )
        )
    )

    # ------------------------------------------------------
    # 8. MAC x-LOCATION
    # ------------------------------------------------------

    wing.x_LE_MAC = (

        wing.y_MAC

        * np.tan(sweep_LE_rad)
    )

    return wing


# ==========================================================
# PRINT RESULTS
# ==========================================================

def print_wing_geometry(wing):

    print("\n========== WING PLANFORM ==========\n")

    print(f"S                      : {wing.S:.2f} m²")
    print(f"AR                     : {wing.AR:.2f}")

    print(f"\nSpan                   : {wing.span:.2f} m")

    print(f"\nRoot chord             : {wing.c_root:.2f} m")
    print(f"Tip chord              : {wing.c_tip:.2f} m")

    print(f"\nMAC                    : {wing.MAC:.2f} m")

    print(f"\ny_MAC                  : {wing.y_MAC:.2f} m")
    print(f"x_LE_MAC               : {wing.x_LE_MAC:.2f} m")

    print(f"\nQuarter-chord sweep    : {wing.sweep_c4_deg:.2f} deg")
    print(f"Leading edge sweep     : {wing.sweep_LE_deg:.2f} deg")
    print(f"Trailing edge sweep    : {wing.sweep_TE_deg:.2f} deg")

    print(f"\nTaper ratio            : {wing.taper_ratio:.2f}")

    print(f"\nt/c ratio              : {wing.tc_ratio:.2f}")

    print(f"\nDihedral               : {wing.dihedral_deg:.2f} deg")

    print(f"Incidence              : {wing.incidence_deg:.2f} deg")

    print(f"Twist                  : {wing.twist_deg:.2f} deg")