
from geometry.wing_planform import (
    load_wing_from_yaml,
    size_wing_planform,
    print_wing_geometry
)

from visualization.wing_plot import plot_wing_planform
from visualization.plot_horizontal_tail import (
    plot_horizontal_tail
)
from empennage.horizontal_stabilizer import (
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


# ==========================================================
# MAIN
# ==========================================================

def main():

    # ------------------------------------------------------
    # LOAD WING CONFIGURATION
    # ------------------------------------------------------

    wing = load_wing_from_yaml(
        "yamls/wing_v2.yaml"
    )

    # ------------------------------------------------------
    # COMPUTE PLANFORM GEOMETRY
    # ------------------------------------------------------

    wing = size_wing_planform(wing)

    # ------------------------------------------------------
    # PRINT RESULTS
    # ------------------------------------------------------

    print_wing_geometry(wing)
    plot_wing_planform(wing)
    # ------------------------------------------------------
    # PLOT HORIZONTAL TAIL
    # ------------------------------------------------------

    plot_horizontal_tail(

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


# ==========================================================
# RUN SCRIPT
# ==========================================================

if __name__ == "__main__":

    main()