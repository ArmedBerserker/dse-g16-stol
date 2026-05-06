
from geometry.wing_planform import (
    load_wing_from_yaml,
    size_wing_planform,
    print_wing_geometry
)

from visualization.wing_plot import plot_wing_planform


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


# ==========================================================
# RUN SCRIPT
# ==========================================================

if __name__ == "__main__":

    main()