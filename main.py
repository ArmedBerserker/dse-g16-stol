import sys
import matplotlib.pyplot as plt

# ==============================================================================
# 1. MODULE IMPORTS
# ==============================================================================
# Wing geometry package imports
from geometry.wing_planform import (
    load_wing_from_yaml,
    size_wing_planform,
    print_wing_geometry
)
from visualization.wing_plot import plot_wing_planform

# Empennage package imports
from empennage.tail_sizing import (
    calculate_empennage,
    print_empennage_results,
    # Fallback variables imported for downstream plotter support
    X_aft_cg, l_fuselage, MAC_w, S_w, AR_w,
    X_h, V_h, AR_h, Sweep_h_LE, Taper_h, t_c_h
)

# Empennage visualization package imports
# Using the updated modular tail plotter logic
try:
    from visualization.tail_plot import plot_empennage
except ImportError:
    print("Warning: Could not load visualization.tail_plot module.")
    plot_empennage = lambda cfg: None

# Optional Legacy Plotter Support (Maintained for assignment checks)
try:
    from visualization.plot_horizontal_tail import plot_horizontal_tail
except ImportError:
    plot_horizontal_tail = lambda *args: None


# ==============================================================================
# 2. MASTER EXECUTION PIPELINE
# ==============================================================================
def main():
    print("\n" + "="*70)
    print("INITIALIZING GROUP 16 DSE CONCEPTUAL DESIGN PIPELINE")
    print("="*70)
    
    # --------------------------------------------------------------------------
    # PHASE 1: WING CONFIGURATION & SIZING
    # --------------------------------------------------------------------------
    wing_config_path = "yamls/wing_v2.yaml"
    print(f"\n[Phase 1] Loading Main Wing Configuration from: '{wing_config_path}'")
    
    try:
        # Load and execute sizing logic
        wing = load_wing_from_yaml(wing_config_path)
        wing = size_wing_planform(wing)
        
        # Display numerical console log
        print_wing_geometry(wing)
        
        # Prepare Matplotlib visualization (Muted plt.show internally)
        plot_wing_planform(wing)
        
    except Exception as e:
        print(f"\n>>> NOTICE: Wing sizing pipeline bypassed or failed: {e}")
        print(">>> Proceeding to Empennage sizing using local baseline fallbacks.")
        wing = None

    # --------------------------------------------------------------------------
    # PHASE 2: EMPENNAGE SIZING
    # --------------------------------------------------------------------------
    print("\n[Phase 2] Executing Coupled Empennage Sizing Physics...")
    
    # Execute calculations (passing dynamic wing data if available)
    emp_results = calculate_empennage(wing_data=wing)
    
    # Output aligned console diagnostics
    print_empennage_results(emp_results)

    # --------------------------------------------------------------------------
    # PHASE 3: EMPENNAGE VISUALIZATION
    # --------------------------------------------------------------------------
    print("[Phase 3] Generating Scaled Visual Planform Vectors...")
    
    # 1. Execute the primary dual-surface dictionary plotter
    plot_empennage(emp_results)
    
    # 2. Execute legacy assignment plotter requirement (if applicable)
    h_data = emp_results["h_tail"]
    plot_horizontal_tail(
        X_aft_cg,
        l_fuselage,
        wing.MAC if (wing and hasattr(wing, 'MAC')) else MAC_w,
        wing.S if (wing and hasattr(wing, 'S')) else S_w,
        h_data["c_root"],  # Pass dynamically derived root chord
        wing.AR if (wing and hasattr(wing, 'AR')) else AR_w,
        X_h,
        V_h,
        AR_h,
        Sweep_h_LE,
        Taper_h,
        t_c_h
    )

    # --------------------------------------------------------------------------
    # PHASE 4: RENDER ALL FIGURES SIMULTANEOUSLY
    # --------------------------------------------------------------------------
    print("\n>>> Pipeline execution complete. Launching Matplotlib figures...")
    print(">>> (Close all plot windows to terminate the script execution)\n")
    
    # Master render call triggers all queued subplots at once
    plt.show()


# ==============================================================================
# ENTRY POINT TRIGGER
# ==============================================================================
if __name__ == "__main__":
    main()