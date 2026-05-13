import numpy as np

# ==========================================
# CONFIGURATION BLOCK (Default Fallbacks)
# ==========================================
X_aft_cg = 1.0        # m from datum
l_fuselage = 8.8      # m

S_w = 23.0            # m^2
MAC_w = 1.5           # m
AR_w = 11.0           # Aspect Ratio

X_h = 3.0             # m from datum to tail quarter-chord
V_h = 0.786           # Volume coefficient (Roskam average)
AR_h = 4.0            # Aspect Ratio
Sweep_h_LE = 25.0     # Leading edge sweep (degrees)
Taper_h = 0.4         # Taper ratio
t_c_h: float = 0.12   # Thickness-to-chord ratio

X_v_factor = 0.9      # Fraction of fuselage length if X_v is estimated
V_v = 0.06275         # Volume coefficient (Roskam average)
AR_v = 1.5            # Aspect Ratio
Sweep_v_LE = 25.0     # Leading edge sweep (degrees)
Taper_v = 0.4         # Taper ratio
t_c_v: float = 0.15   # Thickness-to-chord ratio


def calculate_empennage(wing_data=None):
    """Computes horizontal and vertical tail sizing parameters."""
    # Allow dynamic overrides from wing sizing if provided
    local_S_w = wing_data.S if (wing_data and hasattr(wing_data, 'S')) else S_w
    local_MAC = wing_data.MAC if (wing_data and hasattr(wing_data, 'MAC')) else MAC_w
    local_AR = wing_data.AR if (wing_data and hasattr(wing_data, 'AR')) else AR_w
    
    # Derive main wing span dynamically
    b_w = np.sqrt(local_AR * local_S_w)
    
    # --- HORIZONTAL STABILIZER SIZING ---
    l_h = X_h - X_aft_cg
    S_h = (V_h * local_S_w * local_MAC) / l_h
    b_h = np.sqrt(AR_h * S_h)
    
    c_root_h = (2 * S_h) / (b_h * (1 + Taper_h))
    c_tip_h = Taper_h * c_root_h
    
    tan_sweep_le_h = np.tan(np.radians(Sweep_h_LE))
    tan_sweep_qc_h = tan_sweep_le_h - ((c_root_h - c_tip_h) / b_h)
    sweep_qc_h = np.degrees(np.arctan(tan_sweep_qc_h))

    # --- VERTICAL STABILIZER SIZING ---
    X_v = l_fuselage * X_v_factor
    l_v = X_v - X_aft_cg
    
    S_v = (V_v * local_S_w * b_w) / l_v
    b_v = np.sqrt(AR_v * S_v)
    
    c_root_v = (2 * S_v) / (b_v * (1 + Taper_v))
    c_tip_v = Taper_v * c_root_v
    
    tan_sweep_le_v = np.tan(np.radians(Sweep_v_LE))
    tan_sweep_qc_v = tan_sweep_le_v - ((c_root_v - c_tip_v) / (2 * b_v))
    sweep_qc_v = np.degrees(np.arctan(tan_sweep_qc_v))

    return {
        "b_w": b_w,
        "h_tail": {
            "S_h": S_h, "b_h": b_h, "c_root": c_root_h, 
            "c_tip": c_tip_h, "sweep_qc": sweep_qc_h, "l_h": l_h
        },
        "v_tail": {
            "S_v": S_v, "b_v": b_v, "c_root": c_root_v, 
            "c_tip": c_tip_v, "sweep_qc": sweep_qc_v, "l_v": l_v, "X_v": X_v
        }
    }

def print_empennage_results(res: dict):
    """Outputs cleanly formatted engineering console logs."""
    print("\n" + "-" * 55)
    print("HORIZONTAL STABILIZER RESULTS:")
    print("-" * 55)
    print(f"  Moment Arm (l_h):         {res['h_tail']['l_h']:.2f} m")
    print(f"  Area (S_h):               {res['h_tail']['S_h']:.2f} m²")
    print(f"  Span (b_h):               {res['h_tail']['b_h']:.2f} m")
    print(f"  Root Chord (c_root):      {res['h_tail']['c_root']:.2f} m")
    print(f"  Tip Chord (c_tip):        {res['h_tail']['c_tip']:.2f} m")
    print(f"  Quarter-Chord Sweep:      {res['h_tail']['sweep_qc']:.1f}°")
    
    print("\n" + "-" * 55)
    print("VERTICAL STABILIZER RESULTS:")
    print("-" * 55)
    print(f"  Station Arm (X_v):        {res['v_tail']['X_v']:.2f} m")
    print(f"  Moment Arm (l_v):         {res['v_tail']['l_v']:.2f} m")
    print(f"  Area (S_v):               {res['v_tail']['S_v']:.2f} m²")
    print(f"  Span (b_v):               {res['v_tail']['b_v']:.2f} m")
    print(f"  Root Chord (c_root):      {res['v_tail']['c_root']:.2f} m")
    print(f"  Tip Chord (c_tip):        {res['v_tail']['c_tip']:.2f} m")
    print(f"  Quarter-Chord Sweep:      {res['v_tail']['sweep_qc']:.1f}°")
    print("-" * 55 + "\n")