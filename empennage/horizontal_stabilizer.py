import numpy as np

# ==========================================
# CONFIGURATION BLOCK (Adjust these later)
# ==========================================

# --- Main Aircraft & Loading Parameters ---
X_aft_cg = 1.0        # m from datum
l_fuselage = 8.8      # m

# --- Main Wing Geometry ---
S_w = 23.0            # m^2
MAC_w = 1.5           # m
AR_w = 11.0           # Aspect Ratio

# --- Horizontal Stabilizer Parameters ---
X_h = 3.0             # m from datum to tail quarter-chord
V_h = 0.786           # Volume coefficient (Roskam average)
AR_h = 4.0            # Aspect Ratio
Sweep_h_LE = 25.0     # Leading edge sweep (degrees)
Taper_h = 0.4         # Taper ratio
t_c_h: float = 0.12   # Thickness-to-chord ratio

# --- Vertical Stabilizer Parameters ---
X_v_factor = 0.9      # Fraction of fuselage length if X_v is estimated
V_v = 0.06275         # Volume coefficient (Roskam average)
AR_v = 1.5            # Aspect Ratio
Sweep_v_LE = 25.0     # Leading edge sweep (degrees)
Taper_v = 0.4         # Taper ratio
t_c_v: float = 0.15   # Thickness-to-chord ratio


# ==========================================
# EMPENNAGE CALCULATIONS
# ==========================================

def calculate_empennage():
    # 1. Derive main wing span dynamically
    b_w = np.sqrt(AR_w * S_w)
    
    # --- HORIZONTAL STABILIZER SIZING ---
    l_h = X_h - X_aft_cg
    
    # Area based on V_h formula
    S_h = (V_h * S_w * MAC_w) / l_h
    b_h = np.sqrt(AR_h * S_h)
    
    # Chord distributions
    c_root_h = (2 * S_h) / (b_h * (1 + Taper_h))
    c_tip_h = Taper_h * c_root_h
    
    # Sweep angle conversion (Leading Edge to Quarter-Chord)
    tan_sweep_le_h = np.tan(np.radians(Sweep_h_LE))
    tan_sweep_qc_h = tan_sweep_le_h - ((c_root_h - c_tip_h) / b_h)
    sweep_qc_h = np.degrees(np.arctan(tan_sweep_qc_h))

    # --- VERTICAL STABILIZER SIZING ---
    # Estimate tail quarter-chord arm based on fuselage length fraction
    X_v = l_fuselage * X_v_factor
    l_v = X_v - X_aft_cg
    
    # Area based on V_v formula (CRITICAL FIX: Scaled against main wing span b_w)
    S_v = (V_v * S_w * b_w) / l_v
    b_v = np.sqrt(AR_v * S_v)
    
    # Chord distributions
    c_root_v = (2 * S_v) / (b_v * (1 + Taper_v))
    c_tip_v = Taper_v * c_root_v
    
    # Sweep angle conversion
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


# ==========================================
# EXECUTION & OUTPUT
# ==========================================

if __name__ == "__main__":
    res = calculate_empennage()
    
    print("\n" + "="*55)
    print(f"MAIN WING PARAMETERS")
    print("="*55)
    print(f"  Derived Wing Span (b_w):  {res['b_w']:.2f} m")
    
    print("\n" + "-" * 55)
    print("HORIZONTAL STABILIZER:")
    print("-" * 55)
    print(f"  Moment Arm (l_h):         {res['h_tail']['l_h']:.2f} m")
    print(f"  Area (S_h):               {res['h_tail']['S_h']:.2f} m²")
    print(f"  Span (b_h):               {res['h_tail']['b_h']:.2f} m")
    print(f"  Root Chord (c_root):      {res['h_tail']['c_root']:.2f} m")
    print(f"  Tip Chord (c_tip):        {res['h_tail']['c_tip']:.2f} m")
    print(f"  Quarter-Chord Sweep:      {res['h_tail']['sweep_qc']:.1f}°")
    
    print("\n" + "-" * 55)
    print("VERTICAL STABILIZER:")
    print("-" * 55)
    print(f"  Station Arm (X_v):        {res['v_tail']['X_v']:.2f} m")
    print(f"  Moment Arm (l_v):         {res['v_tail']['l_v']:.2f} m")
    print(f"  Area (S_v):               {res['v_tail']['S_v']:.2f} m²")
    print(f"  Span (b_v):               {res['v_tail']['b_v']:.2f} m")
    print(f"  Root Chord (c_root):      {res['v_tail']['c_root']:.2f} m")
    print(f"  Tip Chord (c_tip):        {res['v_tail']['c_tip']:.2f} m")
    print(f"  Quarter-Chord Sweep:      {res['v_tail']['sweep_qc']:.1f}°")
    print("="*55 + "\n")