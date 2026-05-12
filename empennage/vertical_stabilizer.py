import math

# --- Constants / Design Choices (TODO: link these to yaml files) ---
X_v_factor = 0.9
V_v = 0.06275      # Roskam volume coefficient
AR_v = 1.5         # ADSEE average aspect ratio
Sweep_v = 25       # degrees
lambda_v = 0.4     # ADSEE taper ratio
T_over_C_v = 0.15  # Roskam thickness-to-chord

def calculate_vertical_stabilizer_geometry(l_f, MAC, S, X_aft_cg=1.0):
   
    # --- Calculations ---
    X_v = X_v_factor * l_f
    
    # Vertical tail area based on volume coefficient formula
    # S_v = (V_v * MAC * S) / (X_v - X_aft_cg)
    S_v = (0.36 * MAC) / (X_v - X_aft_cg) * S 
    
    b_v = (AR_v * S_v)**0.5
    Cr_v = (2 * S_v) / (b_v * (1 + lambda_v))
    Ct_v = lambda_v * Cr_v
    
    return {
        "S_v": S_v,
        "b_v": b_v,
        "Cr_v": Cr_v,
        "Ct_v": Ct_v,
        "X_v": X_v,
        "AR_v": AR_v,
        "Sweep_v": Sweep_v
    }

