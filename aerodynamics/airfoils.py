# possible airfoils
# PARAMETERS TO LOOK AT IN AIRFOIL TRADE-OFF:
# Cd @ Cldes, t/c, trim AoA, Clmax, maximum AoA at stall

from dataclasses import dataclass
from typing import List, Tuple

# ==========================================================
# DATA STRUCTURE
# ==========================================================

@dataclass
class Airfoil:
    name: str
    thickness_pct: float     # Maximum thickness as a % of chord
    camber_pct: float        # Maximum camber as a % of chord
    cl_max_2d: float         # Approximate 2D Max Lift Coefficient (clean)
    cm_c4: float             # Pitching moment at the quarter-chord
    cd_min: float            # Approximate Minimum Drag Coefficient (clean, Re ~3M-5M)
    cl_cd_min: float         # The Cl where absolute minimum drag occurs (Ideal Cl)
    drag_bucket_cl: Tuple[float, float] # The Cl range for the low-drag bucket
    stall_type: str          # "Gentle", "Moderate", or "Sharp"
    best_for: str            # Primary use case

# ==========================================================
# AIRFOIL DATABASE
# ==========================================================

AIRFOIL_DB: List[Airfoil] = [
    
    # --- Recommended Additions for DSE STOL Aircraft ---
    Airfoil(
        name="Wortmann FX 63-137",
        thickness_pct=13.7,
        camber_pct=6.0,
        cl_max_2d=1.70,  
        cm_c4=-0.150,
        cd_min=0.008,
        cl_cd_min=0.85, # Highly cambered, absolute lowest drag is at a high Cl
        drag_bucket_cl=(0.5, 1.2), 
        stall_type="Gentle",
        best_for="Dedicated STOL, extremely high lift at low speeds (motorgliders)"
    ),
    Airfoil(
        name="NACA 2412",
        thickness_pct=12.0,
        camber_pct=2.0,
        cl_max_2d=1.60,
        cm_c4=-0.045,
        cd_min=0.006,
        cl_cd_min=0.25, # Excellent alignment with a fast cruise
        drag_bucket_cl=(0.1, 0.5), 
        stall_type="Moderate",
        best_for="Baseline GA comparison (Cessna 172), lower pitching moment"
    ),
    Airfoil(
        name="NASA LS(1)-0413 (GAW-2)",
        thickness_pct=13.0,
        camber_pct=4.0,
        cl_max_2d=1.80,
        cm_c4=-0.100,
        cd_min=0.007,
        cl_cd_min=0.45, # Sits slightly above cruise Cl, but drag bucket covers cruise well
        drag_bucket_cl=(0.2, 0.8), 
        stall_type="Moderate",
        best_for="Modern STOL, thinner alternative to GAW-1, excellent flap response"
    ),

    # --- Classic High-Lift / STOL Airfoils ---
    Airfoil(
        name="NACA 4412",
        thickness_pct=12.0,
        camber_pct=4.0,
        cl_max_2d=1.67,
        cm_c4=-0.095,
        cd_min=0.0065,
        cl_cd_min=0.40,
        drag_bucket_cl=(0.2, 0.7),
        stall_type="Gentle",
        best_for="Light utility, classic STOL trainers"
    ),
    Airfoil(
        name="NACA 4415",
        thickness_pct=15.0,
        camber_pct=4.0,
        cl_max_2d=1.60,
        cm_c4=-0.095,
        cd_min=0.007,
        cl_cd_min=0.40,
        drag_bucket_cl=(0.2, 0.7),
        stall_type="Gentle",
        best_for="Heavier STOL, strut-braced high wings"
    ),
    Airfoil(
        name="Clark Y",
        thickness_pct=11.7,
        camber_pct=3.4,
        cl_max_2d=1.55,
        cm_c4=-0.080,
        cd_min=0.0075,
        cl_cd_min=0.35, # The classic flat-bottom foil hits your cruise Cl perfectly
        drag_bucket_cl=(0.1, 0.6),
        stall_type="Gentle",
        best_for="Vintage aircraft, bush planes, flat-bottom manufacturing"
    ),
    Airfoil(
        name="NASA LS(1)-0417 (GAW-1)",
        thickness_pct=17.0,
        camber_pct=4.0,
        cl_max_2d=1.85,
        cm_c4=-0.110,
        cd_min=0.008,
        cl_cd_min=0.45,
        drag_bucket_cl=(0.3, 0.9),
        stall_type="Moderate",
        best_for="Modern heavy STOL, advanced GA aircraft"
    ),
    
    # --- Low Drag / General Aviation Airfoils ---
    Airfoil(
        name="NACA 23012",
        thickness_pct=12.0,
        camber_pct=1.8,
        cl_max_2d=1.50,
        cm_c4=-0.015,
        cd_min=0.006,
        cl_cd_min=0.15, # Very low ideal Cl due to forward camber placement
        drag_bucket_cl=(0.0, 0.3), 
        stall_type="Sharp",
        best_for="Fast cruisers, aerobatics, Cessna 208 Caravan"
    ),
    Airfoil(
        name="NACA 23015",
        thickness_pct=15.0,
        camber_pct=1.8,
        cl_max_2d=1.45,
        cm_c4=-0.015,
        cd_min=0.0065,
        cl_cd_min=0.15,
        drag_bucket_cl=(0.0, 0.35),
        stall_type="Sharp",
        best_for="Fast cruisers with higher structural loads"
    )
]

# ==========================================================
# PRINT DATABASE (Example Usage)
# ==========================================================

if __name__ == "__main__":
    print("=========================================================")
    print("                 AIRFOIL CHARACTERISTICS                 ")
    print("=========================================================")
    
    for foil in AIRFOIL_DB:
        print(f"Airfoil: {foil.name}")
        print(f"  - Thickness      : {foil.thickness_pct}%")
        print(f"  - Camber         : {foil.camber_pct}%")
        print(f"  - Max Lift (cl)  : {foil.cl_max_2d:.2f}")
        print(f"  - Pitch Mom (cm) : {foil.cm_c4:.3f}")
        print(f"  - Min Drag (cd)  : {foil.cd_min:.4f} @ Cl = {foil.cl_cd_min:.2f}")
        print(f"  - Drag Bucket Cl : {foil.drag_bucket_cl[0]:.2f} to {foil.drag_bucket_cl[1]:.2f}")
        print(f"  - Stall Quality  : {foil.stall_type}")
        print(f"  - Best App       : {foil.best_for}")
        print("-" * 57)