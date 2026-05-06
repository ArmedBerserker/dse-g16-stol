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
    drag_bucket_cl: Tuple[float, float] # Cl range for low-drag bucket
    stall_type: str          # "Gentle", "Moderate", or "Sharp"

# ==========================================================
# AIRFOIL DATABASE
# ==========================================================

AIRFOIL_DB: List[Airfoil] = [
    
    # --- Original Airfoils ---
    Airfoil(
        name="Wortmann FX 63-137",
        thickness_pct=13.7,
        camber_pct=6.0,
        cl_max_2d=1.70,
        cm_c4=-0.150,
        drag_bucket_cl=(0.5, 1.2),
        stall_type="Gentle"
    ),
    Airfoil(
        name="NACA 2412",
        thickness_pct=12.0,
        camber_pct=2.0,
        cl_max_2d=1.60,
        cm_c4=-0.045,
        drag_bucket_cl=(0.1, 0.5),
        stall_type="Moderate"
    ),
    Airfoil(
        name="NASA LS(1)-0413 (GAW-2)",
        thickness_pct=13.0,
        camber_pct=4.0,
        cl_max_2d=1.80,
        cm_c4=-0.100,
        drag_bucket_cl=(0.2, 0.8),
        stall_type="Moderate"
    ),
    Airfoil(
        name="NACA 4412",
        thickness_pct=12.0,
        camber_pct=4.0,
        cl_max_2d=1.67,
        cm_c4=-0.095,
        drag_bucket_cl=(0.2, 0.7),
        stall_type="Gentle"
    ),
    Airfoil(
        name="NACA 4415",
        thickness_pct=15.0,
        camber_pct=4.0,
        cl_max_2d=1.60,
        cm_c4=-0.095,
        drag_bucket_cl=(0.2, 0.7),
        stall_type="Gentle"
    ),
    Airfoil(
        name="Clark Y",
        thickness_pct=11.7,
        camber_pct=3.4,
        cl_max_2d=1.55,
        cm_c4=-0.080,
        drag_bucket_cl=(0.1, 0.6),
        stall_type="Gentle"
    ),
    Airfoil(
        name="NASA LS(1)-0417 (GAW-1)",
        thickness_pct=17.0,
        camber_pct=4.0,
        cl_max_2d=1.85,
        cm_c4=-0.110,
        drag_bucket_cl=(0.3, 0.9),
        stall_type="Moderate"
    ),
    Airfoil(
        name="USA 35B",
        thickness_pct=11.6,
        camber_pct=4.2,
        cl_max_2d=1.55,
        cm_c4=-0.090,
        drag_bucket_cl=(0.1, 0.7),
        stall_type="Gentle"
    ),
    Airfoil(
        name="Riblett GA30-413.5",
        thickness_pct=13.5,
        camber_pct=4.0,
        cl_max_2d=1.70,
        cm_c4=-0.080,
        drag_bucket_cl=(0.2, 0.8),
        stall_type="Gentle"
    ),
    Airfoil(
        name="NACA 63(3)-418",
        thickness_pct=18.0,
        camber_pct=2.2,
        cl_max_2d=1.55,
        cm_c4=-0.070,
        drag_bucket_cl=(0.2, 0.6),
        stall_type="Moderate"
    ),
    Airfoil(
        name="Selig S1223",
        thickness_pct=12.1,
        camber_pct=8.7,
        cl_max_2d=2.10,          
        cm_c4=-0.200,         
        drag_bucket_cl=(0.8, 1.6),  
        stall_type="Moderate"    
    ),
    Airfoil(
        name="Eppler 423",
        thickness_pct=12.5,
        camber_pct=5.5,
        cl_max_2d=1.90,
        cm_c4=-0.120,
        drag_bucket_cl=(0.5, 1.2),
        stall_type="Gentle"
    ),

   
    Airfoil(
        name="NACA 23012",
        thickness_pct=12.0,
        camber_pct=1.8,
        cl_max_2d=1.60,
        cm_c4=-0.015,                
        drag_bucket_cl=(0.1, 0.4),
        stall_type="Sharp"           
    ),
    Airfoil(
        name="NACA 63A418",
        thickness_pct=18.0,
        camber_pct=2.2,
        cl_max_2d=1.55,
        cm_c4=-0.070,
        drag_bucket_cl=(0.2, 0.6),
        stall_type="Moderate"
    ),
    Airfoil(
        name="NACA 64-514",
        thickness_pct=14.0,
        camber_pct=2.8,
        cl_max_2d=1.60,
        cm_c4=-0.090,
        drag_bucket_cl=(0.3, 0.7),
        stall_type="Moderate"
    ),
    Airfoil(
        name="NACA 63A516",
        thickness_pct=16.0,
        camber_pct=2.8,
        cl_max_2d=1.55,
        cm_c4=-0.090,
        drag_bucket_cl=(0.3, 0.7),
        stall_type="Moderate"
    )
]