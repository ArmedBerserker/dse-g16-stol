import sys
import os

# Fix path FIRST, before any local imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import math
import numpy as np
from classes.isa import *
from lookups.consts import *
import matplotlib.pyplot as plt
from pathlib import Path


def stall_speed(mtow, C_l_max, rho, S):
    return math.sqrt((2 * mtow)/(rho * S * C_l_max))


def power_calc(cd0, e, AR, mtow, V, S,alt_m, eta_p, P_shaft,Delta_T):
    #one altitude
    atmos_model = Atmosphere(alt_m, Delta_T)
    rho = atmos_model.density[0]
    C_l = mtow / (0.5 * rho * V ** 2 * S)
    C_d = cd0 + 0.95**2*C_l ** 2 / (np.pi * AR * e)
    P_r = C_d * 0.5 * rho * V ** 3 * S #DV
    P_a = eta_p * P_shaft #eta p will be variable
    excess_power = P_a - P_r
    return excess_power,P_r, P_a

# V = np.linspace(1.3 * stall_speed((1920*9.81), 1.38, 1.25, 25),110, 5000)
# excess_power, P_r, P_a = power_calc(0.03482, 0.752, 10.2, (1920*9.81), V,24.4,2500, 0.84, 202000,0)
# print("power",excess_power)

def RoC_vs_V(excess_power, V, mtow):
    RoC = excess_power / mtow
    idx_max = np.argmax(RoC)

    V_max_RoC = V[idx_max]
    max_RoC = RoC[idx_max]

    print(f"Velocity for maximum RoC: {V_max_RoC}")
    print(f"Maximum RoC: {max_RoC}")
    save_path=Path(__file__).parent / "performance_figures/roc_v.png"
    plt.xlabel("Velocity [m/s]")
    plt.ylabel("Rate of Climb [m/s]")
    plt.grid(True)
    plt.plot(V, RoC)
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    return RoC
    


def AoC_vs_V(P_a, P_r, mtow, V):
    AoC = np.arcsin((P_a - P_r) / (mtow * V))
    AoC = np.rad2deg(AoC)
    idx_max = np.argmax(AoC)

    V_max_AoC = V[idx_max]
    max_AoC = AoC[idx_max]

    #print(f"Velocity for maximum AoC: {V_max_AoC}")
    # print(f"Maximum AoC: {max_AoC}")
    save_path=Path(__file__).parent / "performance_figures/aoc_v.png"
    plt.xlabel("Velocity [m/s]")
    plt.ylabel("Angle of Climb [deg]")
    plt.grid(True)
    plt.plot(V, AoC)
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return AoC


def power_curves_altitude(cd0, e, AR, mtow, S, P_a, alt_m, delta_T):
    
    P_a_alts = []
    P_r_alts = []
    V_alts = []


    for alt in alt_m:

        atmos_model = Atmosphere(alt, delta_T)
        rho = atmos_model.density[0]
        V = np.linspace(1.3*stall_speed(mtow, C_l_max, rho, S), 80, 5000)
        V_alts.append(V)

        C_l = mtow / (0.5 * rho * V**2 * S)
        C_d = cd0 + 0.95**2*C_l**2 / (np.pi * AR * e)

        P_r = C_d * 0.5 * rho * V**3 * S
        P_r_alts.append(P_r)
        P_a = np.ones_like(V) * P_a

        P_a_alts.append(P_a)
        

        line, = plt.plot(V, P_r, label=f"{alt: .0f} m")
        plt.plot(V, P_a, color=line.get_color())
    save_path=Path(__file__).parent / "performance_figures/power_curves_alt.png"
    plt.xlabel("Velocity [m/s]")
    plt.ylabel("Power [W]")
    plt.grid(True)
    plt.legend()
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    return P_a_alts, P_r_alts, V_alts

def climb_gradient(roc, v):
    Climb_grad=roc/V
    return Climb_grad

def RoC_multiple_alts(P_a_alts, P_r_alts, V_alts, mtow,alt_m):
    
    for i, alt in enumerate(alt_m):

        V = V_alts[i]

        roc = (P_a_alts[i] - P_r_alts[i]) / mtow

        plt.plot(V, roc, label=f"{alt} m")
    save_path=Path(__file__).parent / "performance_figures/roc_alt.png"
    plt.xlabel("Velocity [m/s]")
    plt.ylabel("Rate of Climb [m/s]")
    plt.grid(True)
    plt.legend()
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()


def max_RoC_altitude(mtow, P_a, cd0, S, AR, e, delta_T):
    alts = np.arange(0, 10000, 1)  # [m]

    max_RoCs = []

    for a in alts:
        atmos_model = Atmosphere(a, delta_T)
        rho = atmos_model.density[0]
        V = np.linspace(1.3*stall_speed(mtow, C_l_max, rho, S), 80, 5000)

        C_l = mtow / (0.5 * rho * V ** 2 * S)
        C_d = cd0 + 0.95**2*C_l ** 2 / (np.pi * AR * e)

        P_r = C_d * 0.5 * rho * V ** 3 * S
        P_a = np.ones_like(V) * P_a

        RoC = (P_a - P_r) / mtow
        idx_max = np.argmax(RoC)
        max_RoC = RoC[idx_max]

        max_RoCs.append(max_RoC)

    save_path=Path(__file__).parent / "performance_figures/max_roc_alt.png"
    plt.xlabel("Maximum Rate of Climb [m/s]")
    plt.ylabel("Altitude [m]")
    plt.grid(True)
    plt.plot(max_RoCs, alts)
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    return max_RoCs, alts


def envelope(C_l_max, mtow, S, cd0, AR, e, P_a, rho_0, n, delta_T):
    alts = np.arange(0, 10000, 1)  # [m]

    V_s_list = []
    V_max_list = []
    V_max_ROC_list = []
    V_max_AOC_list = []

    for a in alts:
        atmos_model = Atmosphere(a, delta_T)
        rho = atmos_model.density[0]
        V = np.linspace(1.3*stall_speed(mtow, C_l_max, rho, S), 110, 5000)

        V_s = math.sqrt((2 * mtow)/(rho * S * C_l_max))
        V_s_list.append(V_s)

        C_l = mtow / (0.5 * rho * V ** 2 * S)
        C_d = cd0 + C_l ** 2 / (np.pi * AR * e)

        P_r = C_d * 0.5 * rho * V ** 3 * S


        P_a = np.ones_like(V) * P_a

        RoC = (P_a - P_r) / mtow
        idx_max = np.argmax(RoC)
        Vmax_RoC = V[idx_max]
        V_max_ROC_list.append(Vmax_RoC)

        AoC = np.arcsin((P_a - P_r) / (mtow * V))
        AoC = np.rad2deg(AoC)
        idx_max = np.argmax(AoC)
        V_max_AoC = V[idx_max]
        V_max_AOC_list.append(V_max_AoC)

        diff = P_a - P_r
        crossing_indices = np.where(np.diff(np.sign(diff)) != 0)[0]
        intersections = []
        
        for i in crossing_indices:
            # Linear interpolation
            x1, x2_ = V[i], V[i + 1]
            d1, d2 = diff[i], diff[i + 1]
        
            # intersection x
            xi = x1 - d1 * (x2_ - x1) / (d2 - d1)
        
            # intersection y
            yi = P_a[i] + (P_a[i + 1] - P_a[i]) * (xi - x1) / (x2_ - x1)
        
            intersections.append((xi, yi))
        
        #print("Intersection points:")
        # for pt in intersections:
        #     print(pt)
        
        if len(intersections) == 1:
            i = crossing_indices[0]
        
            # Check sign change direction
            before = diff[i]
            after = diff[i + 1]
        
            if before < 0 and after > 0:
                kind = "minimum"
            elif before > 0 and after < 0:
                kind = "maximum"
            else:
                kind = "undetermined"
        
            #print(f"\nSingle intersection corresponds to a local {kind}.")
        idx = np.where(np.diff(np.sign(P_a - P_r)))[0]
        #print(idx)
        
        V_max = intersections[-1][0]
        V_max_list.append((V_max))

    # print(V_max_list)
    #print(V_s_list)\
    save_path=Path(__file__).parent / "performance_figures/envelope.png"
    plt.xlabel("Velocity [m/s]")
    plt.ylabel("Altitude [m]")
    plt.grid(True)
    plt.plot(V_s_list, alts)
    plt.plot(V_max_AOC_list, alts)
    plt.plot(V_max_ROC_list, alts)
    plt.plot(V_max_list,alts)
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    



if __name__ == '__main__':
    alt = 0  # [m]
    delta_T = 20  # [K]

    
    # ---- Design parameters ----
    e = 0.752
    AR = 10.2
    mtow = 1912 * g  # [N]

    atmos_model = Atmosphere(alt, delta_T)
    rho_0 = atmos_model.density[0]

    C_l_max = 1.38 #cruise
    S = 24.3
    n = 1
    cd0 = 0.0342

    V = np.linspace(1.3 * stall_speed(mtow, C_l_max, rho_0, S),90, 5000)
    eta_p = 0.84
    cruise_alt = 2500 #m cruise alt
    alt_m = [0,500,1000,1500,2000,2500,3000]
    P_shaft = 202000
    excess_power, P_r, P_a = power_calc(cd0, e, AR,mtow, V, S,cruise_alt,eta_p,P_shaft,0)
    roc = RoC_vs_V(excess_power, V, mtow) #cruise
    Aoc = AoC_vs_V(P_a, P_r, mtow, V) #cruise
    P_a_alts, P_r_alts, V_alts = power_curves_altitude(cd0, e, AR, mtow, S, P_a,alt_m, 0)
    RoC_multiple_alts(P_a_alts, P_r_alts, V_alts, mtow,alt_m)
    max_RoC_altitude(mtow, P_a,cd0, S, AR, e, 0)
    envelope(C_l_max, mtow, S, cd0, AR, e, P_a, rho_0, n, delta_T)
    


    


    


    




#     #plt.figure(figsize=(12, 6))

#     # ---- Loop through all combinations ----
#     for design_name, data in designs.items():

#         cd0 = data["cd0"]
#         P_a_0 = data["P_a"]
#         n = data["n"]

#         # Power available
#         P_a = P_a_0 * (rho_0 / 1.225) ** n #changes
#         P_a_list = P_a * np.ones_like(V)

#         # Drag coefficient
#         C_d = calculate_cd(cd0, e, AR, mtow, rho_0, V, S)

#         # Power required
#         P_r_list = calculate_power_required(C_d, rho_0, V, S)
        
#         excess_power = power_curves(P_a_list, P_r_list, V)
#         # Rate of climb
#         # RoC = RoC_vs_V(excess_power, V, mtow)

#         Aoc = AoC_vs_V(P_a_list, P_r_list, mtow,V)
        
#         # P_a_lists, P_r_lists, V_alts = power_curves_altitude(cd0, e, AR, mtow, rho_0, S, P_a, n, delta_T)

#         #max_RoCs, alts = RoC_altitude(mtow, P_a_0, rho_0, n, cd0, S, AR, e, delta_T)

#         # RoC_multiple_alts(P_a_lists, P_r_lists, V_alts, mtow)

#         #envelope(C_l_max, mtow, S, cd0, AR, e, P_a_0, rho_0, n, delta_T)

        
#     #     plt.plot(V, RoC, label=design_name)
#     # plt.xlabel("Velocity [m/s]")
#     # plt.ylabel("Rate of Climb [m/s]")
#     # plt.grid(True)
#     # plt.legend()
#     # plt.show()
#         plt.plot(V, Aoc, label=design_name)
#     plt.xlabel("Velocity [m/s]")
#     plt.ylabel("Angle of Climb [deg]")
#     plt.grid(True)
#     plt.legend()
#     plt.show()
#     #     plt.plot(max_RoCs, alts, label=design_name)
#     # plt.xlabel("Maximum Rate of Climb [m/s]")
#     # plt.ylabel("Altitude [m]")
#     # plt.grid(True)
#     # plt.legend()
#     # plt.show()


    