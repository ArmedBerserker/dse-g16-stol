import sys
import os
from scipy.integrate import quad
# Fix path FIRST, before any local imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import math
import numpy as np
from classes.isa import *
from lookups.consts import *
import matplotlib.pyplot as plt
from pathlib import Path
from Propeller_performance import *
from scipy.optimize import brentq
from scipy.interpolate import interp1d

def stall_speed(mtow, C_l_max, rho, S):
    return math.sqrt((2 * mtow)/(rho * S * C_l_max))
print(stall_speed(1821*9.81, 1.38, 1.225, 24.23))

def power_calc(cd0, e, AR, mtow, V, S,alt_m, P_shaft,Delta_T, D_ft):
    #one altitude
    
    atmos_model = Atmosphere(alt_m, Delta_T)
    rho = atmos_model.density[0]
    V_ms,P_useful_w= curve(D_ft, rho, P_shaft)
    C_l = mtow / (0.5 * rho * V ** 2 * S)
    C_d = cd0 + 0.95**2*C_l ** 2 / (np.pi * AR * e)
    P_r = C_d * 0.5 * rho * V ** 3 * S #DV
    P_a = np.interp(V, V_ms, P_useful_w)
    excess_power = P_a - P_r
    # print("Pa", P_a)
    # print("P_r", P_r)
    # print("excess", excess_power)
    # plt.plot(V, P_a)
    return excess_power,P_r, P_a

# V = np.linspace(1.3 * stall_speed((1920*9.81), 1.38, 1.25, 25),110, 5000)
# excess_power, P_r, P_a = power_calc(0.03482, 0.752, 10.2, (1920*9.81), V,24.4,2500, 0.84, 202000,0)
# print("power",excess_power)

def RoC_vs_V(excess_power, V, mtow):
    RoC = excess_power / mtow
    idx_max = np.argmax(RoC)

    V_max_RoC = V[idx_max]
    max_RoC = RoC[idx_max]

    # print(f"Velocity for maximum RoC: {V_max_RoC}")
    # print(f"Maximum RoC: {max_RoC}")
    save_path=Path(__file__).parent / "performance_figures/roc_v.png"
    plt.xlabel("Velocity [m/s]")
    plt.ylabel("Rate of Climb [m/s]")
    plt.ylim(-4,6)
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
    plt.ylim(-2,9)
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(max_AoC)
    
    return AoC


def power_curves_altitude(cd0, e, AR, mtow, S, alt_m, delta_T,P_shaft,D_ft):
    
    P_a_alts = []
    P_r_alts = []
    V_alts = []


    for alt in alt_m:

        atmos_model = Atmosphere(alt, delta_T)
        rho = atmos_model.density[0]
        V = np.linspace(1.3*stall_speed(mtow, C_l_max, rho, S), 90, 1000)
        V_alts.append(V)

        C_l = mtow / (0.5 * rho * V**2 * S)
        C_d = cd0 + 0.95**2*C_l**2 / (np.pi * AR * e)

        P_r = C_d * 0.5 * rho * V**3 * S
        P_r_alts.append(P_r)
        # print(rho)
        V_ms,P_a= curve(D_ft, rho, P_shaft)
        # print(f"alt={alt}, max P_a={np.max(P_a):.1f} W")

        P_a_alts.append((V_ms, P_a))
        # P_a_interp = np.interp(V, V_ms, P_a)

        line, = plt.plot(V, P_r, label=f"{alt: .0f} m")
        #plt.ylim(150000,200000)
        plt.plot(V_ms, P_a, linestyle='--',color=line.get_color())
        plt.xlim(30,90)
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
        V = V_alts[i]  # P_r grid
        V_ms, P_a = P_a_alts[i]  # unpack both (store them as tuples when you compute them)
        
        P_a_interp = np.interp(V, V_ms, P_a)  # interpolate P_a onto P_r's V grid
        
        roc = (P_a_interp - P_r_alts[i]) / mtow
        plt.plot(V, roc, label=f"{alt} m")
    save_path=Path(__file__).parent / "performance_figures/roc_alt.png"
    plt.xlabel("Velocity [m/s]")
    plt.ylabel("Rate of Climb [m/s]")
    plt.ylim(-2,6)
    plt.grid(True)
    plt.legend()
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

def AoC_multiple_alts(P_a_alts, P_r_alts, V_alts, mtow,alt_m):
    
    for i, alt in enumerate(alt_m):
        V = V_alts[i]  # P_r grid
        V_ms, P_a = P_a_alts[i]  # unpack both (store them as tuples when you compute them)
        
        P_a_interp = np.interp(V, V_ms, P_a)  # interpolate P_a onto P_r's V grid
        
        aoc = np.arcsin((P_a_interp[i] - P_r_alts[i]) / (mtow*V))
        aoc = np.rad2deg(aoc)
        plt.plot(V, aoc, label=f"{alt} m")
    save_path=Path(__file__).parent / "performance_figures/aoc_alt.png"
    plt.xlabel("Velocity [m/s]")
    plt.ylabel("Angle of Climb [deg]")
    plt.ylim(-2,10)
    plt.grid(True)
    plt.legend()
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

def max_RoC_altitude(mtow, cd0, S, AR, e, delta_T,P_shaft,D_ft):
    alts = np.arange(0, 2590.5, 0.5)  # [m]

    max_RoCs = []

    for a in alts:
        atmos_model = Atmosphere(a, delta_T)
        rho = atmos_model.density[0]
        #V = np.linspace(1.3*stall_speed(mtow, C_l_max, rho, S), 80, 1000)
        # V = np.linspace(30, 110, 10000)
        V_ms,P_a= curve(D_ft, rho, P_shaft)
        V = V_ms
        C_l = mtow / (0.5 * rho * V ** 2 * S)
        C_d = cd0 + 0.95**2*C_l ** 2 / (np.pi * AR * e)

        P_r = C_d * 0.5 * rho * V ** 3 * S
        # T_engine, P_a, V_ms = curve(D_ft, rho, P_shaft)
        # P_a_interp = np.interp(V, V_ms, P_a)
        
        

        RoC = (P_a - P_r) / mtow
        idx_max = np.argmax(RoC)
        max_RoC = RoC[idx_max]

        max_RoCs.append(max_RoC)

    # import numpy as np

    integrand = [1.0 / roc for roc in max_RoCs]

    time_to_climb = np.trapz(integrand, alts)
    print(time_to_climb)

    #print(f"Time to climb = {time_to_climb:.1f} s")
    #print(f"= {time_to_climb/60:.2f} minutes")
    save_path=Path(__file__).parent / "performance_figures/max_roc_alt.png"
    plt.xlabel("Maximum Rate of Climb [m/s]")
    plt.ylabel("Altitude [m]")
    # plt.ylim(-4,6)
    plt.grid(True)
    plt.plot(max_RoCs, alts)
    plt.xlim(1, 8)
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()


    return max_RoCs, alts


def envelope(C_l_max, mtow, S, cd0, AR, e, P_a, rho_0, n, delta_T, P_shaft, D_ft):
    alts = np.arange(0, 3000, 1)  # [m]

    V_s_list = []
    V_max_list = []
    V_max_ROC_list = []
    V_max_AOC_list = []
    V_min_list = []

    for a in alts:
        atmos_model = Atmosphere(a, delta_T)
        rho = atmos_model.density[0]

        V_s = math.sqrt((2 * mtow)/(rho * S * C_l_max))
        V_s_list.append(V_s)
        V = np.linspace(0.8*V_s, 110, 1000)
        V_ms,P_a = curve(D_ft, rho, P_shaft)

        C_l = mtow / (0.5 * rho * V ** 2 * S)
        C_d = cd0 + C_l ** 2 / (np.pi * AR * e)

        P_r = C_d * 0.5 * rho * V ** 3 * S


        P_a = np.interp(V, V_ms, P_a)

        RoC = (P_a - P_r) / mtow
        idx_max = np.argmax(RoC)
        Vmax_RoC = V[idx_max]
        V_max_ROC_list.append(Vmax_RoC)

        V_aoc = np.linspace(1.3 * V_s, 150, 10000)
        P_r_aoc = (cd0 + (mtow / (0.5 * rho * V_aoc ** 2 * S)) ** 2 / (np.pi * AR * e)) * 0.5 * rho * V_aoc ** 3 * S
        P_a_aoc = np.interp(V_aoc, V_ms, P_a)  # interpolate on V_aoc, not V

        arg = np.clip((P_a_aoc - P_r_aoc) / (mtow * V_aoc), -1, 1)  # clip to avoid arcsin domain errors
        AoC = np.rad2deg(np.arcsin(arg))
        idx_max = np.argmax(AoC)
        V_max_AoC = V_aoc[idx_max]
        V_max_AOC_list.append(V_max_AoC)

        diff = P_a - P_r
        crossing_indices = np.where(np.diff(np.sign(diff)) != 0)[0]
        intersections = []
        
        for i in crossing_indices:
            x1, x2_ = V[i], V[i + 1]
            d1, d2 = diff[i], diff[i + 1]
            if d2 == d1:
                continue  # avoid division by zero
            xi = x1 - d1 * (x2_ - x1) / (d2 - d1)
            intersections.append(xi)
        
        # Ceiling detection: no intersections means P_a < P_r everywhere
        if len(intersections) == 0:
            V_min_list.append(np.nan)
            V_max_list.append(np.nan)

        elif len(intersections) == 1:
            i = crossing_indices[0]
            if diff[i] > 0 and diff[i + 1] < 0:
                # Only a Vmax crossing exists → aircraft is stall-limited at low speed
                V_min_list.append(np.nan)      # no power-limited Vmin
                V_max_list.append(intersections[0])
            else:
                # Only a Vmin crossing exists → no Vmax found (above ceiling?)
                V_min_list.append(intersections[0])
                V_max_list.append(np.nan)
        else:
            i0 = crossing_indices[0]
            if diff[i0] > 0 and diff[i0 + 1] < 0:
                # First crossing is already a downward one → stall-limited, no power Vmin
                V_min_list.append(np.nan)
                V_max_list.append(intersections[-1])
            else:
                # Genuine two-sided power envelope
                V_min_list.append(intersections[0])
                V_max_list.append(intersections[-1])

        
    # print(V_max_list)
    #print(V_s_list)\
    save_path=Path(__file__).parent / "performance_figures/envelope.png"
    plt.xlabel("Velocity [m/s]")
    plt.ylabel("Altitude [m]")
    plt.grid(True)
    plt.plot(V_s_list, alts, label = "stall speed")
    plt.plot(V_max_AOC_list, alts, label = "max aoc")
    plt.plot(V_max_ROC_list, alts, label = "max roc")
    plt.plot(V_max_list,alts, label = "Vmax")
    plt.plot(V_min_list, alts, label = "Vmin")
    plt.legend()
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    






if __name__ == '__main__':
    alt = 0  # [m]
    delta_T = 0  # [K]

    
    # ---- Design parameters ----
    e = 0.747
    AR = 9
    mtow = 1871 * 9.81  # [N]

    atmos_model = Atmosphere(alt, delta_T)
    rho_0 = atmos_model.density[0]
    
    C_l_max = 1.44+0.287 #cruise
    S = 31.4
    n = 1
    # print(rho_0)
    # print(mtow)
    # print(S)
    # print(C_l_max)
    cd0 = 0.0314
    D_ft=5.675853 #propeller diameter
    V = np.linspace(1.3 * stall_speed(mtow, C_l_max, rho_0, S), 90, 1000)
    
    # print(V)
    cruise_alt = 2590.8 #m cruise alt
    alt_m = [0,500,1000,1500,2000,2500,3000]
    P_shaft = 135.443 #hp

    excess_power, P_r, P_a = power_calc(cd0, e, AR,mtow, V, S,0,P_shaft,0,D_ft)

    roc = RoC_vs_V(excess_power, V, mtow) #cruise
    Aoc = AoC_vs_V(P_a, P_r, mtow, V) #cruise
    P_a_alts, P_r_alts, V_alts = power_curves_altitude(cd0, e, AR, mtow, S,alt_m, 0,P_shaft, D_ft)
    RoC_multiple_alts(P_a_alts, P_r_alts, V_alts, mtow,alt_m)
    AoC_multiple_alts(P_a_alts, P_r_alts, V_alts, mtow,alt_m)
    max_roc, alt = max_RoC_altitude(mtow,cd0, S, AR, e, 0,P_shaft, D_ft)
    envelope(C_l_max, mtow, S, cd0, AR, e, P_a, rho_0, n, delta_T,P_shaft,D_ft)
    print(max_roc)

    


    


    


    




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


    