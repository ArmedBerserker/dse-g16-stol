import sys
import os

# Fix path FIRST, before any local imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import math
import numpy as np
from classes.isa import *
from lookups.consts import *
import matplotlib.pyplot as plt
from scipy.optimize import brentq


def stall_speed(mtow, C_l_max, rho, S):
    return math.sqrt((2 * mtow)/(rho * S * C_l_max))


# def calculate_cd(cd0, e, AR, mtow, rho_0, V, S):
#     C_l = mtow / (0.5 * rho_0 * V ** 2 * S)
#     C_d = cd0 + C_l ** 2 / (np.pi * AR * e)

#     return C_d


# def calculate_power_required(C_d, rho_0, V, S):
#     P_r = C_d * 0.5 * rho_0 * V ** 3 * S
#     return P_r


def power_curves(P_a_list, P_r_list, V):
    #delta_t 20 temperature not taken into account
    # plt.xlabel("Velocity [m/s]")
    # plt.ylabel("Power [W]")
    # plt.grid(True)
    # plt.plot(V, P_a_list)
    # plt.plot(V, P_r_list)
    # plt.show()

    return P_a_list - P_r_list

def power_curves(cd0, e, AR, mtow, V, S,alt_m, eta_p, P_shaft):

    atmos_model = Atmosphere(alt_m, 0)
    rho = atmos_model.density[0]
    C_l = mtow / (0.5 * rho * V ** 2 * S)
    C_d = cd0 + C_l ** 2 / (np.pi * AR * e)
    P_r = C_d * 0.5 * rho * V ** 3 * S #DV
    P_a = 231000 #eta_p * P_shaft #eta p will be variable
    excess_power = P_a - P_r
    return excess_power

excess_power = power_curves(0.03487, 0.752, 10.2, 1865, 68,23.7,3048, 0.84, 298000)
print("power",excess_power)




def RoC_vs_V(excess_power, V, mtow):
    RoC = excess_power / mtow
    idx_max = np.argmax(RoC)

    V_max_RoC = V[idx_max]
    max_RoC = RoC[idx_max]

    print(f"Velocity for maximum RoC: {V_max_RoC}")
    print(f"Maximum RoC: {max_RoC}")
    return RoC
    # plt.xlabel("Velocity [m/s]")
    # plt.ylabel("Rate of Climb [m/s]")
    # plt.grid(True)
    # plt.plot(V, RoC)
    # plt.show()


def AoC_vs_V(P_a_list, P_r_list, mtow, V):
    AoC = np.arcsin((P_a_list - P_r_list) / (mtow * V))
    AoC = np.rad2deg(AoC)
    idx_max = np.argmax(AoC)

    V_max_AoC = V[idx_max]
    max_AoC = AoC[idx_max]

    #print(f"Velocity for maximum AoC: {V_max_AoC}")
    # print(f"Maximum AoC: {max_AoC}")
    # plt.xlabel("Velocity [m/s]")
    # plt.ylabel("Angle of Climb [deg]")
    # plt.grid(True)
    # plt.plot(V, AoC)
    # plt.show()
    return AoC


def power_curves_altitude(cd0, e, AR, mtow, rho_0, S, P_a_0, n, delta_T):
    altitudes_ft = [0, 2000, 4000, 6000, 8000, 10000]

    P_a_alts = []
    P_r_alts = []
    V_alts = []

    for alt_ft in altitudes_ft:

        alt_m = alt_ft * FT_TO_M

        atmos_model = Atmosphere(alt_m, delta_T)
        rho = atmos_model.density[0]
        #V = np.arange(int(1.3*stall_speed(mtow, C_l_max, rho, S)), 110, 1)
        V = np.linspace(1.3*stall_speed(mtow, C_l_max, rho, S), 110, 5000)
        V_alts.append(V)

        C_l = mtow / (0.5 * rho * V**2 * S)
        C_d = cd0 + C_l**2 / (np.pi * AR * e)

        P_r = C_d * 0.5 * rho * V**3 * S
        P_r_alts.append(P_r)
        P_a = P_a_0 * (rho / rho_0) ** n
        P_a = np.ones_like(V) * P_a

        P_a_alts.append(P_a)

        line, = plt.plot(V, P_r, label=f"{alt_m: .0f} m")
        plt.plot(V, P_a, '--', color=line.get_color())

    plt.xlabel("Velocity [m/s]")
    plt.ylabel("Power [W]")
    plt.grid(True)
    plt.legend()
    plt.show()

    return P_a_alts, P_r_alts, V_alts

def climb_gradient(roc, v):
    Climb_grad=roc/V
    return Climb_grad

def RoC_multiple_alts(P_a_lists, P_r_lists, V_alts, mtow):
    
    altitudes_ft = [0, 2000, 4000, 6000, 8000, 10000]

    for i, alt in enumerate(altitudes_ft):

        V = V_alts[i]

        roc = (P_a_lists[i] - P_r_lists[i]) / mtow

        plt.plot(V, roc, label=f"{alt} ft")

    plt.xlabel("Velocity [m/s]")
    plt.ylabel("Rate of Climb [m/s]")
    plt.grid(True)
    plt.legend()
    plt.show()


def RoC_altitude(mtow, P_a_0, rho_0, n, cd0, S, AR, e, delta_T):
    alts = np.arange(0, 10000, 1)  # [ft]
    alts = alts * FT_TO_M

    max_RoCs = []

    for a in alts:
        atmos_model = Atmosphere(a, delta_T)
        rho = atmos_model.density[0]
        #V = np.arange(int(1.3*stall_speed(mtow, C_l_max, rho, S)), 110, 1)
        V = np.linspace(1.3*stall_speed(mtow, C_l_max, rho, S), 110, 5000)

        C_l = mtow / (0.5 * rho * V ** 2 * S)
        C_d = cd0 + C_l ** 2 / (np.pi * AR * e)

        P_r = C_d * 0.5 * rho * V ** 3 * S
        P_a = P_a_0 * (rho / rho_0) ** n
        P_a = np.ones_like(V) * P_a

        RoC = (P_a - P_r) / mtow
        idx_max = np.argmax(RoC)
        max_RoC = RoC[idx_max]

        max_RoCs.append(max_RoC)

    # print(max_RoCs[-2000])
    # print(alts[-2000])
    # plt.xlabel("Maximum Rate of Climb [m/s]")
    # plt.ylabel("Altitude [m]")
    # plt.grid(True)
    # plt.plot(max_RoCs, alts)
    # plt.show()
    return max_RoCs, alts


def envelope(C_l_max, mtow, S, cd0, AR, e, P_a_0, rho_0, n, delta_T):
    alts = np.arange(0, 10000, 1)  # [ft]
    alts = alts * FT_TO_M

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

        P_a = P_a_0 * (rho / rho_0) ** n
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
    #print(V_s_list)
    plt.xlabel("Velocity [m/s]")
    plt.ylabel("Altitude [m]")
    plt.grid(True)
    plt.plot(V_s_list, alts)
    plt.plot(V_max_AOC_list, alts)
    plt.plot(V_max_ROC_list, alts)
    plt.plot(V_max_list,alts)
    plt.show()



if __name__ == '__main__':
    alt = 0  # [m]
    delta_T = 20  # [K]

    
    # ---- Design parameters ----
    e = 0.8
    AR = 9
    mtow = 2000 * g  # [N]

    atmos_model = Atmosphere(alt, delta_T)
    rho_0 = atmos_model.density[0]

    C_l_max = 1.7 #cruise
    S = 25.6

    V = np.linspace(1.3 * stall_speed(mtow, C_l_max, rho_0, S),110, 5000)
    eta_p = 0.8

    designs = {

    # "Taildragger - Piston": {
    #     "cd0": 0.0255,
    #     "P_a": 202000 * eta_p,
    #     "n": 1.0
    # },

    # "Taildragger - Turboprop": {
    #     "cd0": 0.0255,
    #     "P_a": 550000 * eta_p,
    #     "n": 0.8
    # },

    "Boosted Piston": {
        "cd0": 0.04,
        "P_a": 202000 * eta_p,
        "P_a": 202000 * eta_p, # take-off power: 234000
        "n": 1.0
    },

    "Boosted Turboprop": {
        "cd0": 0.04,
        "P_a": 550000 * eta_p,
        "P_a": 550000 * eta_p,  # take-off power: 626000
        "n": 0.8
    }
}

    #plt.figure(figsize=(12, 6))

    # ---- Loop through all combinations ----
    for design_name, data in designs.items():

        cd0 = data["cd0"]
        P_a_0 = data["P_a"]
        n = data["n"]

        # Power available
        P_a = P_a_0 * (rho_0 / 1.225) ** n #changes
        P_a_list = P_a * np.ones_like(V)

        # Drag coefficient
        C_d = calculate_cd(cd0, e, AR, mtow, rho_0, V, S)

        # Power required
        P_r_list = calculate_power_required(C_d, rho_0, V, S)
        
        excess_power = power_curves(P_a_list, P_r_list, V)
        # Rate of climb
        # RoC = RoC_vs_V(excess_power, V, mtow)

        Aoc = AoC_vs_V(P_a_list, P_r_list, mtow,V)
        
        # P_a_lists, P_r_lists, V_alts = power_curves_altitude(cd0, e, AR, mtow, rho_0, S, P_a, n, delta_T)

        #max_RoCs, alts = RoC_altitude(mtow, P_a_0, rho_0, n, cd0, S, AR, e, delta_T)

        # RoC_multiple_alts(P_a_lists, P_r_lists, V_alts, mtow)

        #envelope(C_l_max, mtow, S, cd0, AR, e, P_a_0, rho_0, n, delta_T)

        
    #     plt.plot(V, RoC, label=design_name)
    # plt.xlabel("Velocity [m/s]")
    # plt.ylabel("Rate of Climb [m/s]")
    # plt.grid(True)
    # plt.legend()
    # plt.show()
        plt.plot(V, Aoc, label=design_name)
    plt.xlabel("Velocity [m/s]")
    plt.ylabel("Angle of Climb [deg]")
    plt.grid(True)
    plt.legend()
    plt.show()
    #     plt.plot(max_RoCs, alts, label=design_name)
    # plt.xlabel("Maximum Rate of Climb [m/s]")
    # plt.ylabel("Altitude [m]")
    # plt.grid(True)
    # plt.legend()
    # plt.show()


    