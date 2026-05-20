import math

import numpy as np
from classes.isa import *
from lookups.consts import *
import matplotlib.pyplot as plt


def calculate_cd(cd0, e, AR, mtom, rho, V, S):
    C_l = mtom / (0.5 * rho * V ** 2 * S)
    C_d = cd0 + C_l ** 2 / (np.pi * AR * e)

    return C_d


def calculate_power_required(C_d, rho, V, S):

    P_r = C_d * 0.5 * rho * V ** 3 * S

    return P_r


def power_curves(P_a, P_r, V):

    P_a = P_a * np.ones_like(V)
    # plt.plot(V, P_a)
    # plt.plot(V, P_r)
    # plt.show()

    return P_a - P_r


def RoC_vs_V(excess_power, V, mtom):
    RoC = excess_power / mtom
    idx_max = np.argmax(RoC)

    V_max_RoC = V[idx_max]
    max_RoC = RoC[idx_max]

    # print(f"Velocity for maximum RoC: {V_max_RoC}")
    # print(f"Maximum RoC: {max_RoC}")
    # plt.plot(V, RoC)
    # plt.show()


def AoC_vs_V(P_a, P_r, mtom, V):
    AoC = np.arcsin((P_a - P_r) / (mtom * V))
    AoC = np.rad2deg(AoC)
    idx_max = np.argmax(AoC)

    V_max_AoC = V[idx_max]
    max_AoC = AoC[idx_max]

    # print(f"Velocity for maximum AoC: {V_max_AoC}")
    # print(f"Maximum RoC: {max_AoC}")
    # plt.plot(V, AoC)
    # plt.show()


def power_curves_altitude(cd0, e, AR, W, rho_0, V, S, P_a_0, n):

    altitudes_ft = [0, 2000, 4000, 6000, 8000, 10000]

    P_a_alts = []
    P_r_alts = []

    for alt_ft in altitudes_ft:

        alt_m = alt_ft * FT_TO_M

        atmos_model = Atmosphere(alt_m, 0)
        rho = atmos_model.density[0]

        # --- Power Required ---

        C_l = W / (0.5 * rho * V**2 * S)

        C_d = cd0 + C_l**2 / (np.pi * AR * e)

        P_r = C_d * 0.5 * rho * V**3 * S

        P_r_alts.append(P_r)

        # --- Power Available ---

        P_a = P_a_0 * (rho / rho_0) ** n

        # Make horizontal line
        P_a_curve = np.ones_like(V) * P_a

        P_a_alts.append(P_a_curve)

        # Plot power required first
        # line, = plt.plot(V, P_r, label=f"{alt_ft} ft")

    #     # Use same color for power available
    #     plt.plot(
    #         V,
    #         P_a_curve,
    #         '--',
    #         color=line.get_color()
    #     )
    #
    #
    # plt.xlabel("Velocity [m/s]")
    # plt.ylabel("Power [W]")
    # plt.title("Power Required and Available vs Velocity")
    # plt.grid(True)
    # plt.legend()
    # plt.show()

    return P_a_alts, P_r_alts


def RoC_multiple_alts(P_a, P_r, V, mtom):

    RoC = [(P_a[i] - P_r[i]) / mtom for i in range(len(P_a))]

    for r in RoC:
        plt.plot(V[0:50], r[0:50])

    plt.show()


def RoC_altitude(W, P_a_0, rho_0, n, V, cd0, S, AR, e):
    alts = np.arange(0, 10000, 1)
    max_RoCs = []

    for a in alts:
        atmos_model = Atmosphere(a, 0)
        rho = atmos_model.density[0]

        C_l = W / (0.5 * rho * V ** 2 * S)

        C_d = cd0 + C_l ** 2 / (np.pi * AR * e)

        P_r = C_d * 0.5 * rho * V ** 3 * S

        P_a = P_a_0 * (rho / rho_0) ** n
        P_a = np.ones_like(V) * P_a

        RoC = (P_a - P_r) / W
        idx_max = np.argmax(RoC)
        max_RoC = RoC[idx_max]

        max_RoCs.append(max_RoC)

    # print(max_RoCs[-760])
    # print(alts[-760])
    plt.plot(max_RoCs, alts)
    plt.show()



if __name__ == '__main__':
    alt = 0
    cd0 = 0.015
    e = 0.8
    AR = 9
    mtom = 2000 * g
    atmos_model = Atmosphere(alt, 0)
    rho = atmos_model.density[0]
    V = np.arange(25, 100, 1)
    S = 30 #25.675
    eta_p = 0.8
    P_a = 200000 * eta_p

    C_d = calculate_cd(cd0, e, AR, mtom, rho, V, S)

    P_r = calculate_power_required(C_d, rho, V, S)

    excess_power = power_curves(P_a, P_r, V)

    # RoC_vs_V(excess_power, V, mtom)
    #
    # AoC_vs_V(P_a, P_r, mtom, V)
    n = 0.8 # 1 for piston / 0.8 for turboprop

    P_a, P_r = power_curves_altitude(cd0, e, AR, mtom, rho, V, S, P_a, n)

    # RoC_multiple_alts(P_a, P_r, V, mtom)

    P_a = 200000 * eta_p
    RoC_altitude(mtom, P_a, rho, n, V, cd0, S, AR, e)

# length, width, height, length of nose cone and tail cone, area of base, X location of seats, X location of cargo, volume of cabin
