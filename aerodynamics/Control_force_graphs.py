import numpy as np
import matplotlib.pyplot as plt

''' TO DO:
    - Change control d'tives and Cmalpha and alpha0'''

A = 9

cmde = -1.28
cm0 = -0.009 * A / (A + 2)  # 0.004
cma = -0.004
cna = 0.9

MTOW = 1840  # kg
W = MTOW * 9.81      # N
rho = 1.225
S = 31.4
a0 = 3         # degrees NOTE: fix

def de_v(cmde, cm0, cma, cna, W, rho, S):
    V = np.linspace(15, 120, 100)  # avoid V=0

    de = -(1 / cmde) * (
        cm0 + (cma / cna) * (W / (0.5 * rho * V**2 * S))
    )

    # Convert delta_e to degrees
    de_deg = np.rad2deg(de)

    plt.figure()
    plt.plot(V, de_deg)
    plt.xlabel('Velocity V [m/s]')
    plt.ylabel(r'$\delta_e$ [deg]')
    plt.grid(True)
    plt.gca().invert_yaxis()
    return de


def de_a(cmde, cm0, cma, a0):
    # Alpha range in degrees
    alpha_deg = np.linspace(-5, 15, 100)

    # Convert alpha and a0 to radians for the calculation
    alpha = np.deg2rad(alpha_deg)
    a0_rad = np.deg2rad(a0)

    de = -(1 / cmde) * (cm0 + cma * (alpha - a0_rad))

    # Convert delta_e to degrees for plotting
    de_deg = np.rad2deg(de)

    plt.figure()
    plt.plot(alpha_deg, de_deg)
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$\delta_e$ [deg]')
    plt.grid(True)
    plt.gca().invert_yaxis()
    return de

def Fe_v(ddeltae_se, Se, ce, W, S, Chd, Cmd, rho, Chdt, dte, dte0, CNa, Cha, Cma):
    V = np.linspace(15, 120, 100)  # avoid V=0

    de_da = 0
    Vh_V2 = 1
    X_n_free_frac = Cmd / CNa * Cha / Chd * (1 - de_da) + Cma / CNa
    Fe = ddeltae_se * Se * ce * Vh_V2 * (W / S * Chd / Cmd * X_n_free_frac - 0.5 * rho * V**2 * Chdt * (dte - dte0))

    plt.figure()
    plt.plot(V, Fe)
    plt.xlabel('Velocity V [m/s]')
    plt.ylabel(r'$F_e$ [N]')
    plt.grid(True)
    plt.gca().invert_yaxis()
    return Fe

def Fe_a(ddeltae_se, Se, ce, Chd, Cmd, rho, Chdt, dte, CLa, Cha, Cma, V):
    # Alpha range in degrees
    alpha_deg = np.linspace(-5, 15, 100)

    # Convert alpha and a0 to radians for the calculation
    alpha = np.deg2rad(alpha_deg)
    a0_rad = np.deg2rad(a0)

    CNa = CLa
    de_da = 0
    X_n_free_frac = Cmd / CNa * Cha / Chd * (1 - de_da) + Cma / CNa
    Cmac = -0.147278639
    ih = 0

    Ch0 = -Chd / Cmd * Cmac - Chd / CNhd * CNhafree * (a0_rad + ih) + Chdt * dte
    Chalpha = -Chd / Cmd * CNa * X_n_free_frac

    Fe = -ddeltae_se * 0.5 * rho * V**2 * Se * ce * (Ch0 + Chalpha * (alpha - a0_rad))

    plt.figure()
    plt.plot(alpha_deg, Fe)
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$F_e$ [N]')
    plt.grid(True)
    plt.gca().invert_yaxis()
    return Fe



de_a = de_a(cm0, cm0, cma, a0)
de_v = de_v(cm0, cm0, cma, cna, W, rho, S)
plt.show()

