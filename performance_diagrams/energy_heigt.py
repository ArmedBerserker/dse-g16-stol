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
from performance_diagrams import *


alts = np.arange(0,10001,50)
print((alts))
#V_array = np.arange(1,101,1)
#print(np.shape(V_array))

def stall_speed(mtow, C_l_max, rho, S):
    return math.sqrt((2 * mtow)/(rho * S * C_l_max))

def h_v_plot(alts, mtow,cd0,S,AR,e,P_a,c_l_max):
    ROC = np.zeros((len(alts), 100))
    He = np.zeros((len(alts), 100))
    for i, alt in enumerate(alts):
        rho = Atmosphere(alt, 0).density[0]
        V_array=np.linspace(stall_speed(mtow,c_l_max, rho, S), 101, 100)
        for j, V in enumerate(V_array):
            Cl = mtow / (0.5 * rho * V**2 * S)
            Cd = cd0 + Cl**2/(np.pi*AR*e)
            P_r = 0.5 * rho * V**3 * S * Cd
            P_a = P_a
            ROC[i,j] = (P_a-P_r)/mtow
            He[i, j] = alt + V**2 / (2 * 9.81)
    
   
    print(np.max(ROC))
    cs = plt.contour(V_array,alts, ROC,levels=[0,1,2,3,4,5,6,7,7.5,8,8.5])
    
    plt.clabel(cs)
    he_cs = plt.contour(V_array, alts, He, levels=np.arange(0, 11000, 500),
        linestyles='dashed'
    )
    plt.clabel(he_cs)
    plt.xlabel("Velocity [m/s]")
    plt.ylabel("Altitude [m]")
    plt.show()


    he_flat = He.flatten()
    roc_flat = ROC.flatten()

    he_bins = np.arange(
        np.floor(np.min(he_flat)),
        np.ceil(np.max(he_flat)) + 50,
        50
    )
    
    roc_max = []

    for k in range(len(he_bins) - 1):

        mask = (
            (he_flat >= he_bins[k]) &
            (he_flat < he_bins[k + 1])
        )

        if np.any(mask):
            roc_max.append(np.max(roc_flat[mask]))
        else:
            roc_max.append(np.nan)

    roc_max = np.array(roc_max)

    # Energy-height coordinate corresponding to roc_max
    he_mid = 0.5 * (he_bins[:-1] + he_bins[1:])

    # Remove NaNs and non-positive climb rates
    valid = np.isfinite(roc_max) & (roc_max > 0)

    time_to_climb = np.trapz(
        1 / roc_max[valid],
        he_mid[valid]
    )

    print(f"Minimum time to climb = {time_to_climb:.1f} s")

    return time_to_climb
time_to_climb = h_v_plot(alts,1821*9.81,0.03376,24.23,10.18,0.7521,202000,1.38)   
print(time_to_climb)