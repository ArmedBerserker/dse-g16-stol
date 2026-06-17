import sys
import os

# Fix path FIRST, before any local imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import math
from classes.isa import *
from lookups.consts import *
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

clmax = 1.44
def stall_speed(W, C_l_max, rho, S):
    return math.sqrt((2 * W)/(rho * S * C_l_max))

def descent_perf(delta_T, W, cd0,S,AR,e,T, C_l_max):
    alt_m = [0, 500, 1000, 1500, 2000, 2500, 3000]

    C_L_list = []
    C_D_list = []
    V_list = []
    D_list = []
    aod_list = []
    rod_list =[]

    for alt in alt_m: 
        atmos_model = Atmosphere(alt, delta_T)
        rho = atmos_model.density[0]
        V = np.linspace(5, 110, 5000)
        V_list.append(V)

        V_s = stall_speed(W, C_l_max, rho, S)
        C_L = W / (0.5 * rho * V ** 2 * S)
        C_D = cd0 + C_L ** 2 / (np.pi * AR * e)
        D = C_D*0.5*rho*V**2*S
        L = W

        Cdmin = np.min(C_D)
        print(rho)
        print(Cdmin)
        k = 1/ (np.pi*e*AR)
        LD = C_L/ C_D
        idx_min_angle = np.argmax(LD)
        V_min_angle = np.sqrt(2/rho * np.sqrt(k/cd0)*W/S)
        LD_max = LD[idx_min_angle]
        min_angle = np.degrees(np.arctan(1/LD_max))
        range_m = alt*LD_max
        CL_bg = W / (0.5 *rho *V_min_angle**2*S)
        glide_ratio = range_m/alt

        print(
        f"{alt} m: "
        f"V = {V_min_angle:.2f} m/s, "
        f"theta = {min_angle:.2f} deg " 
        f"range = {range_m} m "
        f"stall = {V_s} m/s "
        f"CL_bg = {CL_bg} "
        f"glide ratio ={glide_ratio} "
        f"Vdif ={V_min_angle-V_s} " )
        
        C_L_list.append(C_L)
        C_D_list.append(C_D)
        D_list.append(D)
        aod_unpowered = np.arctan(C_D/C_L)

        #aod_powered = np.arcsin((D/W)- T/W) #calculate T later
       
        aod_list.append(aod_unpowered)
        
        rod = (D*V)/W
        #rod = V*np.sin(aod_unpowered)
        # rod = V / (C_L/C_D)
        rod_list.append(rod)
        idx_stall = np.argmin(np.abs(V - V_s))


        plt.plot(V, -rod, label=f"{alt} m")
        # plt.scatter(V_s, -rod[idx_stall], marker='x')

    # print(rod_list)
    plt.xlabel("Velocity [m/s]")
    plt.ylabel("Rate of Descent [m/s]")
    plt.grid(True)
    plt.legend()

    save_path=Path(__file__).parent / "performance_figures/descent_performance.png"
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    #plt.show()

    

def required_airspeed(aod_unpowered, C_L, W,S,rho):
    V_req = np.sqrt((2*np.cos(aod_unpowered)*W)/(rho*C_L*S))
    return V_req

# def descent_perf(delta_T, W, cd0,S,AR,e,T, C_l_max)
descent_perf(0, (1871*9.81), 0.0259, 31.4,10.18, 0.783, 200000,1.44)


#high AR negative effect on stall vs optimal
#Cd0 high negative effect on range and stall vs optimal
#weight does not do much
#oswald lower worse for range better for stall vs optimal
