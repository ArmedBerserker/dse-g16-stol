import sys
import os

# Fix path FIRST, before any local imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import math
import numpy as np
from classes.isa import *
from lookups.consts import *
import matplotlib.pyplot as plt

def take_off_vs_alt(delta_T, mtow, S, C_l_max, eta_p, P_s,n):
    alts = np.arange(0, 10000, 1)  # [ft]
    alts = alts * FT_TO_M
    s_g_list = []
    for a in alts:
        atmos_model = Atmosphere(a, delta_T)
        rho = atmos_model.density[0]
        V_s = math.sqrt((2 * mtow)/(rho * S * C_l_max))
        V_inf = 0.77*V_s
        P_a = eta_p*P_s*(rho/1.225)**n        
        T_a = P_a/V_inf
        s_g=(1.21*(mtow/S))/(g*rho*C_l_max*(T_a/mtow))
        s_g_list.append(s_g)
    alt_req = 2000*FT_TO_M  #metres 


    alt_req_line = np.full_like(s_g_list, alt_req)
    plt.xlabel("Take-off Ground Run [m]")
    plt.ylabel("Altitude [m]")
    plt.grid(True)
    plt.plot(s_g_list, alts)
    plt.plot(s_g_list, alt_req_line)
    plt.show()

    return s_g_list,alts

def take_off_vs_temp(mtow, S, C_l_max, eta_p, P_s, n):
    alts = np.arange(0, 10000, 1)  # [ft]
    alts = alts * FT_TO_M
    temps = np.linspace(-25,35,5000) #[K]
    s_g_list = []
    for T in temps:
        atmos_model = Atmosphere(0, T)
        rho = atmos_model.density[0]
        V_s = math.sqrt((2 * mtow)/(rho * S * C_l_max))
        V_inf = 0.77*V_s
        P_a = eta_p*P_s*(rho/1.225)**n        
        T_a = P_a/V_inf
        s_g=(1.21*(mtow/S))/(g*rho*C_l_max*(T_a/mtow))
        s_g_list.append(s_g)
    #alt_req = 2000*FT_TO_M  #metres 
    temp_req = 288.15+20 # [K]
    print(max(temps))

    temp_req_line = np.full_like(s_g_list, temp_req)
    plt.xlabel("Take-off Ground Run [m]")
    plt.ylabel("Temperature [K]")
    plt.grid(True)
    plt.plot(s_g_list, temps+288.15)
    plt.plot(s_g_list, temp_req_line)
    plt.show()
    return(s_g_list, temps)

def take_off_alt_temp(mtow, S, C_l_max, eta_p, P_s, n):
    alts = np.linspace(0, 10000, 200) * FT_TO_M  # altitude sweep (m)
    delta_Ts = [-25, -10, 0, 10, 20,35]  # ISA deviation [K]

    plt.figure()
    for dT in delta_Ts:

        s_g_list = []

        for a in alts:

            atmos_model = Atmosphere(a, dT)
            rho = atmos_model.density[0]

            V_s = math.sqrt((2 * mtow) / (rho * S * C_l_max))
            V_inf = 0.77 * V_s

            P_a = eta_p * P_s * (rho / 1.225)**n
            T_a = P_a / V_inf

            s_g = (1.21 * (mtow / S)) / (g * rho * C_l_max * (T_a / mtow))
            s_g_list.append(s_g)

        plt.plot(s_g_list, alts, label=f"ΔT = {dT} K")
    alt_req = 2000 * FT_TO_M  # m

    x_line = np.linspace(0, max([max(s) for s in [s_g_list]]), 200)
    alt_req_line = np.full_like(x_line, alt_req)

    plt.plot(x_line, alt_req_line, '--', label="Required altitude")
    plt.ylabel("Altitude [m]")
    plt.xlabel("Take-off Ground Run [m]")
    plt.grid(True)
    plt.legend()
    plt.show()

def take_off_dist_mass (delta_T, mtow, S, C_l_max, eta_p, P_st, P_sp,n_t, n_p):
    s_g_list_turbo = []
    s_g_list_piston = []

    atmos_model = Atmosphere(609.6, delta_T)
    rho = atmos_model.density[0]

    weights = np.arange(8500, mtow+1, 1)
    for w in weights:
        V_s = math.sqrt((2 * w)/(rho * S * C_l_max))
        V_inf = 0.77*V_s

        P_a_t = eta_p*P_st*(rho/1.225)**n_t        
        T_a_t = P_a_t/V_inf
        s_g_t=(1.21*(w/S))/(g*rho*C_l_max*(T_a_t/w))
        s_g_list_turbo.append(s_g_t)

        P_a_p = eta_p*P_sp*(rho/1.225)**n_p        
        T_a_p = P_a_p/V_inf
        s_g_p=(1.21*(w/S))/(g*rho*C_l_max*(T_a_p/w))
        s_g_list_piston.append(s_g_p)
   
    plt.xlabel("Take-off Ground Run [m]")
    plt.ylabel("Take-off Weight [N]")
    plt.grid(True)
    #plt.plot(s_g_list_piston, weights)
    #plt.plot(s_g_list_turbo, weights)
    plt.plot(s_g_list_piston, weights, color = 'red',label="Boosted Piston")
    plt.plot(s_g_list_turbo, weights, label="Boosted Turboprop")
    # MTOW horizontal line
    plt.axhline(mtow,
            linestyle='--',
            color='black',
            label='MTOW')
    plt.legend()
    plt.show()

     



alt = 609.6 # [m]
delta_T = 20  # [K]
mtow = 2000 * g # [N]
C_l_max = 2.5 #take-off 
S = 25.675  # [m^2]
eta_p = 0.8 # [-]
#n = 0.8 # 1 for piston / 0.8 for turboprop
P_st = 550000 #[W] 550000 shaft power for turboprop, 202000 for piston 
P_sp = 202000 #W
n_t = 0.8
n_p = 1
s_g_list, alts = take_off_vs_alt(delta_T,mtow, S, C_l_max,eta_p, P_st, n_t)
#print(s_g_list)
s_g_list, temps = take_off_vs_temp(mtow, S, C_l_max, eta_p, P_st, n_t)

take_off_alt_temp(mtow, S, C_l_max, eta_p, P_sp, n_p)
take_off_dist_mass(delta_T, mtow, S, C_l_max, eta_p, P_st, P_sp,n_t, n_p)
    