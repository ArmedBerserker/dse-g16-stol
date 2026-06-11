import numpy as np
import matplotlib.pyplot as plt
from Propeller_performance import curve,D_ft
from classes.isa import *
#conversion
kntstofps=1.68781
kgm3_to_slugft3=0.00194032
kgtolbs=2.20462
sqmetertofeetmeter=10.7639
mstofps=3.28084

d_lG=492.126 #[feet] max landing distance
Vs_L=50#stall speed in landing configuration at MTOW and sea level [knots]
rho=1.07896*kgm3_to_slugft3#density at landing altitude 2000t ISA +20°C #slug/feet^3
W_TO=1839*kgtolbs#[lbs]#change
W_LD=0.99*W_TO #[lbs]
g=32.2 #[feet/s**2]
CLmax_L=1.7#change
h_L=50 #[feet] CS23 requirement
V_sl_isa = 50 * kntstofps  # knots since close to MTOW and will go down with density
mu=0.4#assume hard turf
P_TO=314/2 #hp per engine
delta_T=20 #Kelvin
"""
might not be useful but maybe 
S_LG=0.265*Vs_L**2 #ground run after touchdown [feet]

S_L=1.938*S_LG #Total landing distance [feet]
"""

#ROSKAM first estimation for FAR23 all those number come from approximation in book 7
#gamma=(D-T)/W_LD
#gamma=C_D/C_L
#global variables


def Roskam_landing():
    d_l=1.938*d_lG #Total landing distance [feet]
    gamma=0.1
    deltan=0.1
    a=0.33*g
    f_land=1 #for small aircraft
    #V_sl_isa=np.sqrt(d_l/0.265)*kntstofps #knts max stall speed for landing performance at isa+20

    V_A=1.3*V_sl_isa
    V_TD=V_A*(1-(gamma**2/deltan))**(1/2)

    s_AIR=(1/gamma)*((V_A**2-V_TD**2)/(2*g)+h_L)
    s_LG=(V_TD**2)/(2*a)
    s_L=s_AIR+s_LG

    Wingloading_L= (s_L/(h_L*f_land)-10)*(h_L*rho*g*CLmax_L)/(1.52/(a/g)+1.69)#winfg loading at landing lb/ft2
    return s_L,s_LG,s_AIR



#METHOD 2 GORENBEEK
#unknowns
#efficiency=0.65 #efficiency at VB/Sqrt2
#P_BHP=0.07*270.886 #power at VB/Sqrt2
T_static1,V,eff1= curve(D_ft,rho*1/kgm3_to_slugft3,P_TO)
T_static=2*T_static1
T_BR=0.07*T_static

#T_BR=-0.4*T_STATIC reverse thrust
S_w=333.6812 #surface area of wing [feet2]
S_flaps=2.8*sqmetertofeetmeter*2 #surface of two side of flaps not needed because naomi calculated drag from flaps
AR=9
e= 0.783
C_l0=2.5 #not needed
C_lalpha=1 #we dont care anymore
C_d0=0.08777
deltaf =60  # flap deflection
h=2.5 #height of wing above ground [m]
b=16.8 #span of wing [m]
V_STO_L = np.sqrt(2 * W_TO / (rho * S_w * CLmax_L))  # Stall speed during take off [feet/s]
theta_app = np.radians(3)  # radians or 6 check both

#calculations

def CDi_ground_effect(h,b,c_Di):
    hb = h / b
    if hb < 0.033:
        return 0
    elif hb < 0.33:
        return (1 - (1 - 1.32 * hb) / (1.05 + 7.4 * hb)) * c_Di
    else:
        return 0

def GORENBEEK_landing(W_LD,rho):

    V_TD=V_BR=1.1*V_STO_L
    h_f=0.1512*V_STO_L**2*(1-np.cos(theta_app))#feet
    S_A=(h_L-h_f)/np.tan(theta_app)#verified
    S_F=0.1512*V_STO_L**2*(np.sin(theta_app))#verified
    S_FR=V_TD#verified
    #cl and cd estimations
    #alpha assumed 0 after the breaking sequence starts*
    C_Lland=2*W_LD/(rho*S_w*(1.3*V_STO_L)**2)#double check
    #C_L_ldg = C_Lland  # without speed break
    C_L_ldg=C_Lland-0.4 #speed break
    R=0.25 #chord
    delta1=-21.090*R**3+14.091*R**2+3.165*R-0.00103#check with flaps thcicness
    delta2= -3.795e-7*deltaf**3+5.387e-5*deltaf**2 -6.843e-4*deltaf-1.4729e-3 #check again
    delta_Cd=delta1*delta2*(S_flaps/S_w)
    C_di= C_L_ldg**2 / (np.pi*AR*e)
    C_di_ge=CDi_ground_effect(h,b,C_di)
    C_D_ldg = C_d0 + C_di_ge   # without speed breaks
    C_D_ldg= C_d0+C_di_ge+ 0.04#Cd values after touchdown wasnt able to check that wiht the example speed brake add 0.4

    D_lg=1/2*C_D_ldg*rho*(V_BR/np.sqrt(2))**2*S_w #verified
    L_ld=1/2*C_L_ldg*rho*(V_BR/np.sqrt(2))**2*S_w #verified

    S_BR=-(V_BR**2*W_LD)/(2*g*((T_BR-D_lg-mu*(W_LD-L_ld))))
    #a=g*(T_BR - D_lg - mu * (W_LD - L_ld))/W_LD
    S_LDG=S_A+S_F+S_FR+S_BR
    S_GR=S_FR+S_BR
    return S_LDG,S_GR,D_lg,L_ld


print("Results ROSKAM")
s_L,s_LG,s_AIR=Roskam_landing()
print("Ground run", s_LG)

print("Results GORENBEEK")
S_LDG,S_GR,D_lg,L_ld=GORENBEEK_landing(W_LD,rho)
print("Landing distance", S_LDG,S_GR)


def sensitivity_altitude():
    rho = 1.07896 * kgm3_to_slugft3
    # density at landing altitude 2000t ISA +20°C #slug/feet^3

    altitudes = np.arange(0, 10000, 30)

    fig, axs = plt.subplots(1, 3, figsize=(18, 5))

    #steepness
    slopes_deg = [0, 1, 2, 3, 4]

    for deg in slopes_deg:
        S = []
        for h in altitudes:
            atmos_model = Atmosphere(h/mstofps, delta_T)
            rho_ref = atmos_model.density[0]*kgm3_to_slugft3
            rad=np.radians(deg)
            V_0 = 1.1 * V_STO_L
            V = 0
            S_FR = V_0
            #Assume thurst, drag and lift at Vbr/sqrt(2)
            # Static thrust max from where
            # Rough estimate for fixed-pitch/variable-pitch propeller
            T=0.07*T_static  #formula from table 22.4 GA aircraft design for constant speed propeller
            S_LDG, S_GR, D_lg, L_ld= GORENBEEK_landing(W_LD, rho_ref)
            a = g / W_LD * (T - D_lg- mu * (W_LD * np.cos(rad) - L_ld) + W_LD * np.sin(rad))
            V_0 = V_0 * np.sqrt(rho /rho_ref)
            S_BR=(V**2-V_0**2)/(2*a)
            S_GR=S_BR+S_FR
            S.append(S_GR)
        axs[0].plot(S, altitudes, label=f"{deg}°")
    axs[0].axvline(d_lG, linestyle="--", color="red", label=f"Field limit = {d_lG:.0f} ft")
    axs[0].set_title("Slope sensitivity")
    axs[0].set_xlabel("Ground run [ft]")
    axs[0].set_ylabel("Altitude [ft]")
    axs[0].grid()
    axs[0].legend()

    #temperature
    delta_T_list = np.arange(0, 41, 10)

    for Temp in delta_T_list:
        S = []
        for h in altitudes:
            atmos_model = Atmosphere(h/mstofps, Temp)
            rho_ref = atmos_model.density[0]*kgm3_to_slugft3
            V_0 = 1.1 * V_STO_L
            V = 0
            S_FR = V_0
            # Assume thurst, drag and lift at Vbr/sqrt(2)
            # Static thrust max from where
            # Rough estimate for fixed-pitch/variable-pitch propeller
            T = 0.07 * T_static  # formula from table 22.4 GA aircraft design for constant speed propeller
            S_LDG, S_GR, D_lg, L_ld = GORENBEEK_landing(W_LD, rho_ref)
            a = g / W_LD * (T - D_lg - mu * (W_LD - L_ld))
            V_0 = V_0 * np.sqrt(rho / rho_ref)
            S_BR = (V ** 2 - V_0 ** 2) / (2 * a)
            S_GR = S_BR + S_FR
            S.append(S_GR)

        axs[1].plot(S, altitudes, label=f"{Temp:.1f} temperature")
    axs[1].axvline(d_lG, linestyle="--", color="red", label=f"Field limit = {d_lG:.0f} ft")
    axs[1].set_title("Altitude VS temperature [ft]")
    axs[1].set_xlabel("Ground run")
    axs[1].set_ylabel("Altitude [ft]")
    axs[1].grid()
    axs[1].legend()


    weights = [0.8,0.9, 1.0, 1.1,1.2]

    for w in weights:
        S = []
        for h in altitudes:
            W_LDc=w * W_LD
            atmos_model = Atmosphere(h / mstofps, delta_T)
            rho_ref = atmos_model.density[0]*kgm3_to_slugft3
            V_0 = 1.1 * V_STO_L
            V = 0
            S_FR = V_0
            # Assume thurst, drag and lift at Vbr/sqrt(2)
            # Static thrust max from where
            # Rough estimate for fixed-pitch/variable-pitch propeller
            T = 0.07 * T_static  # formula from table 22.4 GA aircraft design for constant speed propeller
            S_LDG, S_GR, D_lg, L_ld = GORENBEEK_landing(W_LDc, rho_ref)
            a = g / W_LDc * (T - D_lg - mu * (W_LDc - L_ld))
            V_0 = V_0 * np.sqrt(rho / rho_ref)
            S_BR = (V ** 2 - V_0 ** 2) / (2 * a)
            S_GR = S_BR + S_FR
            S.append(S_GR)
        axs[2].plot(S, altitudes, label=f"{w:.1f} MTOW")
    axs[2].axvline(d_lG, linestyle="--", color="red", label=f"Field limit = {d_lG:.0f} ft")
    axs[2].set_title("Weight sensitivity")
    axs[2].set_xlabel("ground run [ft]")
    axs[2].set_ylabel("Altitude [ft]")
    axs[2].grid()
    axs[2].legend()

    plt.tight_layout()
    plt.show()

sensitivity_altitude()

