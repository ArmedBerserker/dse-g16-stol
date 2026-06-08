import numpy as np
import matplotlib.pyplot as plt
from Propeller_performance import curve,D_ft

#conversion
kntstofps=1.68781
kgm3_to_slugft3=0.00194032
kgtolbs=2.20462
sqmetertofeetmeter=10.7639
mstofps=3.28084
d_lG=492.126 #[feet] max landing distance
Vs_L=50#stall speed in landing configuration at MTOW and sea level [knots]
rho=1.07896*kgm3_to_slugft3#density at landing altitude 2000t ISA +20°C #slug/feet^3
W_TO=1821*kgtolbs#[lbs]#change
W_LD=0.99*W_TO #[lbs]
g=32.2 #[feet/s**2]
CLmax_L=2.5#change
h_L=50 #[feet] CS23 requirement
V_sl_isa = 50 * kntstofps  # knots since close to MTOW and will go down with density
mu=0.4#assume hard turf
P_TO=270.886/2 #hp per engine
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
T_STATIC = 2*curve(D_ft,rho*1/kgm3_to_slugft3,P_TO)
T_BR=0.07*T_STATIC
#T_BR=-0.4*T_STATIC
S=269.1 #surface area of wing [feet2]
S_flaps=2.4*sqmetertofeetmeter*2 #surface of two side of flaps
AR=9
e= 0.7830160860700998


C_l0=2.5
C_lalpha=1 #we dont care anymore
C_d0=0.07517854464158887

h=(2.7)*mstofps #height of wing above ground [feet]
b=49.2#span of wing [feet]

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
    deltaf = 40  # flap deflection
    V_TD=V_BR=1.1*V_sl_isa
    h_f=0.1512*V_sl_isa**2*(1-np.cos(theta_app))
    S_A=(h_L-h_f)/np.tan(theta_app)#verified
    S_F=0.1512*V_sl_isa**2*(np.sin(theta_app))#verified
    S_FR=V_TD#verified
    #cl and cd estimations
    #alpha assumed 0 after the breaking sequence starts
    alpha_LDG=0
    C_L_ldg=C_l0+C_lalpha*alpha_LDG #Cl values after touchdown
    R=0.3
    delta1=179.32*R**4-111.6*R**3+28.29*R**2+2.3705*R-0.0089#check with flaps thcicness
    delta2= -3.9877e-12*deltaf**6 + 1.1685e-9*deltaf**5 - 1.2846e-7*deltaf**4 + 6.1742e-6*deltaf**3 - 9.89444e-5*deltaf**2 + 6.8324e-4*deltaf - 3.892e-4#check again
    delta_Cd=delta1*delta2*(S_flaps/S)
    C_di= C_L_ldg**2 / (np.pi*AR*e)

    C_di_ge=CDi_ground_effect(h,b,C_di)

    C_D_ldg= C_d0+delta_Cd+C_di_ge#Cd values after touchdown wasnt able to check that wiht the example

    D_lg=1/2*C_D_ldg*rho*(V_BR/np.sqrt(2))**2*S #verified
    L_ld=1/2*C_L_ldg*rho*(V_BR/np.sqrt(2))**2*S #verified

    S_BR=-(V_BR**2*W_LD)/(2*g*((T_BR-D_lg-mu*(W_LD-L_ld))))
    S_LDG=S_A+S_F+S_FR+S_BR
    S_GR=S_FR+S_BR
    return S_LDG,S_GR,D_lg,L_ld


def rho_at_altitude(h_ft):
    return rho * (1 - 6.875e-6 * h_ft) ** 4.2559

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
            rho_ref = rho_at_altitude(h)
            rad=np.radians(deg)
            V_0 = 1.1 * V_sl_isa
            V = 0
            S_FR = V_0
            #Assume thurst, drag and lift at Vbr/sqrt(2)
            # Static thrust max from where
            # Rough estimate for fixed-pitch/variable-pitch propeller
            T_STATIC = curve(D_ft,rho_ref*1/kgm3_to_slugft3,P_TO)
            T=0.07*T_STATIC  #formula from table 22.4 GA aircraft design for constant speed propeller
            S_LDG, S_GR, D_lg, L_ld= GORENBEEK_landing(W_LD, rho_ref)
            a = g / W_LD * (T - D_lg- mu * (W_LD * np.cos(rad) - L_ld) + W_LD * np.sin(rad))
            V_0 = V_0 * np.sqrt(rho /rho_ref)
            S_BR=(V**2-V_0**2)/(2*a)
            S_GR=S_BR+S_FR
            S.append(S_GR)
        axs[0].plot(S, altitudes, label=f"{deg}°")

    axs[0].set_title("Slope sensitivity")
    axs[0].set_xlabel("Ground run [ft]")
    axs[0].set_ylabel("Altitude [ft]")
    axs[0].grid()
    axs[0].legend()

    #temperature but will change
    temps = [0.9, 1.0, 1.1, 1.2] #temperature factor

    for t in temps:
        S = []
        for h in altitudes:
            rho_ref = rho_at_altitude(h)
            V_0 = 1.1 * V_sl_isa
            V = 0
            S_FR = V_0
            # Assume thurst, drag and lift at Vbr/sqrt(2)
            # Static thrust max from where
            # Rough estimate for fixed-pitch/variable-pitch propeller
            T_STATIC = curve(D_ft,rho*1/kgm3_to_slugft3,P_TO) # [lbs]
            T = 0.07 * T_STATIC  # formula from table 22.4 GA aircraft design for constant speed propeller
            S_LDG, S_GR, D_lg, L_ld = GORENBEEK_landing(W_LD, rho_ref)
            a = g / W_LD * (T - D_lg - mu * (W_LD - L_ld))
            V_0 = V_0 * np.sqrt(rho / rho_ref)
            S_BR = (V ** 2 - V_0 ** 2) / (2 * a)
            S_GR = S_BR + S_FR
            S.append(S_GR)

        axs[1].plot(altitudes, S, label=f"{t:.1f} thrust")

    axs[1].set_title("Temperature / thrust sensitivity not yet done")
    axs[1].set_xlabel("Altitude [ft]")
    axs[1].grid()
    axs[1].legend()


    weights = [0.8,0.9, 1.0, 1.1,1.2]

    for w in weights:
        S = []
        for h in altitudes:
            W_LDc=w * W_LD
            rho_ref = rho_at_altitude(h)
            V_0 = 1.1 * V_sl_isa
            V = 0
            S_FR = V_0
            # Assume thurst, drag and lift at Vbr/sqrt(2)
            # Static thrust max from where
            # Rough estimate for fixed-pitch/variable-pitch propeller
            T_STATIC = curve(D_ft,rho*1/kgm3_to_slugft3,P_TO) # [lbs]
            T = 0.07 * T_STATIC  # formula from table 22.4 GA aircraft design for constant speed propeller
            S_LDG, S_GR, D_lg, L_ld = GORENBEEK_landing(W_LDc, rho_ref)
            a = g / W_LDc * (T - D_lg - mu * (W_LDc - L_ld))
            V_0 = V_0 * np.sqrt(rho / rho_ref)
            S_BR = (V ** 2 - V_0 ** 2) / (2 * a)
            S_GR = S_BR + S_FR
            S.append(S_GR)
        axs[2].plot(S, altitudes, label=f"{w:.1f} MTOW")

    axs[2].set_title("Weight sensitivity")
    axs[2].set_xlabel("ground run [ft]")
    axs[2].grid()
    axs[2].legend()

    plt.tight_layout()
    plt.show()

sensitivity_altitude()

