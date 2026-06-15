import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import pandas as pd
from classes.isa import *
from scipy.optimize import curve_fit
#propeller chosen https://www.propellor.com/ap431hapf-snl68e
graphdata = pd.read_csv("ctcp")

#check units
kgtoslug = 0.00194032
ktastofeet= 1.68781
inchestofeet=1/12
desnityconversion=0.00194032
fpskph=1.09728
mstofps=3.28084
kgtolbs=2.20462
g=32.2
rho_TO=1.07896*desnityconversion
rho_cruise=0.7812*desnityconversion
T_TO=304.11 #take off altitude 2000t ISA +20°C
T_CR=271.3 #isa 8000ft
Vcruise= 132*ktastofeet #[f/s]
W_TO=1809*kgtolbs#[lbs]#change
S=236.80603#surface area of wing [feet2]
C_LmaxTO=2.55#max take off cl
V_STO=np.sqrt(2*W_TO/(rho_TO*S*C_LmaxTO))
Vto=1.1*V_STO*ktastofeet #[f/s]
P_to=254.7942/2#hp change per engine
Mtipmax=0.80
z=(1.7+0.65-0.21) #distance from center of engine to ground meter

#depends on engine
P_hbpTO=160  #of engine at 3258
P_hbpCR=137 #of engine at 3100

Pmax=160 #[hp] maximum power per engine



graphtorrenbeek=np.sqrt(Vcruise*fpskph*P_to)
maxdiameter=2*mstofps
# print(graphtorrenbeek)

value=4.1 #find on graph toreenbeedk hp/feet2



def speed_of_sound(T_K):
    a=np.sqrt(1.4 * 287.0 * T_K)
    return a*mstofps #fps
def D_max_from_clearance(z):
    return 2.0 * (z - 0.18)*mstofps #ft

def prop_diameter_Torenbeek (P_bhp,value):
    return (np.sqrt(P_bhp/value)) #ft

def prop_diameter_power_limit(Pmax,n):
    Pblmax=4.8
    Pblmin=2.8
    return (((4*Pmax)/(np.pi*n*Pblmin))**(1/2),((4*Pmax)/(np.pi*n*Pblmax))**(1/2)) #ft pMax is the maximum power from enfine during take off and n is number of propeller

#not really useful in the end
def prop_diameter_blade(P_bhp):
    #Phbp engine power
    return {
        2: 20.4 * P_bhp**0.25*inchestofeet,
        3: 19.2 * P_bhp**0.25*inchestofeet,
        4: 18 * P_bhp**0.25*inchestofeet
    } #feet


#find RPM

def RPM_from_tip_limit(V_fps,T_K,D_ft):
    afps=speed_of_sound(T_K)
    Vtip_max_fps=Mtipmax*afps
    Vrot = np.sqrt((Vtip_max_fps**2 - V_fps**2))
    return 60 * Vrot / (np.pi * D_ft) #RPM


def advance_ratio(RPM,D_ft,V):
    return 60*V/(RPM*D_ft)





"""
Sc=np.pi*(0.63/2)**2 #max area of nacelle
J=J*(1-0.329*Sc/Dtor**2) #change diameter based on the one selected
"""

def P_to_RPM(P):
    #if need add the capacitor
    Power = np.array([48,60,78,104,126,142,160])
    RPM = np.array([3000,3500,4000,4500,5000,5500,5800])

    RPMnew = RPM/2.54

    # RPM as a function of Power
    spline = CubicSpline(Power, RPMnew)
    xfine = np.linspace(Power.min(), Power.max(), 1000)
    yfine = spline(xfine)


    return float(spline(P))

def cp_calculation(D_ft, rho_kgm3, P_bhp,RPM):
    rho_slug = rho_kgm3*desnityconversion    # slug/ft³
    n = RPM / 60.0                          # rev/s
    P_fps = P_bhp * 550.0                   # ft·lbf/s
    return P_fps / (rho_slug * n**3 * D_ft**5)


def eff(D_ft,rho_kgm3,P_bhp,J):
    RPM=P_to_RPM(P_bhp)
    Cp=cp_calculation(D_ft, rho_kgm3, P_bhp,RPM)

    x = graphdata.iloc[:, 0] #cp
    y = graphdata.iloc[:, 1] #ct/cp
    spline = CubicSpline(x, y)
    ctcp = spline(Cp)


    """
    check the graphs and interpolation
    xfine = np.linspace(min(x), max(x), 1000)
    yfine = spline(xfine)
    plt.scatter(x, y, label="Data")
    plt.plot(xfine, yfine, label="Spline")
    plt.legend()
    plt.grid()
    plt.show()
    """
    # print("Cp range in chart:", x.min(), x.max())
    # print("Cp used:", Cp)

    # print("Y range in chart:", y.min(), y.max())
    efficiency=ctcp*J
    return efficiency,ctcp,Cp

def curve( D_ft,rho_kgm3,P_bhp):
    #power shaft in horsepower per singel engine
    V = np.linspace(20, 150, 1000)  # knots
    Vnew = V * ktastofeet

    RPM = P_to_RPM(P_bhp)
    n = RPM / 60
    J=Vnew/(n*D_ft)
    Sc = np.pi * (0.63*mstofps / 2) ** 2  # max area of nacelle
    J = J * (1 - 0.329 * Sc / D_ft ** 2)  # change diameter based on the one selected

    lam=J
    eff2,ctcp,Cp=eff(D_ft,rho_kgm3,P_bhp,J)
    T=ctcp*Cp*rho_kgm3*desnityconversion*n**2*(D_ft)**4

    q=1/2*rho_kgm3*desnityconversion*(Vnew)**2

    A=np.pi *(D_ft/2)**2

    #https://www.fzt.haw-hamburg.de/pers/Scholz/transfer/Airport2030_TN_Propeller-Efficiency_13-08-12_SLZ.pdf
    num =2*(1 -lam ** 2 * np.log(1 + 1/(lam ** 2)) )
    den = 1 + np.sqrt(1+T/(q*A))-2* lam **2* np.log(1 + 1/(lam ** 2))
    eff11 =num / den#new because other gave me linear relationships
    eff12=2/(1+np.sqrt(1+T/(q*A))) #other flight mechanics book
    eff1=(eff11+eff12)/2
    P_useful = eff1 * P_bhp
    #T_static= ctcp*550*P_bhp/(n*D_ft)#t_static according to Raymer
    T_static=0.85*(P_bhp*550)**(2/3)*(2*rho_kgm3*desnityconversion*A)**(1/3) #one engine and overestimate due to the fact they dont take into account blockage
    P_useful_w = T*Vnew* 745.7*2 #both engines in watts
    P1=P_useful*745.7*2
    V_ms = Vnew * 0.3048

    # print("RPM =", RPM)
    # print("J range =", J.min(), J.max())
    # print("Cp =", Cp)
    #print("ctcp =", ctcp)
    #print("efficiency =", eff2)
    #print("eff1 max =", np.max(eff1))
    """
    V = np.linspace(20, 150, 1000) #knots
    Vnew=V*ktastofeet
    plt.figure()
    plt.plot( V_ms , eff1)
    plt.xlabel("Velocity [ms]")
    plt.ylabel("Propeller efficiency")
    plt.grid()
    print(V_ms,eff1)

    plt.figure()
    plt.plot(V_ms, P_useful)
    plt.xlabel("Velocity [ms]")
    plt.ylabel("Useful power [w]")
    plt.grid()

    plt.show()
    """
    return T_static,Vnew,eff1
    #return V_ms,P_useful_w

D_ft=5.66667
T_static,Vnew,eff1=curve(D_ft,1.02,160)

delta_T=20
altitudes = np.arange(0, 8000, 50)
P_available = []

"""
for alt in altitudes:
    atmos_model = Atmosphere(alt / mstofps, delta_T)
    rho_local = atmos_model.density[0]
    T_static, Vnew, eff1, P_useful_w,T= curve(D_ft, rho_local, 160)
    P_available.append(np.max(P_useful_w))

plt.figure()
plt.plot(P_available, altitudes)
plt.xlabel("Available power [W]")
plt.ylabel("Altitude [ft]")
plt.grid()
plt.show()
T_available = []
for alt in altitudes:
    atmos_model = Atmosphere(alt / mstofps, delta_T)
    rho_local = atmos_model.density[0]
    T_static, Vnew, eff1, P_useful_w,T= curve(D_ft, rho_local, 160)
    T_available.append(np.max(T))

plt.figure()
plt.plot(T_available, altitudes)
plt.xlabel("Available power [W]")
plt.ylabel("Altitude [ft]")
plt.grid()
plt.show()
"""
