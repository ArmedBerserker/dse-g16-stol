import matplotlib.pyplot as plt
from Propeller_performance import curve, D_ft
import numpy as np
from classes.isa import *
# UNIT CONVERSION
kntstofps = 1.68781
kgtoslug = 0.00194032
kgtolbs = 2.20462

# GLOBAL VARIABLE
rho = 1.07896 * kgtoslug  # density at take off altitude 2000ft ISA +20°C #slug/feet^3
rho_SL = 1.225 * kgtoslug
T_ISA_K = 304.188  # temperature at take off altitude 2000ft ISA +20°C #Kelvin
theta = rho / rho_SL  # ratio of density
W_TO = 1988.8 * kgtolbs  # [lbs]#change
S = 269.1  # surface area of wing [feet2]
D2 = 378  # get value from fomula from naomi
C_LmaxTO = 2.5  # max take off lift coefficient
C_LTO = C_LmaxTO / 1.21  # Lift coefficient at lift off
h = 2.7  # height of wing above ground [feet]
b = 18.2  # span of wing [feet]
c_Di = 0.04  # lift induced drag coefficient out of ground effect
C_D0 = 0.04  # cD0 take off run
delta_T = 20  # [K]
# propulsion
P_TO = 289.661  # max shaft power in horsepower during take off all engines operating
efficiency = 0.60  # efficiency of motor at V_LOF/Sqrt2
efficiencyarray = 0.60  # give actual array in function of speed change this
P_bhp = 160  # max hp per engine during take off
# performance
climbangle = np.radians(3.4)  # climb angle in radians
d_TO = 656.168  # [feet] max take-off distance
S_TOG = d_TO  # distance take off
h_TO = 50  # feet
g = 32.2  # [feet/s**2]
mu = 0.08  # ground roll coefficient on wet grass without braking

V_STO = np.sqrt(2 * W_TO / (rho * S * C_LmaxTO))  # Stall speed during take off [feet/s]


# ROSKAM METHOD
def Roskam_TO(C_D0, C_LmaxTO, P_bhp, S_TOG):
    # S_TO = 1.66 * S_TOG  # total take off distance [feet]
    # solve S_TOG=4.9*TOP+0.009*TOP**2 and i checked and deleted the negative solution
    coeff = [0.009, 4.9, -S_TOG]
    TOP = np.max(np.roots(coeff))  # lbs**2/(ft**2hp) only useful for the graph

    # ground run approach
    ug = 0.15  # roskam for soft ground assume worse case middle
    muR = ug + 0.72 * (C_D0 / C_LmaxTO)
    V_LOF = V_STO * 1.2
    # since variable pitch use 5.75, 5.4 since calculated value for graph in propeller sizing value in hp/feet2
    T = 5.75 * P_bhp * ((theta * 2 * D_ft ** 2 / P_bhp) ** (1 / 3))  # lbs average thrust during take off
    TW = T / (W_TO)
    S_GR = (V_LOF ** 2 / (2 * g)) / (TW - muR)
    return S_GR, TOP, T


# ground effcet
def CDi_ground_effect(h, b, c_Di):
    hb = h / b
    if hb < 0.033:
        return 0
    elif hb < 0.33:
        phi = (1 - (1 - 1.32 * hb) / (1.05 + 7.4 * hb))
        return phi * c_Di
    else:
        print("switch method")
        return 0


# gorenbeek if propeller or jet
def GORENBEEK_TO_1(P_TO, b, h, c_Di, S,
                   C_D0):  # didnt amnage to check everything but thrust was checked and probs won't use it
    # for small powered plane VR znd VLOF are assumed to be the same since it lifts off as soon as it rotates
    V_R = 1.1 * V_STO
    V_LOF = 1.556 * np.sqrt(W_TO / (rho * S * C_LmaxTO))
    V = V_LOF / np.sqrt(2)  # speed in feet/s at Vlof/sqrt2 to not over or under estimate
    # estimate lift induced drag due to ground effect
    CDi_IGE = CDi_ground_effect(h, b, c_Di)

    L = 1 / 2 * rho * V ** 2 * S * C_LTO
    D = 1 / 2 * rho * V ** 2 * S * (C_D0 + CDi_IGE)
    T = efficiency * (550 * P_TO) / (V)
    s_GR = (V_LOF ** 2 * W_TO) / (2 * g * (T - D - mu * (W_TO - L)))
    return T, s_GR


# if piston engine faster approach
def GORENBEEK_TO_2():
    # for small powered plane VR znd VLOF are assumed to be the same since it lifts off as soon as it rotates
    V_R = 1.1025 * V_STO
    V_R = 1.1 * V_STO
    V_LOF = 1.556 * np.sqrt(W_TO / (rho * S * C_LmaxTO))

    # estimate lift induced drag due to ground effect
    CDi_IGE = CDi_ground_effect(h, b, c_Di)
    C_DTO = C_D0 + CDi_IGE
    s_GR = V_LOF ** 2 * W_TO / (50051 * efficiency * P_TO / V_LOF + 16.09 * rho * V_LOF ** 2 * S * (mu * C_LTO - C_DTO))
    return s_GR


# serious method numerical method
def GORENBEEK_TO_3(dt=0.05, max_time=200):
    # initial condition
    V = 0.1
    s = 0.0
    t = 0.0

    V_LOF = 1.556 * np.sqrt(W_TO / (rho * S * C_LmaxTO))  # [ft/s]
    V_TR = 1.15 * V_STO  # [ft/s] transition speed

    CDi_IGE = CDi_ground_effect(h, b, c_Di)

    # Static thrust max from where in propeller
    T_STATIC = curve(D_ft, 1.07896, P_TO)  # [lbs]²    4

    s_LO = None
    V_LO = None
    T_TR = None
    V1 = None
    lifted = False

    trajectory = []

    while t < max_time:
        q = 0.5 * rho * V ** 2
        L = q * S * C_LTO
        D = q * S * (C_D0 + CDi_IGE)
        # stop it from going to infinity
        if V < 10.0:
            T = T_STATIC
        else:
            T = efficiencyarray * 550.0 * P_TO / V

        a = g / W_TO * (T - D - mu * (W_TO - L))

        trajectory.append([t, V, s, a, L, T])

        if not lifted and V >= V_LOF:
            s_LO = s
            V_LO = V
            lifted = True

        if T_TR is None and V >= V_TR:
            T_TR = T

        # Stop integration at 1.3 v stall since should not be ground roll anymore
        if V >= 1.3 * V_STO:
            break

        # Integration
        V = max(V + a * dt, 0.1)
        s = s + V * dt + 0.5 * a * dt ** 2
        t += dt

    return s, V, t, np.array(trajectory), T_TR, V1, s_LO, V_LO


def rest(T_TR):
    V_LOF = 1.556 * np.sqrt(W_TO / (rho * S * C_LmaxTO))
    S_ROT = V_LOF  # rotation distance if small aircraft is assumed to be 1sec
    V_TR = 1.15 * V_STO
    T = T_TR
    CL_TR = 2 * W_TO / (rho * V_TR ** 2 * S)
    CDi_IGE = CDi_ground_effect(h, b, c_Di)  # still there because less than 1 wing span
    CD_TR = C_D0 + CDi_IGE
    LD = CL_TR / CD_TR

    anglemethod = np.arcsin(T / W_TO - 1 / (LD))  # CS 23 approach
    # check with climb angle is restrictive

    if climbangle > anglemethod:
        angle = climbangle
        print("OUR constrain dominant for obstacle clearance ")

    else:
        angle = anglemethod
        print("CS23 constrain dominant for obstacle clearance")

    S_TR = 0.2156 * V_STO ** 2 * np.sin(angle)  # distance before transition to straight flight
    h_TR = 0.2156 * V_STO ** 2 * (1 - np.cos(angle))  # transition height  [feet]
    S_C = (h_TO - h_TR) / (np.tan(angle))  # climb disantce until top of obstacle [ft]

    if S_C < 0:
        S_C = 0
        print("obstacle cleared during transition")
    S_OBS = S_TR + S_ROT + S_C
    return S_OBS, h_TR


def ground_run_for(rho_local, w_local, slope_rad, dt=0.05, max_time=200.0):
    c_lmax_local = C_LmaxTO
    v_lof_local = 1.556 * np.sqrt(w_local / (rho_local * S * c_lmax_local))
    rho_SI = rho_local / kgtoslug
    CDi_IGE = CDi_ground_effect(h, b, c_Di)
    T_STATIC = curve(D_ft, rho_SI, P_TO)  # [lbs] #check if  propeller change

    V, s, t = 0.1, 0.0, 0.0
    s_LO = np.nan

    while t < max_time:
        q = 0.5 * rho_local * V ** 2
        L = q * S * C_LTO
        D = q * S * (C_D0 + CDi_IGE)
        T = T_STATIC if V < 10.0 else efficiencyarray * 550.0 * P_TO / V

        a = (g / w_local) * (T - D - mu * (w_local - L) - w_local * np.sin(slope_rad))

        V_new = max(V + a * dt, 0.1)
        s_new = s + V * dt + 0.5 * a * dt ** 2

        if V_new >= v_lof_local:
            # interpolate exact lift-off distance instead of snapping to step
            frac = (v_lof_local - V) / (V_new - V) if V_new != V else 1.0
            s_LO = s + frac * (s_new - s)
            break

        V, s = V_new, s_new
        t += dt
    return s_LO


def rho_at_altitude(h_ft, delta_T_K=20):
    # Density at h_ft above sea level (ISA)
    rho_isa = rho_SL * (1 - 6.875e-6 * h_ft) ** 4.2559
    # Temperature correction (hot day)
    T_isa = 288.15 * (1 - 6.875e-6 * h_ft)  # ISA temperature at h_ft [K]
    rho_hot = rho_isa * (T_isa / (T_isa + delta_T_K))
    return rho_hot


# other approach to find BFL is when the distabce to stop and go equals the distance to stop and keeo going
def acceleratego(V1, dt=0.05, max_time=200):
    V, s, t = 0.1, 0.0, 0.0
    V_LOF = 1.556 * np.sqrt(W_TO / (rho * S * C_LmaxTO))  # [ft/s]
    CDi_IGE = CDi_ground_effect(h, b, c_Di)
    T_STATIC = curve(D_ft, 1.07896, P_TO)  # lbs
    V_TR = 1.15 * V_STO
    T_TR = None
    # accelerate
    while V < V1 and t < max_time:
        q = 0.5 * rho * V ** 2
        L = min(q * S * C_LTO, W_TO)
        D = q * S * (C_D0 + CDi_IGE)
        # stop it from going to infinity
        if V < 10.0:
            T = T_STATIC
        else:
            T = efficiencyarray * 550.0 * P_TO / V

        a = g / W_TO * (T - D - mu * (W_TO - L))
        # Integration
        V = max(V + a * dt, 0.1)
        s = s + V * dt + 0.5 * a * dt ** 2
        t += dt

    # keep going with half thrust until clear the obstacle is the right distance
    while V < V_LOF and t < max_time:
        q = 0.5 * rho * V ** 2
        L = min(q * S * C_LTO, W_TO)
        D = q * S * (C_D0 + CDi_IGE)
        # stop it from going to infinity
        if V < 10.0:
            T = T_STATIC / 2
        else:
            T = efficiencyarray * 550.0 * P_TO / (V * 2)

        a = g / W_TO * (T - D - mu * (W_TO - L))
        # Integration
        if T_TR is None and V >= V_TR:  # check Transisiton thrust
            T_TR = T
        V = max(V + a * dt, 0.1)
        s = s + V * dt + 0.5 * a * dt ** 2
        t += dt

    if T_TR is None:
        T_TR = efficiencyarray * 550.0 * (P_TO / 2) / (V_TR)
    S_OBS, h_TR = rest(T_TR)  # ask paul
    s += S_OBS
    return s


def acceleratestop(V1, dt=0.05, max_time=200):
    # accelerate
    V, s, t, t_rec = 0.1, 0.0, 0.0, 0.0
    V_LOF = 1.556 * np.sqrt(W_TO / (rho * S * C_LmaxTO))  # [ft/s]
    CDi_IGE = CDi_ground_effect(h, b, c_Di)
    T_STATIC = curve(D_ft, 1.07896, P_TO)  # lbs
    reaction = 2
    mu_break = 0.4  # assume hard turf for braking
    # accelerate
    while V < V1 and t < max_time:
        q = 0.5 * rho * V ** 2
        L = min(q * S * C_LTO, W_TO)
        D = q * S * (C_D0 + CDi_IGE)
        # stop it from going to infinity
        if V < 10.0:
            T = T_STATIC
        else:
            T = efficiencyarray * 550.0 * P_TO / V

        a = g / W_TO * (T - D - mu * (W_TO - L))
        # Integration
        V = V + a * dt
        s = s + V * dt + 0.5 * a * dt ** 2
        t += dt

    # reaction time 3 s t rec + activation of breaking device following CS23
    while t_rec < reaction and t < max_time:
        q = 0.5 * rho * V ** 2
        L = min(q * S * C_LTO, W_TO)
        D = q * S * (C_D0 + CDi_IGE)
        # stop it from going to infinity
        if V < 10.0:
            T = T_STATIC / 2  # still half the thurst becasue hasnt reacted
        else:
            T = efficiencyarray * 550.0 * P_TO / (2 * V)

        a = g / W_TO * (T - D - mu * (W_TO - L))
        # Integration
        V = V + a * dt
        s = s + V * dt + 0.5 * a * dt ** 2
        t += dt
        t_rec += dt

    # break applied
    while V > 0.0 and t < max_time:
        q = 0.5 * rho * V ** 2
        L = min(q * S * C_LTO, W_TO)
        D = q * S * (C_D0 + CDi_IGE)
        # stop it from going to infinity
        T = 0

        a = g / W_TO * (T - D - mu_break * (W_TO - L))
        # Integration
        V = V + a * dt
        s = s + V * dt + 0.5 * a * dt ** 2
        if V < 0:
            break
        t += dt
    return s


def V1calculation(d_TO, dt=0.05):
    V1_candidates = np.linspace(0.5 * V_STO, 1.15 * V_STO, 300)
    s_stop_list = []
    s_go_list = []

    for V1 in V1_candidates:
        s_stop_list.append(acceleratestop(V1, dt=dt))
        s_go_list.append(acceleratego(V1, dt=dt))

    s_stop = np.array(s_stop_list)
    s_go = np.array(s_go_list)

    # figure out where they intersect
    diff = s_go - s_stop
    idx = np.argmin(np.abs(diff))
    V1 = V1_candidates[idx]
    BFL = 0.5 * (s_stop[idx] + s_go[idx])

    # plot
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(V1_candidates / kntstofps, s_stop, label="Accelerate-stop")
    ax.plot(V1_candidates / kntstofps, s_go, label="Accelerate-go")
    ax.axhline(d_TO, linestyle="--", color="red", label=f"Field limit = {d_TO:.0f} ft")
    ax.axvline(V1 / kntstofps, linestyle="--", color="gray", label=f"V1 = {V1 / kntstofps:.1f} kts")
    ax.set_xlabel("V1 [kts]")
    ax.set_ylabel("Distance [ft]")
    ax.set_title("Balanced Field Length ")
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.show()

    return V1, BFL


def sensitivity_analysis(steep=False, temperature=False, weight=False, dt=0.05, max_time=200.0):
    altitudes = np.arange(0, 5001, 50)  # altitude above 0 [ft]

    # steepness change
    if steep:
        slopes_deg = [0, 1, 2, 3]
        fig, ax = plt.subplots(figsize=(9, 5))

        for slope_deg in slopes_deg:
            slope_rad = np.radians(slope_deg)
            s_LO_list = []
            for alt in altitudes:
                atmos_model = Atmosphere(alt, delta_T)
                rho_local = atmos_model.density[0]*kgtoslug
                s_LO = ground_run_for(rho_local, W_TO, slope_rad, dt=dt, max_time=max_time)
                s_LO_list.append(s_LO)

            ax.plot(s_LO_list, altitudes, label=f"Slope = {slope_deg}°")

        ax.axvline(d_TO, linestyle="--", color="red", label=f"Field limit = {d_TO:.0f} ft")
        ax.set_xlabel("Ground run to lift-off  [ft]")
        ax.set_ylabel("Altitude above sea level [ft]")
        ax.set_title("Sensitivity: runway slope")
        ax.grid(True)
        ax.legend()
        plt.tight_layout()
        plt.show()

    # temperature change
    if temperature:
        delta_T_list = np.arange(0, 41, 10)  # ISA deviation we need 20K
        fig, ax = plt.subplots(figsize=(9, 5))

        for dT in delta_T_list:
            s_LO_list = []
            for alt in altitudes:
                atmos_model = Atmosphere(alt, dT)
                rho_local = atmos_model.density[0]*kgtoslug
                s_LO = ground_run_for(rho_local, W_TO, 0, dt=dt, max_time=max_time)
                s_LO_list.append(s_LO)

            ax.plot(s_LO_list, altitudes, label=f"ISA + {dT:.0f} °C")

        ax.axvline(d_TO, linestyle="--", color="red", label=f"Field limit = {d_TO:.0f} ft")
        ax.set_xlabel("Ground run to lift-off  [ft]")
        ax.set_ylabel("Altitude above sea level  [ft]")
        ax.set_title("Sensitivity: temperature deviation (ISA + ΔT)")
        ax.grid(True)
        ax.legend()
        plt.tight_layout()
        plt.show()

    # weight change
    if weight:
        weight_fractions = np.linspace(0.80, 1.20, 5)
        weights_lbs = W_TO * weight_fractions
        fig, ax = plt.subplots(figsize=(9, 5))

        for w in weights_lbs:
            s_LO_list = []
            for alt in altitudes:
                atmos_model = Atmosphere(alt, delta_T)
                rho_local = atmos_model.density[0]*kgtoslug
                s_LO = ground_run_for(rho_local, w, 0, dt=dt, max_time=max_time)
                s_LO_list.append(s_LO)

            frac = w / W_TO
            ax.plot(s_LO_list, altitudes, label=f"{frac * 100:.0f}% MTOW ({w:.0f} lbs)")

        ax.axvline(d_TO, linestyle="--", color="red", label=f"Field limit = {d_TO:.0f} ft")
        ax.set_xlabel("Ground run to lift-off  [ft]")
        ax.set_ylabel("Altitude above sea level  [ft]")
        ax.set_title("Sensitivity: take-off weight (80–120% MTOW)")
        ax.grid(True)
        ax.legend(fontsize=8)
        plt.tight_layout()
        plt.show()


# set up everything
if __name__ == "__main__":
    print("Values")
    print(f"  W_TO  = {W_TO:.1f} lbs")
    print(f"  S     = {S:.0f} ft²")
    print(f"  P_TO  = {P_TO:.1f} hp")
    print(f"  V_STO = {V_STO:.2f} ft/s ")

    # Roskam
    print("Roskam Method")
    s_gr_roskam, TOP, T_roskam = Roskam_TO(C_D0, C_LmaxTO, P_bhp, S_TOG)
    print(f"    Ground run  = {s_gr_roskam:.1f} ft")
    print(f"    Avg thrust  = {T_roskam:.1f} lbs")
    print(f"    TOP         = {TOP:.2f} lbs²/(ft²·hp)")

    # Torenbeek analytical
    print("Gorenbeek Method (analytical)")
    T2, s_gr_tbk = GORENBEEK_TO_1(P_TO, b, h, c_Di, S, C_D0)
    print(f"    Ground run  = {s_gr_tbk:.1f} ft")
    print(f"    Thrust      = {T2:.1f} lbs")

    # Torenbeek numerical
    print("Gorenbeek Method (numerical integration)")
    s, V, t, traj, T_TR, V1, s_LO, V_LO = GORENBEEK_TO_3()
    print(f"    Lift-off distance = {s_LO:.1f} ft")
    print(f"    Lift-off speed    = {V_LO:.2f} ft/s  ({V_LO / kntstofps:.1f} kts)")
    print(f"    V1     = {V1:.2f} ft/s" if V1 is not None else "    V1     = not reached before lift-off")

    # Rest of the phase
    S_OBS, h_TR = rest(T_TR)
    S_TO_total = s_LO + S_OBS
    print(f"    Airborne distance = {S_OBS:.1f} ft ")
    print(f"    Taake-off distance = {S_TO_total:.1f} ft  ")

    # BFL
    V1_speed, BFL = V1calculation(d_TO)
    print(f"    V1  = {V1_speed / kntstofps:.1f} kts")
    print(f"    BFL = {BFL:.1f} ft")

    # -- Sensitivity analyses --
    sensitivity_analysis(steep=True)
    sensitivity_analysis(temperature=True)
    sensitivity_analysis(weight=True)

"""
not applicable in the end too big airplanes so other check
def BFL() : #checked assumpiton is the only check left
    S_GR1, TOP, T = Roskam_TO(C_D0, C_LmaxTO, P_bhp, S_TOG)
    T_OEI=1/2*T #check if this assumption is valid lbs
    CDi_IGE = CDi_ground_effect(h, b,c_Di)
    V2=1.2*V_STO
    D2=1/2*rho*V2**2*S*(CDi_IGE+C_D0) # increase in drag due to ground effetc

    gamma2=np.arcsin((T_OEI-D2)/W_TO) #this is found in toreenbeekk P.161
    gamma2min=0.024
    Dgamma=gamma2-gamma2min
    W_TOS=W_TO/S

    TW=T/W_TO
    mup=0.01*C_LmaxTO+0.02
    dS_to=1.1*V_STO*2 #[feet] inertia distance was fixed with torenbeek at 655 but doesn't make sense so be conservative and assume reaction time of 2 second and speed of 1.1*Vstol
    CL2=0.694*C_LmaxTO #lift coefficient at V2 assume V2=1.2Vs

    balanced=0.863/(1+2.3*Dgamma)*(W_TOS/(rho*g*CL2)+h_TO)*(2.7+1/(TW-mup))+(dS_to/np.sqrt(theta))


    return balanced

def Sensitivity_analysis(Steep=False,altitude=False,Temperature=False,Weight=False,dt=0.05, max_time=200):
    #steeper in function of speed
    if Steep:
        plt.figure()
        angles = np.radians([0, 1, 2, 3])

        for angle in angles:
            V = 0.1
            s = 0.0
            t = 0.0
            trajectory = []
            CDi_IGE = CDi_ground_effect(h, b)

            while t < max_time:

                q = 0.5 * rho * V ** 2
                L = q * S * C_LTO
                D = q * S * (C_D0 + CDi_IGE)
                liftoff = None
                V_LOF = 1.556 * np.sqrt(W_TO / (rho * S * C_LmaxTO))



                if V > 0:
                    T = efficiencyarray * (550 * P_TO) / V
                a = (g / W_TO) * (
                        T
                        - D
                        - mu * (W_TO - L)
                        - W_TO * np.sin(angle)
                )

                V_new = V + a * dt
                s_new = s + V * dt + 0.5 * a * dt **2
                trajectory.append([t, V, s, a, L, T])

                # stop condition: liftoff
                if liftoff is None and V >= V_LOF:
                    break

                V = max(V_new, 0.1)
                s = s_new
                t += dt
            traj = np.array(trajectory)
            plt.plot(traj[:, 2], traj[:, 1], label=f"θ = {np.degrees(angle):.1f}°")

        # V_LOF line
        V_LOF=traj[-1][1]
        plt.axhline(V_LOF, linestyle="--", color="black", label="V_LOF")

        plt.xlabel("Takeoff distance [ft]")
        plt.ylabel("Speed [ft/s]")
        plt.title("Sensitivity: thrust effect on takeoff acceleration curve")
        plt.grid()
        plt.legend()
        plt.show()


    #altitude sensitivity
    #Temperature sensitivity
    #weight sensitivity
    return

"""

# set up everything

"""
print("Numerical TO distance:", S_TO_num)
print("Final speed:", V_final)
print("Time:", t_final)
"""
"""
plt.figure()
plt.plot(traj[:,1], traj[:,5])
#plt.plot(traj[:,0], traj[:,-1])
plt.xlabel("Airspeed [ft/s]")
plt.ylabel("Thrust [lbs]")
plt.title("Thurst vs Airspeed")
plt.grid()

#airspeed vs groyudn run
plt.figure()
plt.plot(traj[:,2], traj[:,1])

plt.xlabel("Distance [ft]")
plt.ylabel("Airspeed [ft/s]")
plt.title("Airspeed Vs Ground run")
plt.grid()
plt.show()
a=Sensitivity_analysis(Steep=False,altitude=False,Temperature=False,Weight=False,dt=0.05, max_time=200)


maybe useful plots
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].plot(traj[:, 2], traj[:, 1])
    axes[0].axvline(s_LO, linestyle="--", color="red", label=f"Lift-off @ {s_LO:.0f} ft")
    axes[0].set_xlabel("Ground run [ft]")
    axes[0].set_ylabel("Speed [ft/s]")
    axes[0].set_title("Speed vs ground run (numerical)")
    axes[0].grid(True)
    axes[0].legend()

    axes[1].plot(traj[:, 1], traj[:, 5])
    axes[1].set_xlabel("Speed [ft/s]")
    axes[1].set_ylabel("Thrust [lbs]")
    axes[1].set_title("Thrust vs speed (numerical)")
    axes[1].grid(True)

    plt.tight_layout()
    plt.show()
"""

# created a graph of wing power per wing loading
# higher is bad than line
"""
CLmax_array = np.arange(1.0, 2.5, 0.2)
WS = np.linspace(5, 60, 500)

plt.figure(figsize=(10,6))
for CLmax in CLmax_array:
    WP = TOP * theta * CLmax / WS
    plt.plot(WS, WP, label=f"CLmax = {CLmax:.1f}")

plt.xlabel("Wing Loading W/S [lb/ft²]")
plt.ylabel("Power Loading W/P [lb/hp]")
plt.title("Takeoff Constraint")

plt.grid(True, which='major', linestyle='-')
plt.grid(True, which='minor', linestyle=':', alpha=0.5)
plt.minorticks_on()

plt.legend(ncol=2)
plt.xlim(0,60)
plt.show()
"""