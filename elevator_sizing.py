import matplotlib
matplotlib.use("TkAgg")

import numpy as np
import matplotlib.pyplot as plt
from classes.isa import Atmosphere
import scipy.interpolate as intp

atm = Atmosphere(609.6, 20)
atm_c = Atmosphere(2590.8)
# Parameters
rho = atm.density[0]               # at 2000 feet, ISA + 20 conditions
rho_c = atm_c.density[0]           # at 8500 feet

angular_acc = 8 * np.pi/180               # rad/s^2, assumption
Iyy_mg = 8745.44

Vs = 50*0.51444444   # m/s
Vs_c = 28.8
Vc = 130*0.51444444  # m/s
Vmax = 150*0.5144444            # m/s

Wmto = 1840              # kg
T_to = 300000 * 0.8 / (1.2*Vs)   # N, 2*Tmax

Sw = 31.4            # m^2
cw_mac = 1.87        # m
AR = 9
e = 0.783
Hfus = 1.7           # m, fuselage height
Hlg = 0.8            # m, height of landing gear (ground to bottom of fuselage)

alpha = 0                          # aoa at the moment of rotation, assumption
ih = -2.5 * np.pi/180
iw = 0 * np.pi/180
alpha_w = iw                            # assumption

CL_alpha = 5.35                          # For clean wing
#Cm_alpha1 = ...                         # most fwd cg
Cm_alpha2 = -1                           # most aft cg, OpenVSP
CL_to = Wmto*9.81/(0.5*rho*Vs**2*Sw)
CD_0 = 0.03144
CD_to = CD_0 + CL_to**2/(np.pi*AR*e)
Cm_ac_wf = -0.01

# Locations (from the nose)
x_cg_fwd = 4.9
x_cg_aft = 5.4
x_mg = 5.56
z_mg = 0
zcg = 0.4 * Hfus + Hlg
z_D = Hfus + Hlg
z_T = 1.7 - 0.2 + 0.8
x_ac_wf = 4.4 + 1.87 * 0.25
x_ac_h = 10.4

Sh = 6.7
Ah = 6
de_max = -25 * np.pi/180           # degrees, assumption
CL_alpha_h = 4.13                   # it exists
vol_coef_h1 = Sh * (x_ac_h - x_cg_fwd)/(Sw * cw_mac)   # most fwd c.g.
vol_coef_h2 = Sh * (x_ac_h - x_cg_aft)/(Sw * cw_mac)   # most aft c.g.
b_h = np.sqrt(Sh * Ah)
b_e = 0.9*b_h                          # assumption that elevator span is = horizontal tail span
alpha_h_stall_noE = 9              # depends on the airfoil NACA0014
delta_alpha_hE = 8.1               # assumption from Table 12.19


Cm0 = Cm_ac_wf - CL_alpha_h * ih * vol_coef_h2

friction = 0.07                   # assumption from Table 9.7


# TAKE_OFF ROLL REQUIREMENT
D_to = 0.5*rho*Vs**2*Sw*CD_to
L_wf = 0.5*rho*Vs**2*Sw*CL_to
Mac_wf = 0.5*rho*Vs**2*Sw*Cm_ac_wf*cw_mac
ma = T_to - D_to - friction*(Wmto*9.81 - L_wf)

num_lh = L_wf*(x_mg - x_ac_wf) + Mac_wf + ma*(zcg - z_mg) - Wmto*9.81*(x_mg - x_cg_fwd) + D_to*(z_D - z_mg) - T_to*(z_T - z_mg) - Iyy_mg*angular_acc
den_lh = x_ac_h - x_mg
Lh = num_lh / den_lh
print(f'force by tail is {Lh}')

CLh = Lh / (0.5*rho*Vs**2*Sh)
print(f'CLh is {CLh}')

epsilon0 = (2*CL_to)/(np.pi * AR)
de_da = 0
downwash = 0
alpha_h = alpha + ih - downwash


tau_e = (CLh - CL_alpha_h*alpha_h) / (CL_alpha_h*de_max)
print(f'tau_e = {tau_e}')
tau_e = 0.60

# LONGITUDINAL TRIM REQUIREMENTS
CL_de = CL_alpha_h * Sh/Sw * b_e/b_h * tau_e

print(f'CL_de is: {CL_de}')
V_list = np.linspace(Vs_c, Vmax, 50)
de_list = []

# LONGITUDINAL TRIM REQUIREMENT (most aft c.g., sea level)
Cm_de2 = - CL_alpha_h * vol_coef_h2 * b_e/b_h * tau_e

for V in V_list:
    T = 159000/V
    CL_1 = 0.95*Wmto*9.81/(0.5*rho*V**2*Sw)
    CL_0 = 0
    de_num = ((T * (-z_T))/(0.5*rho*V**2*Sw*cw_mac) + Cm0) * CL_alpha + (CL_1 - CL_0) * Cm_alpha2
    de_den = CL_alpha * Cm_de2 - Cm_alpha2 * CL_de
    de = - de_num/de_den
    de_list.append(de * 180/np.pi)



# LONGITUDINAL TRIM REQUIREMENT (most aft c.g., cruise)
de_list_2 = []
for V in V_list:
    T = 159000/V
    CL_1 = 0.95*Wmto*9.81/(0.5*rho_c*V**2*Sw)
    CL_0 = 0
    de_num = ((T * (-z_T))/(0.5*rho_c*V**2*Sw*cw_mac) + Cm0) * CL_alpha + (CL_1 - CL_0) * Cm_alpha2
    de_den = CL_alpha * Cm_de2 - Cm_alpha2 * CL_de
    de = - de_num/de_den
    de_list_2.append(de * 180/np.pi)

print(f'de_max_up with most aft c.g. [degrees] = {de_list_2[0]}')
print(f'de_max_down with most aft c.g. [degrees] = {de_list_2[-1]}')

#plt.plot(V_list, de_list_c1, marker='o', label='With most fwd c.g.')
# plt.plot(V_list, de_list_c2, marker='s', label='With most aft c.g.')
# plt.title('Elevator deflection at cruise conditions')
# plt.legend()
# plt.show()

# STALL CHECK
alpha_h_stall = (alpha_h_stall_noE - delta_alpha_hE)
alpha_h_to = alpha_w + ih
if abs(alpha_h_to) < abs(alpha_h_stall):
    print('Stall Check Pass')
else:
    print('Stall Check Fail')
    print(f'alpha_h_stall = {alpha_h_stall*180/np.pi}')
    print(f'alpha_h_to = {alpha_h_to*180/np.pi}')

V_func = intp.CubicSpline(de_list_2, V_list)
print(f"trim speed is : {V_func(0) / 0.51444}")
print(f'speed for max deflection is [m/s]: {V_func(-25)}')

fig, ax = plt.subplots(figsize=(8, 5))

# Sea level
ax.plot(
    V_list,
    de_list,
    '-o',
    linewidth=2,
    markersize=4,
    label='Sea level'
)

# Cruise
ax.plot(
   V_list,
    de_list_2,
    '-s',
    linewidth=2,
   markersize=4,
   label='Cruise'
)

# Reference speeds
ax.axvline(x=Vs_c, color='blue', linestyle='--', linewidth=1, label=r'$V_{stall}$')
ax.axvline(x=Vc, color='red', linestyle='--', linewidth=1, label=r'$V_{cruise}$')

ax.set_xlabel('Airspeed [m/s]')
ax.set_ylabel(r'Elevator Deflection $\delta_e$ [deg]')
ax.set_title('Elevator Trim Requirement')

ax.grid(True, linestyle='--', alpha=0.7)
ax.legend()

plt.tight_layout()
plt.savefig("elevator_trim.png", dpi=300)
plt.show()

