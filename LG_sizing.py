import numpy as np

# Misc
Mmto = 1840
Wmto = Mmto * 9.81                          # N
Wland = 0.95 * Mmto                         # N
Vs = 50 * 0.51444444                        # m/s
f_nlg_min = 0.08
f_nlg_max = 0.15
hfus = 1.7                                  # m

# Main Landing Gear
dmlg = 13.25 * 0.0254                       # m
wmlg = 5.05 * 0.0254                         # m
static_loaded_mlg = 5.2 * 0.0254          # m
longitudinal_mlg = 6.01                     # m, from the nose
min_height_mlg = 0.652                      # m
Rated_load_mlg = 9900 * 0.45359237 * 9.81   # N

# Nose Landing Gear
dnlg = 13.25 * 0.0254                       # m
wnlg = 5.05 * 0.0254                        # m
static_loaded_nlg = 5.2 * 0.0254            # m
longitudinal_nlg = 0.877                     # m, from the nose
min_height_nlg = 0.33655                    # m
Rated_load_nlg = 1200 * 0.45359237 * 9.81   # N

# Assumptions
N_gear = 3                                  # Gear load factor for GA aircraft (T11.5 Raymer)
V_vertical = 3                              # Sink speed [m/s] (Raymer)
shock_eff = 0.65                            # Shock absorber efficiency, conservative estimate for fixed oleo absorber (T11.4 Raymer)
tire_eff = 0.47                             # Tire efficiency (T11.4 Raymer)
P = 1800 * 6894.75729                       # Oleo internal pressure [Pa] (Raymer)

'''---SHOCK ABSORBER SIZING--------------------------------------------------------------------------------'''

# Shock absorber sizing
St_mlg = 0.5 * dmlg - static_loaded_mlg
S_mlg = V_vertical**2/(2*9.81*shock_eff*N_gear) - tire_eff/shock_eff * St_mlg + 1 * 0.0254
length_oleo_mlg = 2.5 * S_mlg

S_nlg = S_mlg
length_oleo_nlg = 2.5 * S_nlg

Loleo_mlg = 1.1 * 0.5 * Wmto * (1 - f_nlg_min)   # 110% MTOW
Doleo_mlg = 1.3 * np.sqrt(4*Loleo_mlg/(P*np.pi))

Loleo_nlg = 1.1 * Wmto * f_nlg_max
Doleo_nlg = 1.3 * np.sqrt(4*Loleo_nlg/(P*np.pi))

# Wheel swelling
extra_d_mlg = 0.03 * dmlg
extra_d_nlg = 0.03 * dnlg

mlg_height_total = length_oleo_nlg + dmlg + extra_d_mlg
nlg_height_total = length_oleo_nlg + dnlg + extra_d_nlg

# Dynamic Load check (For braking on concrete)
hcg = length_oleo_nlg + dmlg + extra_d_mlg + 0.5 * hfus
F_dyn = 3 * hcg * Wmto / (9.81 * (longitudinal_mlg - longitudinal_nlg))
F_static = F_dyn/1.45
if F_static < Rated_load_nlg:
    print('Dynamic Load check PASSED')

print('SHOCK ABSORBER')
print('Main Landing Gear')
print(f'Tire diameter: {dmlg}')
print(f'Oleo Length: {length_oleo_mlg}')
print(f'Oleo Diameter: {Doleo_mlg}')
print(f'Height MLG: {mlg_height_total}')

print('Nose Landing Gear')
print(f'Tire diameter: {dnlg}')
print(f'Oleo Length: {length_oleo_nlg}')
print(f'Oleo Diameter: {Doleo_nlg}')
print(f'Height NLG: {nlg_height_total}')

'''-----BRAKE SIZING----------------------------------------------------------------------------------'''
# Brake requirements
V_land = Vs
V_rejected = 1.2 * Vs

KE = 0.5 * V_rejected**2 * Mmto
m = (0.49 * (1000-50))/KE
m_alt = KE/(1000*750)
print(f'm_brake: {m}')
print(f'm_alt: {m_alt}')

# Vb = 1.2 * Vs * 3.2808399      # ft/s
# W1 = 0.95 * Mmto * 2.20462262  # Design Landing Weight [lbs], 100 stops at 10ft/s^2
# W2 = 0.99 * Mmto * 2.20462262  # Maximum Landing Weight [lbs], 5 stops at 10ft/s^2
# W3 = Mmto * 2.20462262         # Rejected TO [lbs], 1 stop at 6ft/s^2
#
#
# KEbrake1 = 0.0443*W1*Vb**2 / 10e6
# KEbrake2 = 0.0443*W2*Vb**2 / 10e6
# KEbrake3 = 0.0443*W3*Vb**2 / 10e6
#
# Wbrake1 = -1.12e-1 * KEbrake1**2 + 16.7 * KEbrake1 + 13.7
# Wbrake3 = -9.90e-3 * KEbrake3**2 + 5.41 * KEbrake3 + 9.97e-1
# Wbrake2 = -2.99e-2 * KEbrake2**2 + 8.46 * KEbrake2 - 2.10
# Wbrake = np.average([Wbrake1, Wbrake2, Wbrake3])
# # Design landing weight, 250 stops
# Material = input("Enter brake material (Steel or Carbon): ").strip().capitalize()
#
# print('BRAKES')
# if Material == 'Steel':
#     V = 2.55 * Wbrake                 # Volume eq. derived from
#     V2 =  y = -5.026722953e-11*(Wbrake**5)+1.062448154e-7*(Wbrake**4)-0.00004980130665*(Wbrake**3)+0.009869164127*(Wbrake**2)+2.089188985*Wbrake+4.366837109
#     print(f'Material: {Material}')
#     print(f'Weight [kg]: {Wbrake * 0.45359237}')
#     print(f'Volume [m^3]: {V * 0.00004916}')
#     print(f'V2: {V2 * 0.00004916}')
#
# elif Material == 'Carbon':
#     Wbrake = 2.5 * Wbrake
#     V = (3.3 * Wbrake - 84.2) * 1.28
#     print(f'Material: {Material}')
#     print(f'Weight [kg]: {Wbrake * 0.45359237}')
#     print(f'Volume [m^3]: {V * 0.00004916}')
#
# else:
#     print("Invalid material. Please enter 'Steel' or 'Carbon'.")

'''----PRELIM WEIGHT--------------------------------------------------------------'''
# From S. Currey
Klg = 0.85 + 0.15 + 0.11
Wlg1 = 0.046 * Klg * Wland
Wlg2 = 20.45 * 0.7 * (Wland * 2.20462262 * 1/1000)**1.17 * 0.45359237


# From Raymer (for GA)
Wmlg = 0.095 * (N_gear * 1.5 * Wland/9.81 * 1/0.45359237)**0.768 * (mlg_height_total *39.3700787/12)**0.409 /0.45359237
Wnlg = 0.125 * (N_gear * 1.5 * Wland/9.81 *1/0.45359237)**0.566 * (nlg_height_total *39.3700787/12)**0.845 /0.45359237
print(f'Raymer Weight estimate: {Wmlg + Wnlg}')


print(f'Wlg1: {Wlg1}')
print(f'Wlg2: {Wlg2}')