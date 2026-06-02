from classes.aircraft_2 import Aircraft
from c2_m import Snet, LE_sweep_deg, sweep_at_x_c_deg, lift_slope, S_wf, closest_value, chord_at_y_span
from lookups.consts import *
import pandas as pd
import numpy as np
from scipy.interpolate import RegularGridInterpolator, interp1d
from scipy.optimize import brentq
from classes.isa import Atmosphere
import matplotlib.pyplot as plt

''' How to reduce drag:

    Nacelle-wing interference:      Smaller local wing chord, 
                                    less wide nacelle, 
                                    engine below wing, 
                                    see optimal incidence and position in eqn 4.63 Roskam VI
    Nose landing gear:              make x_nlg close to 2x tire diameter
    Nacelle drag:                   Smaller base area,
                                    larger length/diameter ratio
    
    '''

# Data and interpolation functions
def read_csv(path):
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()
    return df

def interp_value(df : pd.DataFrame,
                 x_query,
                 x_col : str,
                 y_col : str,
                 log_x = False) -> float:
    data = df[[x_col, y_col]].dropna().sort_values(x_col)
 
    x = data[x_col].to_numpy(dtype = float)
    y = data[y_col].to_numpy(dtype = float)
 
    if log_x:
        x = np.log10(x)
        x_query = np.log10(x_query)
       
 
    return float(np.interp(x_query, x, y))

# Intermediate variable functions
def mu_air(T):
    return 1.81 * 1e-5 * (T / 293.15)**1.5 * (293.15 + 110.4) / (T + 110.4)

def R_wf(M, R_N_fus):
    M = closest_value(M, values=[0.25, 0.4])
    if M == 0.25:
        return interp_value(read_csv('lookups/roskam_p6_fig_4_1_rwf.csv'), R_N_fus, x_col='Fuselage Reynolds Number (M = 0.25)', y_col='Wing-Fuse Interference Factor (M = 0.25)', log_x=True)
    else:
        return interp_value(read_csv('lookups/roskam_p6_fig_4_1_rwf_2.csv'), R_N_fus, x_col='Fuselage Reynolds Number (M = 0.4)', y_col='Wing-Fuse Interference Factor (M = 0.4)', log_x=True)

def R_LS(sweep_t_c_max_deg):
    return interp_value(read_csv('lookups/roskam_p6_fig_4_2_rls.csv'), np.cos(np.deg2rad(sweep_t_c_max_deg)), 'cos(quarter chord)', 'lifting surface correction', log_x=False)

def R_N(rho, V, l, mu):
    return rho * V * l / mu

def C_f(R_N, M):  # DATCOM fig 4.1.5.1-26
    x = np.log10(R_N)
    # CF3 = 3.92725e-6 * x**5 -1.30370e-4 * x**4 -1.65388e-3 * x**3 -9.59519e-3 * x**2 +2.18366e-2 * x 
    # CF0 = 4.12963e-6 * x**5 -1.36204e-4 * x**4 + 1.71620e-3 * x**3 -9.88935e-3 * x**2 +2.23641e-2 * x
    if M > 0.15:
        CF3 = x * (2.18366e-2 + x * (-9.59519e-3 + x * (1.65388e-3 + x * (-1.30370e-4 + x * 3.92725e-6))))
        # print(f' R_N: {R_N}, CF3: {CF3}')
        return CF3
    else: 
        CF0 = x * (2.23641e-2 + x * (-9.88935e-3 + x * (1.7162e-3 + x * (-1.36204e-4 + x * 4.12963e-6))))
        return CF0

def L_dash(x_c_t_c_max):
    if x_c_t_c_max < 0.3:
        return 2.0
    else:
        return 1.2

def S_wet_w(S_exp, surface_type: str = 'wing' # 'wing' or other like 'empennage' or 'pylon'
            ):
    if surface_type == 'wing':
        return 2 * 1.07 * S_exp
    else:
        return 2 * 1.05 * S_exp
    
def S_wet_fus(l_nosecone, l_tot, l_tailcone, d_max):
    l2 = l_tot - l_nosecone - l_tailcone
    return np.pi * d_max / 4 * (1 / (3 * l_nosecone**2) * ((4 * l_nosecone**2 + d_max**2 / 4)**1.5 - d_max**3 / 8) - d_max + 4 * l2 + 2 * np.sqrt(l_tailcone**2 + d_max**2 / 4))

def Snet(S, b_f, taper, b, c_r):
        c_fus_int = chord_at_y_span(c_r, taper, b_f/2, b)
        return S - (c_r + c_fus_int) * b_f / 2

def exposed_wing_area_and_mgc(S, b, taper, c_r, b_f):
    c_fus_int = chord_at_y_span(c_r, taper, b_f/2, b)
    S_net = S - (c_r + c_fus_int) * b_f / 2
    return S_net, S_net / (b - b_f)

def exposed_vt_area_and_mgc(S_fus_base, l_f, l_tc, x_le_vt, S_vt, b_vt, taper_vt, c_r_vt, h_f):
    d_min = np.sqrt(4 * S_fus_base / np.pi)
    d_fus_local = h_f - (x_le_vt - (l_f - l_tc)) / l_tc * (h_f - d_min)
    c_fus_int = chord_at_y_span(c_r_vt, taper_vt, d_fus_local/2, b_vt)
    S_net = S_vt - (c_r_vt + c_fus_int) / 2 * d_fus_local / 2
    return S_net, S_net / (b_vt - d_fus_local / 2)


def D_CD_flap_stuff(ac: Aircraft, flap_type, t_c_max, flap_deflection, wing_drag, flight_condition, cdash_c):
    cf_c = ac.hld_and_ailerons.flaps['cf_c']  # flap chord length / chord length
    if flap_type == 'split':  # 'split' or 'plain' or 'slotted' 'fowler' or 'krueger'
        fd = closest_value(t_c_max, values=[10, 20, 30])
        D_CD_flap_stuff = interp_value(read_csv(f'lookups/t_c_0.{fd}.csv'), cf_c, f'cf/c ({flight_condition})', f'dCdp ({flight_condition})', log_x=False)
    elif flap_type == 'plain':
        fd = closest_value(flap_deflection, values=[15, 60])
        D_CD_flap_stuff = interp_value(read_csv(f'lookups/d_f_{fd}.csv'), cf_c, 'cf/c', 'dCdp', log_x=False)
    elif flap_type == 'slotted':
        fd = closest_value(cf_c, values=[0.1, 0.2, 0.3])
        D_CD_flap_stuff = interp_value(read_csv('lookups/cf_c_comb2.csv'), flap_deflection, 'df', f'dCdp (cf={fd}0)', log_x=False)
    elif flap_type in ['fowler', 'double slotted', 'triple slotted']:
        fd = f"{closest_value(cf_c, values=[0.1, 0.2, 0.3, 0.4]):.1f}"
        # fd = closest_value(cf_c, values=[0.1, 0.2, 0.3, 0.4])
        df = read_csv('lookups/roskam_p6_fig_4_48.csv')
        print(df.columns.tolist())
        print(f'df(cf/c={fd})', f'dCdp(cf/c={fd})')
        D_CD_flap_stuff = interp_value(read_csv('lookups/roskam_p6_fig_4_48.csv'), flap_deflection, f'df(cf/c={fd})', f'dCdp(cf/c={fd})', log_x=False)
        # D_CD_flap_stuff = interp_value(read_csv('lookups/roskam_p6_fig_4_48.csv'), flap_deflection, f'df({fd})', f'dCdp(cf/c={fd})', log_x=False)
    elif flap_type == 'kruger':
        D_CD_flap_stuff = wing_drag * cdash_c   # (cf_c * np.cos(np.deg2rad(flap_deflection)) + 1)
    elif flap_type in ['slat', 'fixed slot', 'leading edge flap']:
        D_CD_flap_stuff = wing_drag * cdash_c
    else: 
        raise ValueError(f"Flap type give: {flap_type}, check possible entries for D_CD_flap_stuff variable")
    return D_CD_flap_stuff

def flap_interference_factor(flap_type: str, # 'split' or 'plain' or 'slotted' 'fowler' or 'krueger'
                             ):
    if flap_type == 'split':
        return -0.15
    elif flap_type == 'plain':
        return 0
    elif flap_type == 'slotted':
        return 0.4
    elif flap_type == 'fowler':
        return 0.25
    else:
        return 0.1
# Component drag eqns

def CD0_wing(M, R_N_fus, sweep_t_c_max_deg, R_N_w, t_c_max, x_c_t_c_max, S_exp, surface_type, S):
    if surface_type == 'wing':
        Rwf = R_wf(M, R_N_fus)
    else:
        Rwf = 1
    RLS = R_LS(sweep_t_c_max_deg)
    Cfw = C_f(R_N_w, M)
    Ldash = L_dash(x_c_t_c_max)
    Swet = S_wet_w(S_exp, surface_type)
    # print(f' \n Wing stuff: \n Rwf: {Rwf}, RLS: {RLS}, Cfw: {Cfw}, Ldash: {Ldash}, Swet: {Swet}')
    return Rwf * RLS * Cfw * (1 + Ldash * t_c_max + 100 * t_c_max**4) * Swet / S

def CD0_fuselage(M, R_N_fus, l_fus, d_fus_max, l_nosecone, l_tailcone, S, d_b, S_fus_max, surface_type:str = 'fuselage'):
    Rwf = R_wf(M, R_N_fus)
    Cf = C_f(R_N_fus, M)
    if surface_type == 'fuselage':
        S_wet = S_wet_fus(l_nosecone, l_fus, l_tailcone, d_fus_max)
    elif surface_type == 'nacelle':
        S_wet = np.pi * (d_fus_max * l_fus + d_b**2 / 4)
    CD0_f_less_b = Rwf * Cf * (1 + 60 / ((l_fus / d_fus_max)**3) + 0.0025 * (l_fus / d_fus_max)) * S_wet / S
    # print(f'CD0_f_less_b: {CD0_f_less_b} \n Rwf: {Rwf}, Cf: {Cf}, S_wet: {S_wet}')
    CD0_b = (0.029 * (d_b / d_fus_max)**3 / (CD0_f_less_b * (S / S_fus_max))**0.5) * S_fus_max / S
    return CD0_f_less_b + CD0_b

def CD0_wing_nacelle_interference(i_n_deg, nacelle_above_wing: bool, chord_at_nacelle, b_n, S):
    if nacelle_above_wing:
        DCl1 = 0.2
    else:
        DCl1 = -0.3
    DCl2 = 0.056 * i_n_deg
    return 0.036 * (chord_at_nacelle * b_n / S) * (DCl1 + DCl2)**2

def CD0_windmilling_prop(V_fps, q_psf, S_sqft, SHP_max):
    return 33 / (q_psf * S_sqft) * SHP_max / V_fps

def CD0_windmilling_prop_inop(n_blades, D_prop, S):
    return 0.00125 * n_blades * D_prop**2 / S

def CD0_flap_profile(ac: Aircraft, flap_type, t_c_max, flap_deflection, wing_drag, flight_condition, slat_or_flap: str, cdash_c):
    D_CD_flap = D_CD_flap_stuff(ac, flap_type, t_c_max, flap_deflection, wing_drag, flight_condition, cdash_c)
    if slat_or_flap == 'flap':
        Swf = ac.hld_and_ailerons.flaps['S_wf']
    elif slat_or_flap == 'slat':
        Swf = ac.hld_and_ailerons.slats['S_wf']
    else:
        raise ValueError(f" slat_or_flap was entered as {slat_or_flap}, expected either 'slat' or 'flap' as entry")
    return D_CD_flap * np.cos(np.deg2rad(ac.wing.sweep)) * Swf / ac.wing.area

def CD0_flap_interference(CD0flap_profile, flap_type):
    return CD0flap_profile * flap_interference_factor(flap_type)

def CD0_gear(x_nlg, height_nlg, n_main_wheels, n_nose_wheels, nlg_tire_diameter, nlg_tire_width, mlg_tire_diameter, mlg_tire_width, S, fairing_type: str = 'C'):
    # Type 1 landing gear (fig 4.54 Roskam VI)
    CD = 0.71
    if fairing_type == 'A':
        CD = 1.15
    elif fairing_type == 'B':
        CD = 1.05
    CD_mlg = CD * mlg_tire_diameter * mlg_tire_width / S
    a_Dt = x_nlg / nlg_tire_diameter
    e_Dt = x_nlg / height_nlg
    a_Dt_data = closest_value(a_Dt, values=[0.75, 2.15, 3.6])
    if a_Dt_data == 0.75:
        x = [1.1, 1.5, 2, 2.4]
        y = [0.28, 0.35, 0.47, 0.55]
    elif a_Dt_data == 2.15:
        x = [0.5, 0.75, 1.1, 1.5, 2, 2.5, 2.75]
        y = [0.1, 0.25, 0.49, 0.375, 0.4, 0.49, 0.52]
    else:
        x = [1.1, 1.5, 2]
        y = [0.4, 0.48, 0.63]
    CD_nlg = float(np.interp(e_Dt, x, y)) * nlg_tire_diameter * nlg_tire_width / S
    return n_main_wheels * CD_mlg + n_nose_wheels * CD_nlg

def CD0(ac: Aircraft, 
        n_engine_inoperative: int, 
        flight_condition: str = 'cruise', # 'cruise' or 'landing' or 'take-off'
        update_ac: bool = False):
    S = ac.wing.area
    if flight_condition == 'cruise':
        temp_shift = 0
        alt = ac.requirements.cruise['cr_altitude'] * FT_TO_M
        speed = ac.requirements.cruise['cr_speed'] * KTS_TO_MS
    elif flight_condition == 'take-off':
        C_L = ac.hld_and_ailerons.take_off_lift['CL_max'] / 1.21
        temp_shift = ac.requirements.take_off['to_temp_shift']
        alt = ac.requirements.take_off['to_altitude'] * FT_TO_M
        mass_frac = ac.requirements.take_off['to_mass_frac']
        Atm = Atmosphere(alt, temp_shift)
        temp = float(Atm.temp)
        density = float(Atm.density)
        speed = np.sqrt(mass_frac * ac.weights.m_takeoff * 9.81 / (0.5 * density * C_L * S))
        flap_deflection = ac.hld_and_ailerons.flaps['to_deflection']
        cdash_c_flap = ac.hld_and_ailerons.flaps['cdash_c_to']
        cdash_c_slat = ac.hld_and_ailerons.slats['cdash_c_to']
    elif flight_condition == 'landing':
        C_L = ac.hld_and_ailerons.landing_lift['CL_max']
        temp_shift = ac.requirements.landing['la_temp_shift']
        alt = ac.requirements.landing['la_altitude'] * FT_TO_M
        Atm = Atmosphere(alt, temp_shift)
        temp = float(Atm.temp)
        density = float(Atm.density)
        mass_frac = ac.requirements.landing['la_mass_frac']
        # speed = 1.3 * ac.requirements.general['stall_speed'] * KTS_TO_MS
        speed = 1.3 * np.sqrt(mass_frac * ac.weights.m_takeoff * 9.81 / (0.5 * density * C_L * S))
        flap_deflection = ac.hld_and_ailerons.flaps['ld_deflection']
        cdash_c_flap = ac.hld_and_ailerons.flaps['cdash_c_ld']
        cdash_c_slat = ac.hld_and_ailerons.slats['cdash_c_ld']
    else:
        raise ValueError(f' Flight condition given: {flight_condition}, not "cruise" or "landing" or "take-off" ')
    
    # Atmosphere
    Atm = Atmosphere(alt, temp_shift)
    temp = float(Atm.temp)
    density = float(Atm.density)
    mu = float(mu_air(temp))
    M = speed / np.sqrt(1.4 * 287 * temp)

    # Define ac params
    f = ac.fuselage
    w = ac.wing
    ht = ac.empennage.horizontal_tail
    vt = ac.empennage.vertical_tail
    eng = ac.engine
    flap = ac.hld_and_ailerons.flaps
    slat = ac.hld_and_ailerons.slats
    lg = ac.landing_gear

    # Geometry
    b = ac.wing.span
    LE_sweep_deg = ac.wing.sweep_LE_deg
    taper = ac.wing.taper_ratio
    c_r = ac.wing.c_root
    b_f = ac.fuselage.width
    h_f = ac.fuselage.height
    l_f = ac.fuselage.length
    d_fus_max = max(h_f, b_f)
    S_base = ac.fuselage.base_area
    d_base = np.sqrt(S_base * 4 / np.pi)

    # Reynolds numbers
    R_N_fuselage = R_N(density, speed, l_f, mu)

    # Wing
    wing = CD0_wing(M, R_N_fuselage, 
                    sweep_t_c_max_deg = sweep_at_x_c_deg(LE_sweep_deg, c_r, b, taper, w.x_c_t_c_max), 
                    R_N_w = R_N(density, speed, exposed_wing_area_and_mgc(S, b, taper, c_r, b_f)[1], mu), 
                    t_c_max = w.t_c_max, 
                    x_c_t_c_max = w.x_c_t_c_max, 
                    S_exp=exposed_wing_area_and_mgc(S, b, taper, c_r, b_f)[0], 
                    surface_type='wing', S=S)

    # Fuselage
    fuselage = CD0_fuselage(M, R_N_fuselage, l_f, d_fus_max, f.nose_cone_length, f.tail_cone_length, S, d_base, f.max_cross_section_area, surface_type='fuselage')

    # HT
    h_tail = CD0_wing(M, R_N_fuselage, 
                      sweep_t_c_max_deg = sweep_at_x_c_deg(ht['sweep_LE_deg'], ht['c_r_h'], ht['b_h'], ht['taper_ratio'], ht['loc_t_c_max']), 
                      R_N_w = R_N(density, speed, exposed_vt_area_and_mgc(S_base, l_f, f.tail_cone_length, ht['x_le'], ht['area']/2, ht['b_h']/2, ht['taper_ratio'], ht['c_r_h'], h_f)[1], mu), 
                      t_c_max = ht['t_c_max'], 
                      x_c_t_c_max = ht['loc_t_c_max'], 
                      S_exp = exposed_vt_area_and_mgc(S_base, l_f, f.tail_cone_length, ht['x_le'], ht['area']/2, ht['b_h']/2, ht['taper_ratio'], ht['c_r_h'], h_f)[0], 
                      surface_type = 'ht', 
                      S = ht['area'])
    
    # VT
    v_tail = CD0_wing(M, R_N_fuselage, 
                      sweep_t_c_max_deg = sweep_at_x_c_deg(vt['sweep_LE_deg'], vt['c_r_v'], vt['b_v'], vt['taper_ratio'], vt['loc_t_c_max']), 
                      R_N_w = R_N(density, speed, exposed_vt_area_and_mgc(S_base, l_f, f.tail_cone_length, vt['x_le'], vt['area'], vt['b_v'], vt['taper_ratio'], vt['c_r_v'], h_f)[1], mu), 
                      t_c_max = vt['t_c_max'], 
                      x_c_t_c_max = vt['loc_t_c_max'], 
                      S_exp = exposed_vt_area_and_mgc(S_base, l_f, f.tail_cone_length, vt['x_le'], vt['area'], vt['b_v'], vt['taper_ratio'], vt['c_r_v'], h_f)[0], 
                      surface_type = 'vt', 
                      S = vt['area'])
    
    # Nacelle
    nacelle_isolated = CD0_fuselage(M, R_N_fus=R_N(density, speed, eng.length_nac, mu),
                                    l_fus=eng.length_nac,
                                    d_fus_max=eng.nac_diameter,
                                    l_nosecone=None,
                                    l_tailcone=None,
                                    S=S,
                                    d_b=0.5*eng.nac_diameter,
                                    S_fus_max=eng.nac_diameter**2 / 4 * np.pi,
                                    surface_type='nacelle')
    
    nacelle_interference = CD0_wing_nacelle_interference(i_n_deg=eng.i_n, 
                                                         nacelle_above_wing=eng.eng_above_wing,
                                                         chord_at_nacelle=chord_at_y_span(c_r, taper, y=eng.eng_y_pos_fuselage+b_f/2, b=b),
                                                         b_n=eng.nac_diameter,
                                                         S=S)
    
    # Propeller
    if n_engine_inoperative != 0:
        wind_milling_inop = CD0_windmilling_prop_inop(eng.n_prop_blades, eng.prop_diameter, S)
        # wind_milling = CD0_windmilling_prop(V_fps=speed*MpS_TO_FpS,
        #                                     q_psf=0.5*density*speed**2/PSF_TO_PA,
        #                                     S_sqft=S*M2_TO_F2,
        #                                     SHP_max=eng.SHP_max)
        wind_milling = 0
    else:
        wind_milling_inop = 0
        wind_milling = 0
    
    propeller = wind_milling * eng.count + wind_milling_inop * n_engine_inoperative
    
    # Flaps and Slats
    if flight_condition != 'cruise':
        flap_profile = CD0_flap_profile(ac, flap_type=flap['flap_type'], 
                                        t_c_max=w.t_c_max, 
                                        flap_deflection=flap_deflection, 
                                        wing_drag=wing, 
                                        flight_condition=flight_condition, 
                                        slat_or_flap='flap', 
                                        cdash_c=cdash_c_flap)
        flap_interference = CD0_flap_interference(flap_profile, flap['flap_type'])
        flap_drag = flap_profile + flap_interference

        slat_profile = CD0_flap_profile(ac, flap_type=slat['slat_type'],
                                        t_c_max=w.t_c_max,
                                        flap_deflection=None, 
                                        wing_drag=wing,
                                        flight_condition=flight_condition,
                                        slat_or_flap='slat',
                                        cdash_c=cdash_c_slat)
        slat_interference = CD0_flap_interference(slat_profile, slat['slat_type'])
        slat_drag = slat_profile + slat_interference
    else:
        flap_drag = 0
        slat_drag = 0
    
    # Landing Gear
    gear = CD0_gear(x_nlg=lg.longitudinal_nlg, height_nlg=lg.height_nlg,
                    n_main_wheels=lg.n_wheels_mlg,
                    n_nose_wheels=lg.n_wheels_nlg,
                    nlg_tire_diameter=lg.selected_nlg_tire["Outside Diameter Max (In)"] * 2.54 / 100,
                    nlg_tire_width=lg.selected_nlg_tire["Section Width Max (In)"] * 2.54 / 100,
                    mlg_tire_diameter=lg.selected_mlg_tire["Outside Diameter Max (In)"] * 2.54 / 100,
                    mlg_tire_width=lg.selected_mlg_tire["Section Width Max (In)"] * 2.54 / 100,
                    S=S)

    CD0 = wing + fuselage + h_tail + v_tail + nacelle_isolated + nacelle_interference + propeller + flap_drag + slat_drag + gear
    # print(f' Speed: {speed}, density: {density}, dynamic pressure: {0.5*density*speed**2}')
    print(f' \n CD0 overview: \n wing: {wing} \n fuselage: {fuselage} \n ht: {h_tail} \n vt: {v_tail} \n isolated nacelle: {nacelle_isolated} \n nacelle interference: {nacelle_interference} \n propeller: {propeller} \n flap drag: {flap_drag} \n slat drag: {slat_drag} \n landing gear: {gear} \n \n total: {CD0}')
    if update_ac:
        w.CD0 = CD0
    return CD0

def C_D_L(ac: Aircraft, 
          CD0,
          flight_condition: str = 'cruise', # 'cruise' or 'landing' or 'take-off'
          update_ac: bool = False,
          wing_tip: bool = False):
    S = ac.wing.area
    if flight_condition == 'cruise':
        temp_shift = 0
        alt = ac.requirements.cruise['cr_altitude'] * FT_TO_M
        speed = ac.requirements.cruise['cr_speed'] * KTS_TO_MS
        Atm = Atmosphere(alt, temp_shift)
        density = float(Atm.density)
        mass_frac = ac.requirements.cruise['cr_mass_frac']
        C_L = mass_frac * ac.weights.m_takeoff * g / (0.5 * density * speed**2 * S)
        flap_deflection = 0
    if flight_condition == 'take-off':
        C_L = ac.hld_and_ailerons.take_off_lift['CL_max']
        temp_shift = ac.requirements.take_off['to_temp_shift']
        alt = ac.requirements.take_off['to_altitude'] * FT_TO_M
        Atm = Atmosphere(alt, temp_shift)
        density = Atm.density
        mass_frac = ac.requirements.take_off['to_mass_frac']
        speed = np.sqrt(mass_frac * ac.weights.m_takeoff / (0.5 * density * C_L * S))
        flap_deflection = ac.hld_and_ailerons.flaps['to_deflection']
    if flight_condition == 'landing':
        temp_shift = ac.requirements.landing['la_temp_shift']
        alt = ac.requirements.landing['la_altitude'] * FT_TO_M
        Atm = Atmosphere(alt, temp_shift)
        density = Atm.density
        mass_frac = ac.requirements.landing['la_mass_frac']
        speed = 1.3 * ac.requirements.general['stall_speed'] * KTS_TO_MS
        C_L = mass_frac * ac.weights.m_takeoff * g / (0.5 * density * speed**2 * S)
        flap_deflection = ac.hld_and_ailerons.flaps['ld_deflection']
    
    # General
    A = ac.wing.aspect_ratio
    A_wing_tip = 0.004  # NOTE: add wing tip effect here
    if wing_tip:
        A_wing_tip = ...
    A_eff = A + A_wing_tip 
    c_r = ac.wing.c_root
    b = ac.wing.span
    taper = ac.wing.taper_ratio
    if ac.wing.sweep == 0:
        e = 1.78 * (1 - 0.045 * A_eff**0.68) - 0.64
    else:
        e = 4.61 * (1 - 0.045 * A_eff**0.68) * (np.cos(np.deg2rad(LE_sweep_deg(ac.wing.sweep, c_r, b, taper))))**0.15 - 3.1

    # Flap
    if flap_deflection != 0:
        e += 0.0046 * flap_deflection

    K = 1 / (np.pi * A_eff * e)
    CDi = C_L**2 * K
    tip_twist = ac.wing.tip_twist  # degrees
    ld = 0.5 * np.sqrt(1 / (K * CD0))
    if tip_twist != 0:
        CDi += 0.00004 * 2 / 3 * tip_twist
    if update_ac:
        ac.wing.e = e
        ac.wing.k = K
        ac.wing.ld = ld
        ac.wing.aspect_ratio = A_eff
    return CDi, e, K, ld