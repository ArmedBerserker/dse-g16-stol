from classes.aircraft_2 import Aircraft
from c2_m import Snet, LE_sweep_deg, sweep_at_x_c_deg, lift_slope, S_wf, closest_value, chord_at_y_span
from lookups.consts import *
import pandas as pd
import numpy as np
from scipy.interpolate import RegularGridInterpolator, interp1d
from scipy.optimize import brentq
from classes.isa import Atmosphere
import matplotlib.pyplot as plt

''' TO DO:
    - change S_b_nac eqn'''

def mu_air(T):
    return 1.81 * 1e-5 * (T / 293.15)**1.5 * (293.15 + 110.4) / (T + 110.4)

def exposed_wing_mgc(S, b, taper, c_r, b_f):
    S_net = Snet(S, b_f, taper, b, c_r)
    return S_net / (b - b_f)

def CD0_w(R_wf, R_LS, C_f_w, l_dash, t_c_max, S_wet_w, S):
    return R_wf * R_LS * C_f_w * (1 + l_dash * t_c_max + 100 * t_c_max**4) * S_wet_w / S

def CD0_fuselage(density, speed, l_f, mu, S_fus, # max fuselage cross section area
                 S_b_fus, S, M, R_N_fus, S_wet_fus): # density, speed, l_nac, mu, S_nac_max, S_b_nac, S, M, R_N_nac, S_wet_fus=S_wet_nac
    R_wf = Rwf(M, R_N_fus)
    C_f_fus = C_f(R_N=density * speed * l_f / mu, M=M)
    print(f'R_wf: {R_wf}, C_f_fus: {C_f_fus}')
    d_b = np.sqrt(4 * S_b_fus / np.pi)
    d_f_eq = np.sqrt(4 * S_fus / np.pi)
    # print(f'R_wf: {R_wf}, C_f_fus: {C_f_fus}, d_b: {d_b}, d_f_eq: {d_f_eq}, density: {density}, speed: {speed}, l_f: {l_f}, mu: {mu}, S_fus: {S_fus}, S_b_fus: {S_b_fus}, S: {S}, M: {M}, R_N_fus: {R_N_fus}, S_wet_fus: {S_wet_fus}')
    C_D0_b_fus = R_wf * C_f_fus * (1 + 60 / (l_f / d_f_eq)**3 + 0.0025 * (l_f / d_f_eq)) * S_wet_fus / S
    C_D_b_fuse = (0.029 * (d_b / d_f_eq)**3 / (C_D0_b_fus * (S / S_fus))**0.5) * (S / S_fus)
    return C_D0_b_fus + C_D_b_fuse

def S_wet_wing(t_c_r, t_c_t, S, b_f, taper, b, c_r):
    tau = t_c_t / t_c_r
    # print(f'Snet: {Snet(S, b_f, taper, b, c_r)}')
    # return 2 * Snet(S, b_f, taper, b, c_r) * (1 + 0.25 * t_c_r * (1 + tau * taper) / (1 + taper))
    return 2 * 1.07 * Snet(S, b_f, taper, b, c_r)

def C_f(R_N, M):  # DATCOM fig 4.1.5.1-26
    x = np.log10(R_N)
    CF3 = 3.92725e-6 * x**5 -1.30370e-4 * x**4 -1.65388e-3 * x**3 -9.59519e-3 * x**2 +2.18366e-2 * x 
    CF0 = 4.12963e-6 * x**5 -1.36204e-4 * x**4 + 1.71620e-3 * x**3 -9.88935e-3 * x**2 +2.23641e-2 * x
    if M > 0.3:
        return CF0
    else: 
        return CF0 # - (CF0 - CF3) * (M / 0.3)

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

def Rwf(M, R_N_fus):
    M = closest_value(M, values=[0.25, 0.4])
    if M == 0.25:
        return interp_value(pd.read_csv('lookups/roskam_p6_fig_4_1_rwf.csv'), R_N_fus, x_col='Fuselage Reynolds Number (M = 0.25)', y_col='Wing-Fuse Interference Factor (M = 0.25)', log_x=True)
    else:
        return interp_value(pd.read_csv('lookups/roskam_p6_fig_4_1_rwf_2.csv'), R_N_fus, x_col='Fuselage Reynolds Number (M = 0.4)', y_col='Wing-Fuse Interference Factor (M = 0.4)', log_x=True)
    

def C_D0(ac: Aircraft, 
         n_engine_operative: int, #  flap_type: str, # 'split' or 'plain' or 'slotted' 'fowler' or 'krueger'
         flight_condition: str = 'cruise', # 'cruise' or 'landing' or 'take-off'
         update_ac: bool = False
         ):
    S = ac.wing.area
    flap_type = ac.hld_and_ailerons.flaps['flap_type']
    nacelle_on_top_of_wing: bool = ac.engine.eng_above_wing
    if flight_condition == 'cruise':
        temp_shift = 0
        alt = ac.requirements.cruise['cr_altitude'] * FT_TO_M
        speed = ac.requirements.cruise['cr_speed'] * KTS_TO_MS
        Atm = Atmosphere(alt, temp_shift)
        temp = float(Atm.temp)
        density = float(Atm.density)
    if flight_condition == 'take-off':
        C_L = ac.hld_and_ailerons.take_off_lift['CL_max'] * 0.95
        temp_shift = ac.requirements.take_off['to_temp_shift']
        alt = ac.requirements.take_off['to_altitude'] * FT_TO_M
        Atm = Atmosphere(alt, temp_shift)
        density = float(Atm.density)
        mass_frac = ac.requirements.cruise['to_mass_frac']
        speed = np.sqrt(mass_frac * ac.weights.m_takeoff / (0.5 * density * C_L * S))
        flap_deflection = ac.hld_and_ailerons.flaps['to_deflection']
    if flight_condition == 'take-off':
        temp_shift = ac.requirements.landing['la_temp_shift']
        alt = ac.requirements.landing['la_altitude'] * FT_TO_M
        Atm = Atmosphere(alt, temp_shift)
        density = float(Atm.density)
        mass_frac = ac.requirements.landing['la_mass_frac']
        speed = 1.3 * ac.requirements.general['stall_speed'] * KTS_TO_MS
        flap_deflection = ac.hld_and_ailerons.flaps['ld_deflection']
        
    # Wing
    mu = float(mu_air(temp))
    M = speed / np.sqrt(1.4 * 287 * temp)
    b = ac.wing.span
    sweep_c_4_deg = ac.wing.sweep
    taper = ac.wing.taper_ratio
    c_r = ac.wing.c_root
    b_f = ac.fuselage.width
    l_f = ac.fuselage.length
    c_w_e = exposed_wing_mgc(S, b, taper, c_r, b_f)
    # R_N_w = density * speed * c_w_e / mu
    R_N_fus = density * speed * l_f / mu
    R_wf = Rwf(M, R_N_fus)
    C_f_w = C_f(R_N=density * speed * c_w_e / mu, M=M)
    l_dash = 1.2
    t_c_max = ac.wing.t_c_max # at mean geometric chord
    x_c_t_c_max = ac.wing.x_c_t_c_max
    if x_c_t_c_max <= 0.3:
        l_dash = 2.0
    t_c_r = t_c_max
    t_c_t = t_c_max
    sweep_t_c_max_deg = sweep_at_x_c_deg(ac.wing.sweep_LE_deg, c_r, b, taper, x_c=x_c_t_c_max)
    R_LS = interp_value(pd.read_csv('lookups/roskam_p6_fig_4_2_rls.csv'), np.cos(np.deg2rad(sweep_t_c_max_deg)), 'cos(quarter chord)', 'lifting surface correction', log_x=False)
    S = ac.wing.area
    print(f'Swet_wing: {S_wet_wing(t_c_r, t_c_t, S, b_f, taper, b, c_r)}')
    print(R_wf, R_LS, C_f_w, l_dash, t_c_max, S_wet_wing(t_c_r, t_c_t, S, b_f, taper, b, c_r), S)
    wing_drag = CD0_w(R_wf, R_LS, C_f_w, l_dash, t_c_max, S_wet_wing(t_c_r, t_c_t, S, b_f, taper, b, c_r), S)  # R_wf * R_LS * C_f_w * (1 + l_dash * t_c_max + 100 * t_c_max**4) * S_wet_w / S

    # Fuselage
    d_f_max_eq = np.sqrt(4 * ac.fuselage.max_cross_section_area / np.pi)
    d_f_b_eq = np.sqrt(4 * ac.fuselage.base_area / np.pi)
    f = ac.fuselage
    S_fus = f.max_cross_section_area  # max fuselage cross section area
    S_b_fus = f.base_area
    S_wet_fus = np.pi * d_f_max_eq / 2 * (1.08 * (ac.fuselage.nose_cone_length + ac.fuselage.tail_cone_length) + 2 * (l_f - ac.fuselage.nose_cone_length - ac.fuselage.tail_cone_length))
    fuselage_drag = CD0_fuselage(density, speed, l_f, mu, S_fus, S_b_fus, S, M, R_N_fus=R_N_fus, S_wet_fus=S_wet_fus)

    # VT
    vt = ac.empennage.vertical_tail
    vt_sweep_c_4_deg = ac.empennage.vertical_tail['sweep']
    c_r_vt = vt['c_r_v']
    c_t_vt = vt['c_t_v']
    taper_vt = vt['taper_ratio']
    b_vt = vt['b_v']
    S_vt = vt['area']
    vt_t_c_max = ac.empennage.vertical_tail['t_c_max']
    vt_x_c_t_c_max = ac.empennage.vertical_tail['loc_t_c_max']
    vt_sweep_t_c_max_deg = sweep_at_x_c_deg(LE_sweep_deg(vt_sweep_c_4_deg, c_r_vt, b_vt, taper_vt), c_r_vt, b_vt, taper_vt, x_c=vt_x_c_t_c_max)
    R_LS_vt = interp_value(pd.read_csv('lookups/roskam_p6_fig_4_2_rls.csv'), np.cos(np.deg2rad(vt_sweep_t_c_max_deg)), 'cos(quarter chord)', 'lifting surface correction', log_x=False)
    b_f_vt = ac.fuselage.height - (d_f_max_eq - d_f_b_eq) * (ac.empennage.vertical_tail['x_v_frac_lf'] - 1 + ac.fuselage.tail_cone_fuselage_ratio) / ac.fuselage.tail_cone_fuselage_ratio  # fuselage width at vt intersection position
    C_f_vt = C_f(R_N=density * speed * exposed_wing_mgc(S_vt, b_vt, taper_vt, c_r_vt, b_f_vt) / mu, M=M)
    l_dash_vt = 1.2
    if vt_x_c_t_c_max < 0.3:
        l_dash_vt = 2.0
    t_c_r_vt = vt_t_c_max
    t_c_t_vt = vt_t_c_max
    vt_drag = CD0_w(R_wf=1.0, R_LS=R_LS_vt, C_f_w=C_f_vt, l_dash=l_dash_vt, t_c_max=vt_t_c_max, S_wet_w=S_wet_wing(t_c_r_vt, t_c_t_vt, S_vt, b_f_vt, taper_vt, b_vt, c_r_vt), S=S_vt)

    # HT
    t_tail_condition: bool = ac.empennage.t_tail_condition  # t-tail or not
    ht = ac.empennage.horizontal_tail
    c_r_ht = ht['c_r_h']
    taper_ht = ht['taper_ratio']
    b_ht = ht['b_h']
    S_ht = ht['area']
    ht_t_c_max = ac.empennage.horizontal_tail['t_c_max']
    ht_x_c_t_c_max = ac.empennage.horizontal_tail['loc_t_c_max']
    ht_sweep_t_c_max_deg = sweep_at_x_c_deg(ac.wing.sweep_LE_deg, c_r_ht, b_ht, taper_ht, x_c=ht_x_c_t_c_max)
    R_LS_ht = interp_value(pd.read_csv('lookups/roskam_p6_fig_4_2_rls.csv'), np.cos(np.deg2rad(ht_sweep_t_c_max_deg)), 'cos(quarter chord)', 'lifting surface correction', log_x=False)
    b_f_ht = c_t_vt * t_c_t_vt / 2  # fuselage width at ht intersection position
    if not t_tail_condition:
        b_f_ht = b_f - (d_f_max_eq - d_f_b_eq) * (ac.empennage.horizontal_tail['x_h_frac_lf'] - 1 + ac.fuselage.tail_cone_fuselage_ratio) / ac.fuselage.tail_cone_fuselage_ratio
    R_N_ht = density * speed * exposed_wing_mgc(S_ht, b_ht, taper_ht, c_r_ht, b_f_ht) / mu
    if not t_tail_condition:
        R_N_ht *= np.sqrt(0.85)
    C_f_ht = C_f(R_N_ht, M)
    l_dash_ht = 1.2
    if ht_x_c_t_c_max < 0.3:
        l_dash_ht = 2.0
    t_c_r_ht = ht_t_c_max
    t_c_t_ht = ht_t_c_max
    ht_drag = CD0_w(R_wf=1.0, R_LS=R_LS_ht, C_f_w=C_f_ht, l_dash=l_dash_ht, t_c_max=ht_t_c_max, S_wet_w=S_wet_wing(t_c_r_ht, t_c_t_ht, S_ht, b_f_ht, taper_ht, b_ht, c_r_ht), S=S_ht)


    # Nacelle/pylon
    n_eng = ac.engine.count
    l_nac = ac.engine.length_nac  # nacelle length
    S_nac_max = ac.engine.nac_diameter**2 / 4 * np.pi  # max nacelle cross section area
    S_b_nac = (ac.engine.nac_diameter * 0.5)**2 / 4 * np.pi  # nacelle base area
    S_wet_nac = ac.engine.nac_diameter * np.pi * l_nac + 2 * S_b_nac
    R_N_nac = density * speed * l_nac / mu
    isolated_nac = CD0_fuselage(density, speed, l_nac, mu, S_nac_max, S_b_nac, S, M, R_N_nac, S_wet_fus=S_wet_nac)
    isolated_nac = 0 # .06 * S_nac_max / S
    # print(f'Isolated nacelle: {isolated_nac}, {density, speed, l_nac, mu, S_nac_max, S_b_nac, S, M, R_N_nac, S_wet_nac}')

    c_nac = chord_at_y_span(c_r, taper, y=b_f/2+ac.engine.eng_y_pos_fuselage, b=b)  # chord at nacelle 
    c_r_nac = 0.25 * c_nac
    b_nac = ac.engine.eng_vdist_from_wing_y_c * c_nac  # nacelle width
    nac_t_c_max = ac.engine.nac_t_c_max
    nac_x_c_t_c_max = ac.engine.nac_x_c_t_c_max
    nac_sweep_t_c_max_deg = sweep_at_x_c_deg(LE_sweep_deg(sweep_c_4=0, c_r=c_r_nac, b=b_nac, taper_ratio=0), c_r_nac, b_nac, taper_ratio=0, x_c=nac_x_c_t_c_max)
    R_LS_nac = interp_value(pd.read_csv('lookups/roskam_p6_fig_4_2_rls.csv'), np.cos(np.deg2rad(nac_sweep_t_c_max_deg)), 'cos(quarter chord)', 'lifting surface correction', log_x=False)
    C_f_nac = C_f(R_N=density * speed * l_nac / mu, M=M)
    l_dash_nac = 1.2
    if nac_x_c_t_c_max < 0.3:
        l_dash_nac = 2.0
    isolated_pylon = CD0_w(R_wf=1.0, R_LS=R_LS_nac, C_f_w=C_f_nac, l_dash=l_dash_nac, t_c_max=nac_t_c_max, S_wet_w=S_wet_wing(nac_t_c_max, nac_t_c_max, S=b_nac*c_r_nac, b_f=t_c_max*(c_r + taper*c_r)/2, taper=0, b=b_nac, c_r=c_r_nac), S=b_nac*c_r_nac)
    isolated_pylon = 0
    print(f'Isolated pylon: {isolated_pylon}')
    i_n = ac.engine.i_n  # nacelle incidence angle [deg]
    D_cl_1 = -0.3
    if nacelle_on_top_of_wing:
        D_cl_1 = 0.2
    D_cl_2 = 0.056 * i_n
    wing_nac_interference = 0.036 * (c_nac * b_nac / S) * (D_cl_1 + D_cl_2)**2
    print(f'Wing nacelle interference: {wing_nac_interference}')

    SHP = ac.engine.power_cr / HP_TO_W  # shaft horse power
    D_prop =  ac.engine.prop_diameter # propeller diameter
    wind_milling = 33 / (0.5 * density * speed**2 * PA_TO_LBSpFT2 * S) * SHP / (speed * MpS_TO_FpS)
    wind_milling = 0
    wind_milling_inoperative = 0
    if n_engine_operative != n_eng:
        wind_milling_inoperative = 0.00125 * ac.engine.eta_prop * D_prop**2 / S
    propulsion_drag = n_eng * (wind_milling + wing_nac_interference + isolated_pylon + isolated_nac) + (n_eng - n_engine_operative) * wind_milling_inoperative
    propulsion_drag = 0
    # Flap:
    if flight_condition != 'cruise':
        y_start_flap = ac.hld_and_ailerons.flaps['y_flap_in']
        y_end_flap = ac.hld_and_ailerons.flaps['y_flap_out']
        cf_c = ac.hld_and_ailerons.flaps['cf_c']  # flap chord length / chord length
        if flap_type == 'split':  # 'split' or 'plain' or 'slotted' 'fowler' or 'krueger'
            fd = closest_value(t_c_max, values=[10, 20, 30])
            D_CD_flap_stuff = interp_value(pd.read_csv(f'lookups/t_c_0.{fd}.csv'), cf_c, f'cf/c ({flight_condition})', f'dCdp ({flight_condition})', log_x=False)
        elif flap_type == 'plain':
            fd = closest_value(flap_deflection, values=[15, 60])
            D_CD_flap_stuff = interp_value(pd.read_csv(f'lookups/d_f_{fd}.csv'), cf_c, 'cf/c', 'dCdp', log_x=False)
        elif flap_type == 'slotted':
            fd = closest_value(cf_c, values=[0.1, 0.2, 0.3])
            D_CD_flap_stuff = interp_value(pd.read_csv('lookups/cf_c_comb2.csv'), flap_deflection, 'df', f'dCdp (cf={fd}0)', log_x=False)
        elif flap_type == 'fowler':
            fd = closest_value(cf_c, values=[0.1, 0.2, 0.3, 0.4])
            D_CD_flap_stuff = interp_value(pd.read_csv('lookups/roskam_p6_fig_4_48.csv'), flap_deflection, f'df(cf/c={fd})', f'dCdp(cf/c={fd})', log_x=False)
        elif flap_type == 'kruger':
            D_CD_flap_stuff = wing_drag * (cf_c * np.cos(np.deg2rad(flap_deflection)) + 1)
        else: 
            raise ValueError(f"Flap type give: {flap_type}, check possible entries for D_CD_flap_stuff variable")
        flap_profile = D_CD_flap_stuff * np.cos(np.deg2rad(sweep_c_4_deg)) * S_wf(y_start_flap, y_end_flap, taper, c_r, b) / S

        # b_fi_b = y_start_flap * 2 / b
        # b_fo_b = y_end_flap * 2 / b
        # K = 
        # Delta_CL_max_flapped = 
        # induced_flap = K**2 * Delta_CL_max_flapped**2 * np.cos(np.deg2rad(sweep_c_4_deg))

        interference_flap = flap_profile * flap_interference_factor(flap_type)

        flap_drag = interference_flap + flap_profile # + induced_flap
    else: 
        flap_drag = 0

    # Landing gear
    w_tire = ac.landing_gear.selected_mlg_tire["Section Width Max (In)"] * 2.54 / 100  # tire width NOTE: check all lg dictionary names
    d_tire = ac.landing_gear.selected_mlg_tire["Outside Diameter Max (In)"] * 2.54 / 100  # tire diameter
    w_strut = 0.5 * w_tire
    l_strut = np.abs(ac.landing_gear.height_mlg) - d_tire / 2
    m = (w_tire * d_tire + l_strut * w_strut) / ((w_tire + w_strut) * (l_strut + 0.5 * d_tire))
    lg_drag = ((w_tire + w_strut) * (l_strut + 0.5 * d_tire) / S) * 0.04955 * np.exp(5.615 * m)

    # Miscelaneous (+5%)
    C_D0 = (wing_drag + fuselage_drag + ht_drag + vt_drag + propulsion_drag + flap_drag + lg_drag) * 1.05
    print(f'CD0 components: {wing_drag, fuselage_drag, ht_drag, vt_drag, propulsion_drag, flap_drag, lg_drag}')

    if update_ac:
        ac.wing.CD0 = C_D0
    return C_D0

def C_D_L(ac:Aircraft, 
          flight_condition: str = 'cruise', # 'cruise' or 'landing' or 'take-off'
          update_ac: bool = False,
          wing_tip: bool = False
          ):
    S = ac.wing.area
    if flight_condition == 'cruise':
        temp_shift = 0
        alt = ac.requirements.cruise['cr_altitude'] * FT_TO_M
        speed = ac.requirements.cruise['cr_speed'] * KTS_TO_MS
        Atm = Atmosphere(alt, temp_shift)
        density = Atm.density
        mass_frac = ac.requirements.cruise['cr_mass_frac']
        C_L = mass_frac * ac.weights.m_takeoff * g / (0.5 * density * speed**2 * S)
        flap_deflection = 0
    if flight_condition == 'take-off':
        C_L = ac.hld_and_ailerons.take_off_lift['CL_TO']
        temp_shift = ac.requirements.take_off['to_temp_shift']
        alt = ac.requirements.take_off['to_altitude'] * FT_TO_M
        Atm = Atmosphere(alt, temp_shift)
        density = Atm.density
        mass_frac = ac.requirements.cruise['to_mass_frac']
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
    A_wing_tip = 0  # NOTE: add wing tip effect here
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

    # Wing

    # Fuselage

    # Empennage

    # Nacelle/pylon

    # Flap
    if flap_deflection != 0:
        e += 0.0046 * flap_deflection

    K = 1 / (np.pi * A_eff * e)
    CDi = C_L**2 * K
    tip_twist = ac.wing.tip_twist  # degrees
    if tip_twist != 0:
        CDi += 0.00004 * 2 / 3 * tip_twist
    if update_ac:
        ac.wing.e = e
        ac.wing.k = K
    return CDi, e, K
