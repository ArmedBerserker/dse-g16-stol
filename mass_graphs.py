from classes.aircraft_2 import Aircraft, loader, Requirements, Mission, Fuselage, Wing, Engine, Weights, Empennage, HLD_and_AIL, Landing_Gear
from lookups.consts import *
import pandas as pd
import numpy as np
from scipy.interpolate import RegularGridInterpolator, interp1d
from scipy.optimize import brentq
from classes.isa import Atmosphere
import matplotlib.pyplot as plt
import os
import matplotlib.cm as cm

def mass_pie_charts(w_fc, w_supercap, w_propeller, w_engine, w_fuel_system, w_hps_els, w_api, w_paint, w_furnishing, w_iae, W_tfo, W_fuel, wwing, w_ht, w_vt, wfus, wnac, w_mlg, w_nlg, base_file_name: str = 'final_pie_chart', show_pie_chart: bool = False, W_PL: float = 662.0):
    w_power = w_propeller + w_engine + w_fuel_system
    w_fxeq = w_hps_els + w_api + w_paint + w_furnishing + w_iae
    w_structure = wwing + w_ht + w_vt + wfus + wnac + w_mlg + w_nlg
    W_oe = w_structure + w_power + w_fxeq
    categories = ['Structural', 'Power', 'Fixed equipment']
    raw_values = [w_structure, w_power, w_fxeq]
    values = np.array(raw_values) / W_oe * 100
    fig, ax = plt.subplots(figsize=(8, 6))
    wedges, _ = ax.pie(values, startangle=90)
    total = sum(values)
    for i, wedge in enumerate(wedges):
        angle = (wedge.theta2 + wedge.theta1) / 2
        x = np.cos(np.radians(angle)) * 0.7
        y = np.sin(np.radians(angle)) * 0.7
        percentage = values[i] / total * 100
        ax.text(x, y,
                f'{categories[i]}\n{percentage:.1f}%\n({raw_values[i]:.1f} kg)',
                ha='center', va='center', fontsize=8)
    ax.axis('equal')
    plt.title('Distribution of OEW')
    plt.savefig(f'outputs/Final_report_graphs/{base_file_name}_OEW.png', dpi=400)
    if show_pie_chart:
        plt.show()

    # Structural 
    categories = ['Wing', 'Horizontal tail', 'Vertical tail', 'Fuselage', 'Nacelles', 'Landing gear']
    raw_values = [wwing, w_ht, w_vt, wfus, wnac, w_mlg+w_nlg]
    values = np.array(raw_values) / w_structure * 100
    fig, ax = plt.subplots(figsize=(8, 6))
    wedges, _ = ax.pie(values, startangle=90)
    total = sum(values)
    for i, wedge in enumerate(wedges):
        angle = (wedge.theta2 + wedge.theta1) / 2
        x = np.cos(np.radians(angle)) * 0.7
        y = np.sin(np.radians(angle)) * 0.7
        percentage = values[i] / total * 100
        ax.text(x, y,
                f'{categories[i]}\n{percentage:.1f}%\n({raw_values[i]:.1f} kg)',
                ha='center', va='center', fontsize=8)
    ax.axis('equal')
    plt.title('Distribution of structural weight')
    plt.savefig(f'outputs/Final_report_graphs/{base_file_name}_struc.png', dpi=400)
    if show_pie_chart:
        plt.show()

    # MTOW
    categories = ['Empty weight', 'Payload', 'Trapped fuel and oil', 'Fuel']
    raw_values = [W_oe, W_PL, W_tfo, W_fuel]
    values = np.array(raw_values) # / Wto * 100
    fig, ax = plt.subplots(figsize=(8, 6))
    wedges, _ = ax.pie(values, startangle=90)
    total = sum(values)
    for i, wedge in enumerate(wedges):
        angle = (wedge.theta2 + wedge.theta1) / 2
        x = np.cos(np.radians(angle)) * 0.7
        y = np.sin(np.radians(angle)) * 0.7
        if categories[i]=='Trapped fuel and oil':
            x = np.cos(np.radians(angle))
        percentage = values[i] / total * 100
        ax.text(x, y,
                f'{categories[i]}\n{percentage:.1f}%\n({raw_values[i]:.1f} kg)',
                ha='center', va='center', fontsize=8)
    ax.axis('equal')
    plt.title('Distribution of MTOW')
    plt.savefig(f'outputs/Final_report_graphs/{base_file_name}_MTOW.png', dpi=400)
    if show_pie_chart:
        plt.show()
    
    # Fixed equipment
    categories = ['Flight controls', 'Supercapacitors', 'HPS & ELS', 'API', 'Paint', 'Furnishing', 'IAE']
    raw_values = [w_fc, w_supercap, w_hps_els, w_api, w_paint, w_furnishing, w_iae]
    values = np.array(raw_values) # / Wto * 100
    fig, ax = plt.subplots(figsize=(8, 6))
    wedges, _ = ax.pie(values, startangle=90)
    total = sum(values)
    for i, wedge in enumerate(wedges):
        angle = (wedge.theta2 + wedge.theta1) / 2
        x = np.cos(np.radians(angle)) * 0.7
        y = np.sin(np.radians(angle)) * 0.7
        if categories[i]=='Trapped fuel and oil':
            x = np.cos(np.radians(angle))
        percentage = values[i] / total * 100
        ax.text(x, y,
                f'{categories[i]}\n{percentage:.1f}%\n({raw_values[i]:.1f} kg)',
                ha='center', va='center', fontsize=8)
    ax.axis('equal')
    plt.title('Distribution of fixed equiment mass')
    plt.savefig(f'outputs/Final_report_graphs/{base_file_name}_FXEQ.png', dpi=400)
    if show_pie_chart:
        plt.show()

    # Power
    categories = ['Propellers', 'Engines', 'Fuel system']
    raw_values = [w_propeller, w_engine, w_fuel_system]
    values = np.array(raw_values) # / Wto * 100
    fig, ax = plt.subplots(figsize=(8, 6))
    wedges, _ = ax.pie(values, startangle=90)
    total = sum(values)
    for i, wedge in enumerate(wedges):
        angle = (wedge.theta2 + wedge.theta1) / 2
        x = np.cos(np.radians(angle)) * 0.7
        y = np.sin(np.radians(angle)) * 0.7
        if categories[i]=='Trapped fuel and oil':
            x = np.cos(np.radians(angle))
        percentage = values[i] / total * 100
        ax.text(x, y,
                f'{categories[i]}\n{percentage:.1f}%\n({raw_values[i]:.1f} kg)',
                ha='center', va='center', fontsize=8)
    ax.axis('equal')
    plt.title('Distribution of power system mass')
    plt.savefig(f'outputs/Final_report_graphs/{base_file_name}_power.png', dpi=400)
    if show_pie_chart:
        plt.show()

def scissor_plot(x_le, x_cg_fwd, x_cg_aft, Cmac, x_ac, mac, l_h, C_L_alpha_h, 
                 C_L_alpha_A_less_h, de_da, C_L_h, C_L_A_less_h, 
                 Vh_V2: float = 1, SM: float = 0.05, show_plot: bool = False, 
                 save_path: str = 'outputs/Final_report_graphs/scissor_plot.png'):
    x_cg_arr = np.linspace(-1.5, 1.5, 150)
    
    x_cg_aft = (x_cg_aft - x_le) / mac
    x_cg_fwd = (x_cg_fwd - x_le) / mac
    controllability = (x_cg_arr - x_ac + Cmac / C_L_A_less_h) / (C_L_h / C_L_A_less_h * l_h / mac * Vh_V2)
    stability = (x_cg_arr - x_ac + SM) / (C_L_alpha_h / C_L_alpha_A_less_h * (1 - de_da) * l_h / mac * Vh_V2)
    nstability = (x_cg_arr - x_ac) / (C_L_alpha_h / C_L_alpha_A_less_h * (1 - de_da) * l_h / mac * Vh_V2)
    Sh_S_cont = (x_cg_fwd - x_ac + Cmac / C_L_A_less_h) / (C_L_h / C_L_A_less_h * l_h / mac * Vh_V2)
    Sh_S_stab = (x_cg_aft - x_ac + SM) / (C_L_alpha_h / C_L_alpha_A_less_h * (1 - de_da) * l_h / mac * Vh_V2)
    
    plt.figure()
    plt.plot(x_cg_arr, controllability, color='blue', label='controllability')
    plt.plot(x_cg_arr, stability, color='cyan', label='stability (SM=0.05)')
    plt.plot(x_cg_arr, nstability, linestyle='--', color='grey', label='neutral stability')
    if Sh_S_cont > Sh_S_stab:
        plt.scatter(x_cg_fwd, Sh_S_cont, color='blue', s=50, zorder=5)
        plt.annotate(f'({x_cg_fwd:.3f}, {Sh_S_cont:.3f})',
                    (x_cg_fwd, Sh_S_cont),
                    xytext=(10, 10),
                    textcoords='offset points',
                    color='blue')
        plt.axhline(Sh_S_cont, color='blue', linestyle=':', alpha=0.6)
        plt.axvline(x_cg_fwd, color='blue', linestyle=':', alpha=0.6)
    else:
        plt.scatter(x_cg_aft, Sh_S_stab, color='cyan', s=50, zorder=5)
        plt.annotate(f'({x_cg_aft:.3f}, {Sh_S_stab:.3f})',
                    (x_cg_aft, Sh_S_stab),
                    xytext=(10, -15),
                    textcoords='offset points',
                    color='cyan')
        plt.axhline(Sh_S_stab, color='cyan', linestyle=':', alpha=0.6)
        plt.axvline(x_cg_aft, color='cyan', linestyle=':', alpha=0.6)

    plt.ylim((0, 0.6))
    plt.legend(loc='upper left')
    plt.xlabel(f'$\\frac{{x_{{cg}}-LEMAC}}{{\\bar{{c}}}}$')
    plt.ylabel(f'$\\frac{{S_h}}{{S}}$')
    plt.title('Scissor plot')
    plt.savefig(save_path, dpi=400)
    if show_plot:
        plt.show()

if __name__ == '__main__':
    Cmac = -0.14728
    VhV2 = 1
    de_da = 0
    x_le = 4.4
    x_cg_fwd = 5.051
    x_cg_aft = 5.19
    x_ac = 0.231752
    mac = 1.87
    l_h = 5.2
    C_L_alpha_h = 4.370966148
    C_L_alpha_A_less_h = 5.348506148
    C_L_h = -0.635992207
    C_L_A_less_h = 0.997

    scissor_plot(x_le, x_cg_fwd, x_cg_aft, Cmac, x_ac, mac, l_h, C_L_alpha_h, C_L_alpha_A_less_h, de_da, C_L_h, C_L_A_less_h, show_plot=True)

    MTOW = 1848
    W_tfo = 0.005 * MTOW
    W_fuel = 168.29
    W_PL = 662 

    W_fus = 152.005
    W_wing = 161.2
    w_ht = 20
    w_vt = 33.48
    wnac = 11.27
    w_mlg = 33
    w_nlg = 8

    w_struc = W_fus + W_wing + w_ht + w_vt + wnac + w_mlg + w_nlg

    w_propeller = 29
    w_engine = 224.4
    w_fuel_system = 22.7

    w_hps_els = 51.44
    w_api = 110 / 2.205
    w_paint = 5.5
    w_furnishing = 43.08
    w_iae = 32.07
    w_supercap = 105
    w_fc = 38

    w_power = w_propeller + w_engine + w_fuel_system
    w_fxeq = w_hps_els + w_api + w_paint + w_furnishing + w_iae + w_supercap + w_fc

    print(f'OEW sum: {w_power + w_fxeq + w_struc}')

    mass_pie_charts(w_fc, w_supercap, w_propeller, w_engine, w_fuel_system, w_hps_els, w_api, w_paint, w_furnishing, w_iae, W_tfo, W_fuel, W_wing, w_ht, w_vt, W_fus, wnac, w_mlg, w_nlg)

