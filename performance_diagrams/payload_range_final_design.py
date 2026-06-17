import sys
import os

# Fix path FIRST, before any local imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import math
import numpy as np
import matplotlib.pyplot as plt
# import os
from pathlib import Path
from lookups.consts import *


def calculate_range(eta_eng, eta_p, L_D_max, e_f, m_e, m_pl, m_f, R_aux):
    R = eta_eng * eta_p * L_D_max * e_f / g * math.log((m_e + m_pl + m_f) / (m_e + m_pl)) - R_aux
    return R


def calculate_payload_range_parameters(m_mto, m_e, m_f_des, m_f_max, m_pl_max, m_pilot, A, e, C_d0,
                                       eta_p, eta_eng, e_f, R_des, V_cr):

    L_D_max = 0.5 * math.sqrt((math.pi * A * e) / (C_d0))
    t_e = 45 * 60
    R_eq_res = t_e * V_cr

    R_eq = R_des + R_eq_res
    R_aux = R_eq - R_des

    m_f_required_cruise_reserve = m_mto * (1 - math.exp((-R_eq) / (eta_eng * eta_p * e_f * L_D_max / g)))
    m_f_required_cruise = m_mto * (1 - math.exp((-R_des) / (eta_eng * eta_p * e_f * L_D_max / g)))
    m_f_required_reserve = m_mto * (1 - math.exp((-R_eq_res) / (eta_eng * eta_p * e_f * L_D_max / g)))
    m_f_required = m_mto * (1 - ((1 - (m_f_required_cruise_reserve) / m_mto) * 0.992 * 0.996 * 0.996 * 0.99 * 0.992 * 0.992))
    m_f_required_told = m_f_required - m_f_required_cruise_reserve
    print(m_f_required_cruise, m_f_required_reserve, m_f_required_cruise_reserve, m_f_required_told, m_f_required)

    if m_f_des > m_f_required:
        print(f"Design fuel mass ({m_f_des} kg) is sufficient to meet fuel requirements ({m_f_required:.2f} kg) for "
              f"design range ({R_des/1000} km)")
    else:
        print(f"Design fuel mass ({m_f_des} kg) is INSUFFICIENT to meet fuel requirements ({m_f_required:.2f} kg) "
              f"for design range ({R_des/1000} km)")

    design_point_1 = (0, m_pl_max)

    m_f_pl_max = m_mto - m_e - m_pl_max - m_f_required_told
    print(m_f_pl_max)
    R_pl_max = calculate_range(eta_eng, eta_p, L_D_max, e_f, m_e, m_pl_max, m_f_pl_max, R_aux)
    design_point_2 = (R_pl_max, m_pl_max)

    m_pl_f_max = m_mto - m_e - m_f_max
    m_f_required_cruise_reserve_pax = m_mto * (1 - math.exp((-(2150000+R_eq_res)) / (eta_eng * eta_p * e_f * L_D_max / g)))
    m_f_required_cruise_pax = m_mto * (1 - math.exp((-2150000) / (eta_eng * eta_p * e_f * L_D_max / g)))
    m_f_required_reserve_pax = m_mto * (1 - math.exp((-R_eq_res) / (eta_eng * eta_p * e_f * L_D_max / g)))
    m_f_required_pax = m_mto * (
                1 - ((1 - (m_f_required_cruise_reserve_pax) / m_mto) * 0.992 * 0.996 * 0.996 * 0.99 * 0.992 * 0.992))
    m_f_required_told_pax = m_f_required_pax - m_f_required_cruise_reserve_pax
    print(m_f_required_cruise_pax, m_f_required_reserve_pax, m_f_required_cruise_reserve_pax, m_f_required_told_pax, m_f_required_pax)

    R_f_max = calculate_range(eta_eng, eta_p, L_D_max, e_f, m_e, m_pl_f_max, m_f_max - m_f_required_told, R_aux)
    design_point_3 = (R_f_max, m_pl_f_max)

    R_ferry = calculate_range(eta_eng, eta_p, L_D_max, e_f, m_e, m_pilot, m_f_max - m_f_required_told, R_aux)
    design_point_4 = (R_ferry, m_pilot)

    design_point_5 = (R_ferry, 0)

    return [design_point_1, design_point_2, design_point_3, design_point_4, design_point_5]


def plot_payload_range(payload_range_points, m_pl_des, R_des,
                       save_path=Path(__file__).parent / "performance_figures/payload_range_diagram.png",
                       figsize=(10, 6)):

    ranges = [point[0] / 1000 for point in payload_range_points]
    payloads = [point[1] for point in payload_range_points]

    plt.plot(ranges[:4], payloads[:4], marker='.', linestyle='-',linewidth=2, color='green')
    plt.plot(ranges[3:5], payloads[3:5], marker='.', linestyle='--', linewidth=1, color='blue')
    plt.text(ranges[3] + 5, payloads[3] / 2, f'Practical ferry\nrange = {ranges[3]:.0f} km', color='blue',
             fontsize=9, horizontalalignment='left', verticalalignment='center')

    point_labels = ['A', 'B', 'C', 'D']

    for i, point_label in enumerate(point_labels):
        plt.text(ranges[i], payloads[i] + 10,f' {point_label}',fontsize=10)

    plt.hlines(y=m_pl_des, xmin=0, xmax=R_des/1000, colors='red',linestyle='--', linewidth=1)
    plt.text(R_des/(2*1000), m_pl_des,f'Design Payload = {m_pl_des:.0f} kg', color='red',fontsize=9,
             horizontalalignment='center', verticalalignment='bottom')

    plt.hlines(y=m_pl_des-200, xmin=0, xmax=2150, colors='pink', linestyle='--', linewidth=1)
    plt.text(2200 / 2, m_pl_des-200, f'Only Passenger Payload = {m_pl_des-200:.0f} kg', color='pink', fontsize=9,
             horizontalalignment='center', verticalalignment='bottom')
    #
    plt.vlines(x=R_des/1000, ymin=0, ymax=m_pl_des, colors='red', linestyle='--', linewidth=1)
    plt.text(R_des/1000, m_pl_des/2, f'Design Range = {R_des/1000:.0f} km', rotation=90, color='red', fontsize=9,
             horizontalalignment='right', verticalalignment='center')

    plt.hlines(y=payloads[3], xmin=0, xmax=ranges[3], colors='blue', linestyle='--', linewidth=1)
    plt.text(ranges[3]/2, m_pilot,f'Practical minimum payload = {m_pilot:.0f} kg', color='blue',fontsize=9,
             horizontalalignment='center', verticalalignment='bottom')


if __name__ == '__main__':
    # WEIGHTS
    m_mto = 1871
    m_e = 1035
    m_pl_des = 662
    m_pl_max = 700
    m_pilot = 77
    m_f_des = m_mto - m_e - m_pl_des
    m_f_max = 250

    # AERODYNAMICS
    A = 9
    e = 0.783
    C_d0 = 0.02591

    # PROPULSION
    eta_p = 0.83
    eta_eng = 0.25
    e_f = 44400000

    # DESIGN CHOICES
    R_des = 500 * 1000
    V_cr = 132 * KTS_TO_MS

    payload_range_points = calculate_payload_range_parameters(m_mto, m_e, m_f_des, m_f_max, m_pl_max, m_pilot,
                                                              A, e, C_d0, eta_p, eta_eng, e_f, R_des, V_cr)

    print(payload_range_points)
    plt.figure(figsize=(10,6))
    plot_payload_range(payload_range_points, m_pl_des, R_des)

    # m_mto = 1871 - 2 * 77
    # m_e = 1034
    # m_pl_des = 662 - 2 * 77
    # m_pl_max = 700 - 2 * 77
    # m_pilot = 77
    # m_f_des = m_mto - m_e - m_pl_des
    # m_f_max = 250
    #
    # payload_range_points = calculate_payload_range_parameters(m_mto, m_e, m_f_des, m_f_max, m_pl_max, m_pilot,
    #                                                           A, e, C_d0, eta_p, eta_eng, e_f, R_des, V_cr)
    #
    # print(payload_range_points)
    # plot_payload_range(payload_range_points, m_pl_des, R_des)

    plt.xlim(0, 3000)
    plt.ylim(0, 800)

    plt.xlabel("Range [km]")
    plt.ylabel("Payload [kg]")
    plt.title("Payload-Range Diagram")

    plt.grid(True)
    # plt.legend(fontsize=8)

    save_path = Path(__file__).parent / "performance_figures/payload_range_diagram.png"
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to: {save_path}")
    plt.close()
