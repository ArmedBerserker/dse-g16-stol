from classes.aircraft_2 import Aircraft, loader
from lookups.consts import *
import math
import matplotlib.pyplot as plt
import os
from class1.prelim_drag import *
from pathlib import Path


def calculate_payload_range_parameters(m_pl, m_pl_max, R_design, eta_eng, eta_p, A, e, C_d0, h_cr, V_cr, e_f, me_mto):

    cl_cd = 0.5 * math.sqrt(math.pi * A * e / C_d0)
    print(f"ClCd: {cl_cd}")
    # cl_cd = 11.9
    R_lost = 1 / 0.7 * cl_cd * (h_cr + V_cr ** 2 / (2 * g))

    t_e = 45 * 60
    R_eq_res = t_e * V_cr

    R_eq = R_design + R_lost + R_eq_res
    mf_mto = 1 - math.exp((-R_eq)/(eta_eng * eta_p * e_f * cl_cd / g))
    m_mto = m_pl / (1 - me_mto - mf_mto)

    m_e = me_mto * m_mto
    m_f = mf_mto * m_mto

    print(f"Fuel mass fraction: {mf_mto}")
    print(f"Empty mass: {m_e}")
    print(f"Fuel mass: {m_f}")
    print(f"Payload mass: {m_pl}")
    print(f"Maximum take-off mass: {m_mto}")

    design_point_1 = (0, m_pl_max)
    design_point_3 = (R_design, m_pl)

    R_aux = (R_eq - R_design)
    print(f"Equivalent range: {R_eq}")
    print(f"Auxiliary range: {R_aux}")

    m_f_at_struct_payload = m_mto - m_e - m_pl_max
    R_at_struct_payload = eta_eng * eta_p * cl_cd * e_f / g * math.log((m_e + m_pl_max + m_f_at_struct_payload) / (m_e + m_pl_max)) - R_aux
    design_point_2 = (R_at_struct_payload, m_pl_max)

    R_ferry = eta_eng * eta_p * cl_cd * e_f / g * math.log((m_e + m_f) / m_e) - R_aux
    design_point_4 = (R_ferry, 0)

    payload_range_points = [design_point_1, design_point_2, design_point_3, design_point_4]

    return payload_range_points


def plot_payload_range(payload_range_points,
                       save_path = Path(__file__).parent / "payload_range_diagram.png",
                       figsize = (10, 6)):

    ranges = [point[0] / 1000 for point in payload_range_points]
    payloads = [point[1] for point in payload_range_points]

    plt.figure(figsize=figsize)

    plt.plot(ranges, payloads, marker='o', linestyle='-')

    labels = [chr(65 + i) for i in range(len(payload_range_points))]
    for i, label in enumerate(labels):
        plt.text(ranges[i], payloads[i], f'  {label}', fontsize=10)

    plt.xlabel("Range")
    plt.ylabel("Payload")
    plt.title("Payload-Range Diagram")

    plt.grid(True)

    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')

    print(f"Figure saved to: {save_path}")


if __name__ == '__main__':
    file_path = '../yamls/aircraft.yaml'
    target_class = Aircraft
    aircraft = loader.load(file_path, target_class)

    m_pl = 704
    m_pl_max = 725  # m_pl * 1.03

    R_design = 500 * 1000

    eta_eng = 0.2 # 0.25 piston / 0.2 turboprop
    eta_p = 0.8
    A = 9
    e = 0.7996
    C_d0 = 0.04
    h_cr = 2590
    V_cr = 68
    e_f = 43000000
    me_mto = 0.497 # 0.62 for boosted piston / 0.52 for boosted turboprop - 0.02 added for supercapacitor

    payload_range_points = calculate_payload_range_parameters(m_pl, m_pl_max, R_design, eta_eng, eta_p, A, e, C_d0, h_cr, V_cr, e_f, me_mto)
    plot_payload_range(payload_range_points)


# def calculate_payload_range_parameters(mtom,
#                                        payload_weight,
#                                        eta_eng,
#                                        eta_p,
#                                        e_f,
#                                        oem,
#                                        AR,
#                                        e,
#                                        C_d0,
#                                        mfw):
#     """
#     design range
#     ferry range
#     range at max structural payload
#     """
#
#
#
#
#     maximum_structural_payload = 1.25 * payload_weight
#     Cl_Cd =  0.5 * math.sqrt((math.pi * AR * e)/C_d0)
#     print(Cl_Cd)
#
#     maximum_range_at_maximum_payload = eta_eng * eta_p * (e_f / g) * Cl_Cd * math.log(mtom / (oem + maximum_structural_payload)) - 288000
#     payload_weight = mtom - oem - mfw
#     maximum_range_at_maximum_fuel = eta_eng * eta_p * (e_f / g) * Cl_Cd * math.log(mtom / (oem + payload_weight)) - 288000
#     ferry_range = eta_eng * eta_p * (e_f / g) * Cl_Cd *  math.log((oem + mfw) / oem) - 288000
#
#     payload_range_points = [(0, maximum_structural_payload), (maximum_range_at_maximum_payload, maximum_structural_payload),
#                             (maximum_range_at_maximum_fuel, payload_weight), (ferry_range, 0)]
#     print(payload_range_points)
#     return payload_range_points