import math
import matplotlib.pyplot as plt
import os
from pathlib import Path
from lookups.consts import *


def calculate_payload_range_parameters(m_e, m_pl, m_mto, A, e, C_d0, h_cr, V_cr, R_design, eta_eng, eta_p, e_f, m_max_fuel,
                                    m_pl_max, me_mto, name):


    m_f_given = m_mto - m_e - m_pl

    cl_cd = 0.5 * math.sqrt(math.pi * A * e / C_d0)
    R_lost = 1 / 0.7 * cl_cd * (h_cr + V_cr ** 2 / (2 * g))

    t_e = 45 * 60
    R_eq_res = t_e * V_cr

    R_eq = R_design + R_lost + R_eq_res

    m_f_required = m_mto * (1 - math.exp((-R_eq)/(eta_eng * eta_p * e_f * cl_cd / g)))

    if m_f_given < m_f_required:
        print(f"Fuel mass ({m_f_given}) for empty mass ({m_e}) is lower than required fuel mass ({m_f_required})")
        return None
    elif m_f_given > m_max_fuel:
        print(f"Fuel mass ({m_f_given}) for empty mass ({m_e}) is higher than fuel tank capacity ({m_max_fuel}) \n Limiting fuel mass to maximum fuel tank volume.")
        m_f = m_max_fuel

    mf_mto = 1 - math.exp((-R_eq)/(eta_eng * eta_p * e_f * cl_cd / g))
    m_mto = m_pl / (1 - me_mto - mf_mto)

    m_e = me_mto * m_mto
    m_f = mf_mto * m_mto

    print(f"Fuel mass for {name}: {m_f:.2f} kg")
    print(f"Empty mass for {name}: {m_e:.2f} kg")
    print(f"Maximum take-off mass for {name}: {m_mto:.2f} kg")

    design_point_1 = (0, m_pl_max)
    design_point_3 = (R_design, m_pl)

    R_aux = (R_eq - R_design)

    m_f_at_struct_payload = m_mto - m_e - m_pl_max
    R_at_struct_payload = eta_eng * eta_p * cl_cd * e_f / g * math.log((m_e + m_pl_max + m_f_at_struct_payload) / (m_e + m_pl_max)) - R_aux
    design_point_2 = (R_at_struct_payload, m_pl_max)

    R_ferry = eta_eng * eta_p * cl_cd * e_f / g * math.log((m_e + m_f) / m_e) - R_aux
    design_point_4 = (R_ferry, 0)

    payload_range_points = [design_point_1, design_point_2, design_point_3, design_point_4]

    return payload_range_points


if __name__ == '__main__':
    m_pl = 704
    m_pl_max = 725
    m_mto = 1870 # TO BE FILLED IN
    m_max_fuel = 213
    pilot_weight = 84

    A = 9 # TO BE FILLED IN
    e = 0.8 # TO BE FILLED IN

    eta_p = 0.8
    e_f = 43000000

    R_design = 500 * 1000 # TO BE FILLED IN
    h_cr = 8500 * FT_TO_M # TO BE FILLED IN
    V_cr = 132 * KTS_TO_MS # TO BE FILLED IN

    # Boosted piston taildragger
    eta_eng = 0.25
    C_d0 = 0.04 # TO BE FILLED IN
    m_e = 1000 # TO BE FILLED IN

    payload_range_points_boosted_piston_taildragger = calculate_payload_range_parameters(m_e, m_pl, m_mto, A, e, C_d0, h_cr, V_cr,
                                                                                         R_design, eta_eng, eta_p, e_f, m_max_fuel,

                                                                                         m_pl, m_pl_max, R_design,
                                                                                         eta_eng, eta_p, A, e, C_d0,
                                                                                         h_cr, V_cr, e_f, me_mto,
                                                                                         "boosted piston taildragger")
    #
    # # Boosted turboprop taildragger
    # eta_eng = 0.2
    # C_d0 = 0.04 # 0.0514 # 0.04  # TBC
    # me_mto = 0.585 #0.502
    # payload_range_points_boosted_turboprop_taildragger = calculate_payload_range_parameters(m_pl, m_pl_max, R_design,
    #                                                                                      eta_eng, eta_p, A, e, C_d0,
    #                                                                                      h_cr, V_cr, e_f, me_mto,
    #                                                                                     "boosted turboprop taildragger")
    #
    # # Boosted piston tricycle
    # eta_eng = 0.25
    # C_d0 = 0.04 # 0.0953 # 0.04  # TBC
    # me_mto = 0.590 #0.611
    # payload_range_points_boosted_piston_tricycle = calculate_payload_range_parameters(m_pl, m_pl_max, R_design,
    #                                                                                         eta_eng, eta_p, A, e, C_d0,
    #                                                                                         h_cr, V_cr, e_f, me_mto,
    #                                                                                   "boosted piston tricycle")
    #
    # # Boosted turboprop tricycle
    # eta_eng = 0.2
    # C_d0 = 0.04 # 0.0953 # 0.04  # TBC
    # me_mto = 0.584 #0.526
    # payload_range_points_boosted_turboprop_tricycle = calculate_payload_range_parameters(m_pl, m_pl_max, R_design,
    #                                                                                         eta_eng, eta_p, A, e, C_d0,
    #                                                                                         h_cr, V_cr, e_f, me_mto,
    #                                                                                 "boosted turboprop tricycle")
    # plot_payload_range(
    #     [
    #         payload_range_points_boosted_piston_taildragger,
    #         payload_range_points_boosted_turboprop_taildragger,
    #         payload_range_points_boosted_piston_tricycle,
    #         payload_range_points_boosted_turboprop_tricycle
    #     ],
    #     labels=[
    #         "Boosted Piston Taildragger",
    #         "Boosted Turboprop Taildragger",
    #         "Boosted Piston Tricycle",
    #         "Boosted Turboprop Tricycle"
    #     ],
    #     design_payload=m_pl,
    #     design_range=R_design,
    #     practical_min_payload=pilot_weight
    # )
