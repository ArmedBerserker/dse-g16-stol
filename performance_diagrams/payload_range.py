import math
import matplotlib.pyplot as plt
import os
from pathlib import Path
from lookups.consts import *


def calculate_payload_range_parameters(m_pl, m_pl_max, R_design, eta_eng, eta_p, A, e, C_d0, h_cr, V_cr, e_f, me_mto, name):

    cl_cd = 0.5 * math.sqrt(math.pi * A * e / C_d0)
    R_lost = 1 / 0.7 * cl_cd * (h_cr + V_cr ** 2 / (2 * g))
    t_e = 45 * 60
    R_eq_res = t_e * V_cr

    R_eq = R_design + R_lost + R_eq_res
    mf_mto = 1 - math.exp((-R_eq)/(eta_eng * eta_p * e_f * cl_cd / g))
    m_mto = m_pl / (1 - me_mto - mf_mto)

    m_e = me_mto * m_mto
    m_f = mf_mto * m_mto

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


def plot_payload_range(payload_range_datasets,
                       labels,
                       design_payload,
                       design_range,
                       practical_min_payload,
                       save_path=Path(__file__).parent / "payload_range_diagram.png",
                       figsize=(10, 6)):

    plt.figure(figsize=figsize)

    for idx, payload_range_points in enumerate(payload_range_datasets):

        ranges = [point[0] / 1000 for point in payload_range_points]
        payloads = [point[1] for point in payload_range_points]

        plt.plot(
            ranges,
            payloads,
            marker='o',
            linestyle='-',
            linewidth=2,
            label=labels[idx],
            zorder=3
        )

        # Point labels
        point_labels = ['A', 'B', 'C', 'D']

        for i, point_label in enumerate(point_labels):
            if point_label == 'B':
                plt.text(
                    386,
                    735,
                    f'{point_label}',
                    fontsize=8
                )
            else:
                plt.text(
                    ranges[i],
                    payloads[i] + 10,
                    f' {point_label}',
                    fontsize=8
                )

        for i in range(len(payload_range_points) - 1):

            x1, y1 = payload_range_points[i]
            x2, y2 = payload_range_points[i + 1]

            # Check if segment crosses practical minimum payload
            if (y1 - practical_min_payload) * (y2 - practical_min_payload) <= 0 and y1 != y2:

                # Linear interpolation
                t = (practical_min_payload - y1) / (y2 - y1)
                x_intersection = x1 + t * (x2 - x1)

                x_intersection_km = x_intersection / 1000

                # Vertical dashed line
                plt.axvline(
                    x=x_intersection_km,
                    color = 'blue',
                    linestyle='--',
                    linewidth=1
                )


                plt.text(
                    x_intersection_km - 20,
                    600,
                    f'{x_intersection_km:.0f} km',
                    rotation=90,
                    color='blue',
                    fontsize=8,
                    va='bottom',
                    zorder=4
                )

                break


    plt.axhline(
        y=design_payload,
        color='red',
        linestyle='--',
        linewidth=1,
        zorder=0,
    )

    plt.text(
        100,
        design_payload - 20,
        f'Design Payload = {design_payload:.0f} kg',
        color='red',
        fontsize=8
    )

    plt.axvline(
        x=design_range / 1000,
        color='red',
        linestyle='--',
        linewidth=1,
        zorder=0,
    )

    plt.text(
        design_range / 1000 - 20,
        200,
        f'Design Range = {design_range / 1000:.0f} km',
        rotation=90,
        color='red',
        fontsize=8,
        va='bottom'
    )

    plt.axhline(
        y=practical_min_payload,
        color='blue',
        linestyle='--',
        linewidth=1,
        zorder=0,
    )

    plt.text(
        100,
        practical_min_payload + 5,
        f'Practical minimum payload = {practical_min_payload:.0f} kg',
        color='blue',
        fontsize=8
    )

    bbox = dict(boxstyle="round", fc="0.8")
    arrowprops = dict(
        arrowstyle="->",
        # connectionstyle="angle,angleA=0,angleB=90,rad=10"
    )

    plt.annotate(
        f"Limited by maximum structural payload",
        (200, 725),
        xytext=(27, 800),
        bbox=bbox, arrowprops=arrowprops, fontsize=8)

    plt.annotate(
        f"Payload mass traded for fuel mass",
        (450, 712),
        xytext=(450, 800),
        bbox=bbox, arrowprops=arrowprops, fontsize=8)

    plt.text(
        550,
        200,
        f"Range limited by fuel tank volume",
        None,
        bbox=bbox, fontsize=8, rotation=-57)

    plt.xlim(0, 1200)
    plt.ylim(0, 900)

    plt.xlabel("Range [km]")
    plt.ylabel("Payload [kg]")
    plt.title("Payload-Range Diagram")

    plt.grid(True)
    plt.legend(fontsize=8)

    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')

    print(f"Figure saved to: {save_path}")


if __name__ == '__main__':
    m_pl = 704
    m_pl_max = 725  # m_pl * 1.03
    R_design = 500 * 1000
    pilot_weight = 84
    eta_p = 0.8
    A = 9
    e = 0.8
    h_cr = 2590.8
    V_cr = 67.9
    e_f = 43000000

    # Boosted piston taildragger
    eta_eng = 0.25
    C_d0 = 0.04 # TBC
    me_mto = 0.587
    payload_range_points_boosted_piston_taildragger = calculate_payload_range_parameters(m_pl, m_pl_max, R_design,
                                                                                         eta_eng, eta_p, A, e, C_d0,
                                                                                         h_cr, V_cr, e_f, me_mto,
                                                                                         "boosted piston taildragger")

    # Boosted turboprop taildragger
    eta_eng = 0.2
    C_d0 = 0.04  # TBC
    me_mto = 0.502
    payload_range_points_boosted_turboprop_taildragger = calculate_payload_range_parameters(m_pl, m_pl_max, R_design,
                                                                                         eta_eng, eta_p, A, e, C_d0,
                                                                                         h_cr, V_cr, e_f, me_mto,
                                                                                        "boosted turboprop taildragger")

    # Boosted piston tricycle
    eta_eng = 0.25
    C_d0 = 0.04  # TBC
    me_mto = 0.611
    payload_range_points_boosted_piston_tricycle = calculate_payload_range_parameters(m_pl, m_pl_max, R_design,
                                                                                            eta_eng, eta_p, A, e, C_d0,
                                                                                            h_cr, V_cr, e_f, me_mto,
                                                                                      "boosted piston tricycle")

    # Boosted turboprop tricycle
    eta_eng = 0.2
    C_d0 = 0.04  # TBC
    me_mto = 0.526
    payload_range_points_boosted_turboprop_tricycle = calculate_payload_range_parameters(m_pl, m_pl_max, R_design,
                                                                                            eta_eng, eta_p, A, e, C_d0,
                                                                                            h_cr, V_cr, e_f, me_mto,
                                                                                    "boosted turboprop tricycle")
    plot_payload_range(
        [
            payload_range_points_boosted_piston_taildragger,
            payload_range_points_boosted_turboprop_taildragger,
            payload_range_points_boosted_piston_tricycle,
            payload_range_points_boosted_turboprop_tricycle
        ],
        labels=[
            "Boosted Piston Taildragger",
            "Boosted Turboprop Taildragger",
            "Boosted Piston Tricycle",
            "Boosted Turboprop Tricycle"
        ],
        design_payload=m_pl,
        design_range=R_design,
        practical_min_payload=pilot_weight
    )
