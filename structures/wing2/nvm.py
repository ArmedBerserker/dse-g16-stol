from loads import *
from geom import *

def calculate_and_plot_vmt(
        # --- Structural/Aero Inputs ---
        half_wing_lift_target, b, y_root, c, q,
        flap_start, flap_end, delta_cl_flap,
        aileron_start, aileron_end, delta_cl_aileron,
        filepath, front_spar_pct, main_spar_pct, num_box_stringers, num_le_booms,
        t_skin, t_web, A_spar, A_stringer,

        # --- Material & Point Load Constants ---
        rho_material=2780.0,  # Aluminum 7075-T6 density (kg/m^3)
        engine_mass = 92.0,  # Engine weight in kg
        engine_y_loc=2.5,  # Engine spanwise position in meters
        engine_x_loc=0.35,  # Engine CG relative to Leading Edge (meters)
        wtd_mass=15.0,  # Wing Tip Device weight in kg
        wtd_x_loc=0.50,  # WTD CG relative to Leading Edge (meters)
        n_load_factor=1,
        plot = True # Load Factor (e.g., Limit Load Factor)
):
    g = 9.81  # m/s^2

    y_stations, lift_dist = schrenk_half_wing_loading(
        half_wing_lift_target, b, y_root, c, q,
        flap_start, flap_end, delta_cl_flap,
        aileron_start, aileron_end, delta_cl_aileron,
        num_points=500
    )

    # Clean data sort order (root to tip)
    idx = np.argsort(y_stations)
    y_stations = y_stations[idx]
    lift_dist = -lift_dist[idx] * n_load_factor
    x_aero = 0.25 * c  # Aerodynamic Center

    nodes, elements, meta, _ = generate_megson_idealization(
        filepath, c, front_spar_pct, main_spar_pct, num_box_stringers, num_le_booms
    )

    X_bar, _, total_mat_area, _, _, _ = analyze_section_properties(
        nodes, elements, t_skin, t_web, A_spar, A_stringer, meta
    )

    # Distributed wing structural weight per meter (N/m)
    weight_dist = -total_mat_area * rho_material * g
    x_cg_dist = X_bar  # Structural centroid

    # Net loading vectors per meter
    w_net = (lift_dist - weight_dist)
    m_le_dist = (lift_dist * x_aero) - (weight_dist * x_cg_dist)

    num_pts = len(y_stations)
    V = np.zeros(num_pts)
    M = np.zeros(num_pts)
    T_le = np.zeros(num_pts)

    W_engine = -engine_mass * g
    W_wtd = -wtd_mass * g

    # Boundary conditions at the free wingtip (i = last element)
    V[-1] = -W_wtd
    M[-1] = 0.0
    T_le[-1] = -W_wtd * wtd_x_loc

    for i in range(num_pts - 2, -1, -1):
        dy = y_stations[i + 1] - y_stations[i]

        w_avg = 0.5 * (w_net[i] + w_net[i + 1])
        m_avg = 0.5 * (m_le_dist[i] + m_le_dist[i + 1])

        V[i] = V[i + 1] + w_avg * dy
        M[i] = M[i + 1]  - V[i + 1] * dy - 0.5 * w_avg * dy ** 2
        T_le[i] = T_le[i + 1] + m_avg * dy

        # Intercept and append engine point load
        if y_stations[i] <= engine_y_loc < y_stations[i + 1]:
            V[i] += W_engine
            T_le[i] += W_engine * engine_x_loc

    # Shift reference axis from Leading Edge to the Elastic Axis (X_bar)
    T_ea = T_le - V * X_bar

    if plot == True:
        fig, axs = plt.subplots(3, 1, figsize=(10, 11), sharex=True)

        # Shear Diagram
        axs[0].plot(y_stations, V / 1000.0, 'b-', linewidth=2)
        axs[0].set_ylabel('Shear Force [kN]', fontweight='bold')
        axs[0].grid(True, linestyle=':')
        axs[0].axvline(x=engine_y_loc, color='orange', linestyle='--', label=f'Engine ({engine_y_loc}m)')
        axs[0].legend(loc='upper right')

        # Bending Moment Diagram
        axs[1].plot(y_stations, M / 1000.0, 'r-', linewidth=2)
        axs[1].set_ylabel('Bending Moment [kN·m]', fontweight='bold')
        axs[1].grid(True, linestyle=':')
        axs[1].axvline(x=engine_y_loc, color='orange', linestyle='--')

        # Torque Diagram
        axs[2].plot(y_stations, T_ea, 'g-', linewidth=2)
        axs[2].set_ylabel('Torsion [N·m]', fontweight='bold')
        axs[2].set_xlabel('Spanwise Coordinate y [m]', fontweight='bold')
        axs[2].grid(True, linestyle=':')
        axs[2].axvline(x=engine_y_loc, color='orange', linestyle='--')

        plt.tight_layout()
        plt.show()

    return y_stations, V, M, T_ea


if __name__ == "__main__":
    # Define Structural Cross-Section Configuration
    dat_file = "NACA23012.dat"
    chord_length = 1.87
    thickness_skin = 0.001
    thickness_web = 0.004
    area_spar_cap = 0.0005
    area_stringer = 0.0001

    # Define External Weights & Point Placements
    engine_spanwise_pos = 2.2 # Engine location along span (m)
    engine_mass_kg = 122.0  #  mass in kg
    wtd_mass_kg = 20.0  # Wing tip device mass in kg

    y_steps, shear, moment, torsion = calculate_and_plot_vmt(
        # Aerodynamic Arguments
        half_wing_lift_target=10000,
        b=16.8,
        y_root=0.725,
        c=chord_length,
        q=571.0,
        flap_start=1.1,
        flap_end=3.8,
        delta_cl_flap=1.25,
        aileron_start=5.6,
        aileron_end=7.5,
        delta_cl_aileron=0.9,

        # Structural Layout Arguments
        filepath=dat_file,
        front_spar_pct=0.16,
        main_spar_pct=0.56,
        num_box_stringers=8,
        num_le_booms=4,
        t_skin=thickness_skin,
        t_web=thickness_web,
        A_spar=area_spar_cap,
        A_stringer=area_stringer,

        # Auxiliary Point Load Configurations
        engine_mass=engine_mass_kg,
        engine_y_loc=engine_spanwise_pos,
        engine_x_loc=0.3 * chord_length,  # e.g., Engine CG at 30% local chord
        wtd_mass=wtd_mass_kg,
        wtd_x_loc=0.25 * chord_length,  # Wingtip device twist point
        n_load_factor=1
    )

    print(f"Root Shear Force:   {shear[0] / 1000:.2f} kN")
    print(f"Root Bending Moment:{moment[0] / 1000:.2f} kN*m")
    print(f"Root Torsion {torsion[0]:.2f} N*m")