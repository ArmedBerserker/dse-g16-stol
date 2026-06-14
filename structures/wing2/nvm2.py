from loads import *
from geom import *
import numpy as np
import matplotlib.pyplot as plt


def compute_wingbox_material_area(nodes, elements, t_skin, t_web, A_spar, A_stringer, meta):
    # 1. Discrete Boom Components (Spar caps + Stringers)
    spar_nodes = {
        meta["idx_front_spar_bottom"], meta["idx_front_spar_top"],
        meta["idx_main_spar_bottom"], meta["idx_main_spar_top"]
    }

    total_boom_area = 0.0
    for idx in range(len(nodes)):
        total_boom_area += A_spar if idx in spar_nodes else A_stringer

    # 2. Continuous Internal Webs Components
    h_front = meta["front_spar_height"]
    h_main = meta["main_spar_height"]
    total_web_area = (h_front * t_web) + (h_main * t_web)

    # 3. True Continuous External Skin Component
    total_skin_perimeter = meta["airfoil_perimeter"]
    total_skin_area = total_skin_perimeter * t_skin

    # Global cross-sectional sum (m^2)
    total_structural_area = total_boom_area + total_web_area + total_skin_area
    return total_structural_area


def calculate_and_plot_vmt(
        # --- Structural/Aero Inputs ---
        half_wing_lift_target, b, y_root, c, q,
        flap_start, flap_end, delta_cl_flap,
        aileron_start, aileron_end, delta_cl_aileron,
        filepath, front_spar_pct, main_spar_pct, bays,

        # --- Material & Point Load Constants ---
        rho_material=2780.0,  # Aluminum 7075-T6 density (kg/m^3)
        engine_mass=92.0,  # Engine weight in kg
        engine_y_loc=2.5,  # Engine spanwise position in meters
        engine_x_loc=0.35,  # Engine CG relative to Leading Edge (meters)
        wtd_mass=15.0,  # Wing Tip Device weight in kg
        wtd_x_loc=0.50,  # WTD CG relative to Leading Edge (meters)
        n_load_factor=1,
        plot=True,

        # --- New Mass Items ---
        hld_actuators=None,      # List of dicts: [{'y': float, 'mass': float, 'x_loc': float}]
        spoiler_actuators=None,  # List of dicts: [{'y': float, 'mass': float, 'x_loc': float}]
        fuel_tank_bays=None      # List of dicts: [{'y_start': float, 'y_end': float, 'mass_empty': float}]
):
    g = 9.81  # m/s^2

    # Avoid mutable default argument issues
    if hld_actuators is None: hld_actuators = []
    if spoiler_actuators is None: spoiler_actuators = []
    if fuel_tank_bays is None: fuel_tank_bays = []

    # 1. Generate aerodynamic lift distribution along the span
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

    # 2. Map coordinates relative to the bottom-left spar cap origin
    x_front_spar = front_spar_pct * c

    # Moment arms relative to the front spar origin
    x_aero_spar = (0.25 * c) - x_front_spar
    engine_x_spar = engine_x_loc - x_front_spar
    wtd_x_spar = wtd_x_loc - x_front_spar

    # Map secondary actuator positions relative to the front spar origin
    for act in hld_actuators + spoiler_actuators:
        act['x_spar'] = act['x_loc'] - x_front_spar
        act['weight'] = act['mass'] * g * n_load_factor

    # 3. Dynamic Structural Geometry & Mass Allocation per Bay
    num_pts = len(y_stations)
    struct_weight_dist = np.zeros(num_pts)   # Isolated pure structural weight (N/m)
    inertial_weight_dist = np.zeros(num_pts) # Isolated non-structural distributed weight (N/m)
    web_weight_dist = np.zeros(num_pts)      # Track web-only weight distribution
    X_bar_array = np.zeros(num_pts)

    y_prev = y_root - 1e-6
    for bay in bays:
        mask = (y_stations > y_prev) & (y_stations <= bay['y_end'])
        if not np.any(mask):
            continue

        # Generate idealization once per segment/bay
        nodes, elements, meta, _ = generate_megson_idealization(
            filepath, c, front_spar_pct, main_spar_pct,
            bay['num_box_stringers'], bay['num_le_booms']
        )

        # X_bar is inherently calculated relative to the bottom-left spar cap origin
        X_bar, _, _, _, _, _, _, _, _, _, _, _, _ = analyze_section_properties(
            nodes, elements, bay['t_skin'], bay['t_web'],
            bay['A_spar'], bay['A_stringer'], meta
        )

        total_mat_area = compute_wingbox_material_area(
            nodes, elements, bay['t_skin'], bay['t_web'],
            bay['A_spar'], bay['A_stringer'], meta
        )

        # Separate calculation for internal webs to calculate lightening holes later
        h_front = meta["front_spar_height"]
        h_main = meta["main_spar_height"]
        total_web_area = (h_front * bay['t_web']) + (h_main * bay['t_web'])

        struct_weight_dist[mask] = total_mat_area * rho_material * g
        web_weight_dist[mask] = total_web_area * rho_material * g
        X_bar_array[mask] = X_bar

        for tank in fuel_tank_bays:
            tank_mask = mask & (y_stations >= tank['y_start']) & (y_stations <= tank['y_end'])
            if np.any(tank_mask):
                tank_length = tank['y_end'] - tank['y_start']
                w_tank_per_meter = (tank['mass_empty'] / tank_length) * g
                inertial_weight_dist[tank_mask] += w_tank_per_meter

        y_prev = bay['y_end']

    raw_wing_mass = np.trapezoid(struct_weight_dist / g, y_stations) * 2
    full_web_mass = np.trapezoid(web_weight_dist / g, y_stations) * 2

    # Apply weight savings from lightening cutouts
    web_hole_savings = full_web_mass * 0.50
    net_sim_mass = raw_wing_mass * 0.8 - web_hole_savings

    # --- Rib Mass Calculation ---
    y_ribs = [y_root] + [bay['y_end'] for bay in bays]
    y_tip = y_ribs[-1]
    total_ribs_mass = 0.0

    for station_y in y_ribs:
        span_fraction = (station_y - y_root) / (y_tip - y_root) if (y_tip - y_root) > 0 else 0
        local_rib_thickness = max(0.0005, 0.0010 * (1.0 - 0.5 * span_fraction))
        rib_lightening_factor = 1.0 - (0.35 + 0.25 * span_fraction)
        approx_rib_area = c * (c * 0.12)  # Refers to local chord 'c' and standard chord thickness

        # Multiply by 2 to account for both left and right wing sides
        total_ribs_mass += approx_rib_area * local_rib_thickness * rho_material * rib_lightening_factor * 2

    true_production_mass = net_sim_mass + total_ribs_mass

    # --- Load Vector Recombination for downstream VMT Integration ---
    # Combine structural and non-structural components here so VMT gets everything
    weight_dist = struct_weight_dist + inertial_weight_dist
    w_net = (lift_dist + weight_dist)

    # Distributed torsion about the front spar origin
    m_spar_dist = (-lift_dist * x_aero_spar)

    # 4. Integrate Shear, Bending Moment, and Torque from Tip to Root
    V = np.zeros(num_pts)
    M = np.zeros(num_pts)
    T_origin = np.zeros(num_pts)

    W_engine = engine_mass * g
    W_wtd = wtd_mass * g

    # Boundary conditions at the free wingtip (i = last element)
    V[-1] = W_wtd
    M[-1] = 0.0
    T_origin[-1] = -W_wtd * wtd_x_spar

    for i in range(num_pts - 2, -1, -1):
        dy = y_stations[i + 1] - y_stations[i]

        w_avg = 0.5 * (w_net[i] + w_net[i + 1])
        m_avg = 0.5 * (m_spar_dist[i] + m_spar_dist[i + 1])

        V[i] = V[i + 1] + w_avg * dy
        M[i] = M[i + 1] - V[i + 1] * dy - 0.5 * w_avg * dy ** 2
        T_origin[i] = T_origin[i + 1] + m_avg * dy

        # Intercept and apply engine point load and its torque about origin
        if y_stations[i] <= engine_y_loc < y_stations[i + 1]:
            V[i] += W_engine
            T_origin[i] += W_engine * engine_x_spar

        # Intercept and apply discrete HLD actuator point weights
        for act in hld_actuators:
            if y_stations[i] <= act['y'] < y_stations[i + 1]:
                V[i] += act['weight']
                T_origin[i] += act['weight'] * act['x_spar']

        # Intercept and apply discrete Spoiler actuator point weights
        for act in spoiler_actuators:
            if y_stations[i] <= act['y'] < y_stations[i + 1]:
                V[i] += act['weight']
                T_origin[i] += act['weight'] * act['x_spar']

    # 5. Shift reference axis from the Front Spar Origin to the Elastic Axis (X_bar)
    T_ea = T_origin

    # 6. Plotting
    if plot:
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

    return y_stations, V, M, T_ea, raw_wing_mass, true_production_mass, total_ribs_mass, web_hole_savings

if __name__ == "__main__":
    # Define External Weights & Point Placements
    engine_spanwise_pos = 2.2
    engine_mass_kg = 122.0
    wtd_mass_kg = 40.0

    wing_bays = [
        {'y_end': 2.4, 'num_box_stringers': 14, 'num_le_booms': 5,
         't_skin': 0.0005, 't_web': 0.001, 'A_spar': 7.5E-4, 'A_stringer': 7E-5},

        {'y_end': 4.0, 'num_box_stringers': 10, 'num_le_booms': 4,
         't_skin': 0.0005, 't_web': 0.001, 'A_spar': 4.5E-4, 'A_stringer': 7E-5},

        {'y_end': 5.8, 'num_box_stringers': 6, 'num_le_booms': 2,
         't_skin': 0.0005, 't_web': 0.001, 'A_spar': 2.0E-4, 'A_stringer': 7E-5},

        {'y_end': 7.4, 'num_box_stringers': 3, 'num_le_booms': 1,
         't_skin': 0.0005, 't_web': 0.001, 'A_spar': 0.6E-4, 'A_stringer': 7E-5},

        {'y_end': 8.4, 'num_box_stringers': 1, 'num_le_booms': 1,
         't_skin': 0.0005, 't_web': 0.001, 'A_spar': 0.1E-4, 'A_stringer': 7E-5}
    ]

    y_steps, shear, moment, torsion, raw_mass, true_mass, rib_mass, savings = calculate_and_plot_vmt(
        # Aerodynamic Arguments
        half_wing_lift_target=10000,
        b=16.8,
        y_root=0.725,
        c=1.87,
        q=571.0,
        flap_start=1.1,
        flap_end=3.8,
        delta_cl_flap=1.25,
        aileron_start=5.6,
        aileron_end=7.5,
        delta_cl_aileron=0.9,

        # Structural Layout Arguments
        filepath="NACA23012.dat",
        front_spar_pct=0.16,
        main_spar_pct=0.56,
        bays=wing_bays,

        # Auxiliary Point Load Configurations
        engine_mass=engine_mass_kg,
        engine_y_loc=engine_spanwise_pos,
        engine_x_loc=-0.1 * 1.87,
        wtd_mass=wtd_mass_kg,
        wtd_x_loc=0.25 * 1.87,
        n_load_factor=1
    )

    print("\n=================== MASS BREAKDOWN ===================")
    print(f"Raw Simulation Weight (Solid Webs, No Ribs): {raw_mass:.2f} kg")
    print(f"Web Lightening Hole Savings (50% reduction): -{savings:.2f} kg")
    print(f"Total Added Ribs Mass:                       +{rib_mass:.2f} kg")
    print(f"TRUE PRODUCTION MASS (Full Wing):            {true_mass:.2f} kg")
    print("======================================================")
    print(f"Root Shear Force:   {shear[0] / 1000:.2f} kN")
    print(f"Root Bending Moment:{moment[0] / 1000:.2f} kN*m")
    print(f"Root Torsion:       {torsion[0]:.2f} N*m")