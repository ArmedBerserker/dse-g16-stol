import numpy as np
from scipy.integrate import cumulative_trapezoid
import matplotlib.pyplot as plt

# Assuming geom.py and nvm2.py are in the same directory
from geom import generate_megson_idealization, analyze_section_properties
from nvm2 import calculate_and_plot_vmt


def calculate_section_properties_array(y_stations, bays, filepath, c, front_spar_pct, main_spar_pct):
    num_pts = len(y_stations)

    # Initialize the properties dictionary including the missing spar area arrays
    props = {
        'Ixx': np.zeros(num_pts), 'Iyy': np.zeros(num_pts), 'Ixy': np.zeros(num_pts),
        'Ae1': np.zeros(num_pts), 'Ae2': np.zeros(num_pts),
        'd1': np.zeros(num_pts), 'd2': np.zeros(num_pts), 'd12': np.zeros(num_pts),
        'X_bar': np.zeros(num_pts), 'Y_bar': np.zeros(num_pts),
        't_skin_arr': np.zeros(num_pts), 't_web_arr': np.zeros(num_pts),
        'A_frontspar': np.zeros(num_pts),
        'A_rearspar': np.zeros(num_pts),
        'nodes': [None] * num_pts,
        'elements': [None] * num_pts
    }

    y_prev = y_stations[0] - 1e-6
    for bay in bays:
        mask = (y_stations > y_prev) & (y_stations <= bay['y_end'])
        if not np.any(mask):
            continue

        nodes, elements, meta, _ = generate_megson_idealization(
            filepath, c, front_spar_pct, main_spar_pct,
            bay['num_box_stringers'], bay['num_le_booms']
        )

        (X_bar, Y_bar, total_boom_area, Ixx, Iyy, Ixy, Ae1, Ae2,
         d1, d2, d12, A_frontspar, A_rearspar) = analyze_section_properties(
            nodes, elements, bay['t_skin'], bay['t_web'],
            bay['A_spar'], bay['A_stringer'], meta
        )

        # Assign values to masks across the span
        props['Ixx'][mask] = Ixx
        props['Iyy'][mask] = Iyy
        props['Ixy'][mask] = Ixy
        props['Ae1'][mask] = Ae1
        props['Ae2'][mask] = Ae2
        props['d1'][mask] = d1
        props['d2'][mask] = d2
        props['d12'][mask] = d12
        props['X_bar'][mask] = X_bar
        props['Y_bar'][mask] = Y_bar
        props['t_skin_arr'][mask] = bay['t_skin']
        props['t_web_arr'][mask] = bay['t_web']
        props['A_frontspar'][mask] = A_frontspar
        props['A_rearspar'][mask] = A_rearspar

        # Rigorously map discrete nodes to spanwise indices
        for idx in range(num_pts):
            if mask[idx]:
                props['nodes'][idx] = nodes
                props['elements'][idx] = elements

        y_prev = bay['y_end']

    return props


def solve_multicell_torsion(T, d1, d2, d12, Ae1, Ae2, G):
    matrix = np.array([
        [d1, -d12, -2.0 * G * Ae1],
        [-d12, d2, -2.0 * G * Ae2],
        [2.0 * Ae1, 2.0 * Ae2, 0.0]
    ])
    rhs = np.array([0, 0, T])

    try:
        solution = np.linalg.solve(matrix, rhs)
        return solution[0], solution[1], solution[2]
    except np.linalg.LinAlgError:
        return 0.0, 0.0, 0.0


def calculate_deflections_and_stresses(y_stations, M, T_ea, props, E=71.7e9, G=26.9e9):
    num_pts = len(y_stations)

    # 1. Bending Deflection (Double Integration of Curvature)
    denominator = (props['Ixx'] * props['Iyy']) - (props['Ixy'] ** 2)
    curvature = (M * props['Iyy']) / (E * denominator)

    slope = cumulative_trapezoid(curvature, y_stations, initial=0.0)
    deflection = cumulative_trapezoid(slope, y_stations, initial=0.0)

    # 2. Torsional Deflection & Pure Torsion Shear Flows
    q1_arr = np.zeros(num_pts)
    q2_arr = np.zeros(num_pts)
    dtheta_dy = np.zeros(num_pts)

    for i in range(num_pts):
        q1, q2, dth = solve_multicell_torsion(
            T_ea[i], props['d1'][i], props['d2'][i], props['d12'][i],
            props['Ae1'][i], props['Ae2'][i], G
        )
        q1_arr[i] = q1
        q2_arr[i] = q2
        dtheta_dy[i] = dth

    twist_angle = cumulative_trapezoid(dtheta_dy, y_stations, initial=0.0)

    # 3. Rigorous Tabular Stress Calculation (Iterating every boom per station)
    max_bending_stress = np.zeros(num_pts)
    max_shear_stress = np.zeros(num_pts)

    for i in range(num_pts):
        # A. Bending Stress (Iterating through the idealized nodes array)
        M_x = M[i]
        X_bar = props['X_bar'][i]
        Y_bar = props['Y_bar'][i]
        denom_i = denominator[i]

        max_sigma_local = 0.0

        for node in props['nodes'][i]:
            # Robustly handle dict or list node representations from geom.py
            x_raw = node['x'] if isinstance(node, dict) else node[0]
            y_raw = node['y'] if isinstance(node, dict) else node[1]

            # Translate coordinates relative to the neutral axis
            x_prime = x_raw - X_bar
            y_prime = y_raw - Y_bar

            # Full Generalized Unsymmetrical Bending Equation (Outputs in Pascals)
            sigma_z = (M_x * (props['Iyy'][i] * y_prime - props['Ixy'][i] * x_prime)) / denom_i

            if abs(sigma_z) > max_sigma_local:
                max_sigma_local = abs(sigma_z)

        max_bending_stress[i] = max_sigma_local

        # B. Shear Stresses from Torsion (tau = q / t)
        t_skin = props['t_skin_arr'][i]
        t_web = props['t_web_arr'][i]

        tau_cell1 = np.abs(q1_arr[i] / t_skin)
        tau_cell2 = np.abs(q2_arr[i] / t_skin)
        tau_web = np.abs((q1_arr[i] - q2_arr[i]) / t_web)

        max_shear_stress[i] = max(tau_cell1, tau_cell2, tau_web)

    return deflection, twist_angle, max_bending_stress, max_shear_stress


def plot_structural_response(y_stations, deflection, twist_angle, max_bending_stress, max_shear_stress):
    """Visualizes the deflections and stresses along the wingspan."""
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))

    axs[0, 0].plot(y_stations, deflection, 'b-', linewidth=2)
    axs[0, 0].set_title('Vertical Bending Deflection', fontweight='bold')
    axs[0, 0].set_ylabel('Deflection [m]')
    axs[0, 0].grid(True, linestyle=':')

    axs[0, 1].plot(y_stations, np.degrees(twist_angle), 'g-', linewidth=2)
    axs[0, 1].set_title('Wing Twist Angle', fontweight='bold')
    axs[0, 1].set_ylabel('Twist [Degrees]')
    axs[0, 1].grid(True, linestyle=':')

    # Handle automatic unit scaling for plot readability if raw arrays are passed
    b_stress_m = max_bending_stress / 1e6 if np.max(max_bending_stress) > 1e4 else max_bending_stress
    s_stress_m = max_shear_stress / 1e6 if np.max(max_shear_stress) > 1e4 else max_shear_stress

    axs[1, 0].plot(y_stations, b_stress_m, 'r-', linewidth=2)
    axs[1, 0].set_title('Maximum Bending Stress (\u03c3_z)', fontweight='bold')
    axs[1, 0].set_ylabel('Stress [MPa]')
    axs[1, 0].set_xlabel('Spanwise Coordinate y [m]')
    axs[1, 0].axhline(y=250, color='red', linestyle='--', label='Yield Limit')
    axs[1, 0].grid(True, linestyle=':')
    axs[1, 0].legend()

    axs[1, 1].plot(y_stations, s_stress_m, 'm-', linewidth=2)
    axs[1, 1].set_title('Maximum Torsional Shear Stress (\u03c4)', fontweight='bold')
    axs[1, 1].set_ylabel('Shear Stress [MPa]')
    axs[1, 1].set_xlabel('Spanwise Coordinate y [m]')
    axs[1, 1].axhline(y=200, color='purple', linestyle='--', label='Yield Limit')
    axs[1, 1].grid(True, linestyle=':')
    axs[1, 1].legend()

    plt.tight_layout()
    plt.show()


def evaluate_wing_design(wing_bays):
    chord = 1.87
    front_spar_pct = 0.15
    main_spar_pct = 0.50
    filepath = "NACA23012.dat"

    # --- DEVICE MASS DEFINITIONS ---
    # Flap Actuators: Placed within flap span (1.1m to 3.8m) near trailing edge
    flap_actuators = [
        {'y': 1.5, 'mass': 8.5, 'x_loc': 0.85 * chord},
        {'y': 2.5, 'mass': 8.5, 'x_loc': 0.85 * chord},
        {'y': 3.5, 'mass': 8.5, 'x_loc': 0.85 * chord}
    ]

    # Aileron Actuators: Placed within aileron span (5.6m to 7.5m) near trailing edge
    aileron_actuators = [
        {'y': 6.0, 'mass': 6.0, 'x_loc': 0.85 * chord},
        {'y': 7.0, 'mass': 6.0, 'x_loc': 0.85 * chord}
    ]

    # Combine into unified function variable expected by nvm2 pipeline
    hld_actuators_combined = flap_actuators + aileron_actuators

    # Spoiler Actuators: Placed outboard of propwash (>3.2m), clear of flaps (<5.6m) at upper mid-chord
    spoiler_actuators = [
        {'y': 4.2, 'mass': 5.0, 'x_loc': 0.55 * chord},
        {'y': 5.0, 'mass': 5.0, 'x_loc': 0.55 * chord}
    ]

    # Empty Fuel Tanks: Distributed over structural inboard wing regions
    fuel_tanks = [
        {'y_start': 0.725, 'y_end': 4.0, 'mass_empty': 35.0}
    ]

    # 1. Run Aerodynamic/Inertial Load Pipeline
    y_steps, shear, moment, torsion, fullwingmass = calculate_and_plot_vmt(
        half_wing_lift_target=10000, b=16.8, y_root=0.725, c=chord, q=571.0,
        flap_start=1.1, flap_end=3.8, delta_cl_flap=1.25,
        aileron_start=5.6, aileron_end=7.5, delta_cl_aileron=0.9,
        filepath=filepath, front_spar_pct=front_spar_pct, main_spar_pct=main_spar_pct, bays=wing_bays,
        engine_mass=122.0, engine_y_loc=2.2, engine_x_loc=-0.1 * chord,
        wtd_mass=20.0, wtd_x_loc=0.25 * chord, n_load_factor=3.8,
        plot=False,

        # Injected Device Elements passed dynamically into backward compatibility layer
        hld_actuators=hld_actuators_combined,
        spoiler_actuators=spoiler_actuators,
        fuel_tank_bays=fuel_tanks
    )

    # 2. Section Properties Extraction
    props_arrays = calculate_section_properties_array(
        y_steps, wing_bays, filepath, chord, front_spar_pct, main_spar_pct
    )

    # 3. Parse A_frontspar and A_rearspar directly from the calculated arrays
    A_front_list = props_arrays["A_frontspar"]
    A_rear_list = props_arrays["A_rearspar"]

    # 4. Deflection and Stress solver
    deflection, twist, max_bend_stress, max_shr_stress = calculate_deflections_and_stresses(
        y_steps, moment, torsion, props_arrays
    )

    # 5. Physics engine output conversion (Pascals -> MegaPascals)
    max_bend_stress_mpa = max_bend_stress / 1e6
    max_shr_stress_mpa = max_shr_stress / 1e6

    total_mass = float(fullwingmass)
    peak_bending = float(np.max(np.abs(max_bend_stress_mpa)))
    peak_shear = float(np.max(np.abs(max_shr_stress_mpa)))

    # 6. Return mapped variables cleanly matching your optimization objective expectations
    return {
        "mass": total_mass,
        "y_steps": y_steps,
        "deflection": deflection,
        "twist": twist,
        "max_bend_stress": max_bend_stress_mpa,
        "max_shr_stress": max_shr_stress_mpa,
        "peak_bending_mpa": peak_bending,
        "peak_shear_mpa": peak_shear,
        "A_frontspar": A_front_list,
        "A_rearspar": A_rear_list
    }

if __name__ == "__main__":
    chord = 1.87
    front_spar_pct = 0.15
    main_spar_pct = 0.50
    filepath = "NACA23012.dat"

    wing_bays = [
        {'y_end': 1.5, 'num_box_stringers': 14, 'num_le_booms': 2,
         't_skin': 0.0005, 't_web': 0.0012, 'A_spar': 8.5E-4, 'A_stringer': 7E-5},

        {'y_end': 3.0, 'num_box_stringers': 13, 'num_le_booms': 2,
         't_skin': 0.0005, 't_web': 0.0010, 'A_spar': 6.0E-4, 'A_stringer': 7E-5},

        {'y_end': 4.5, 'num_box_stringers': 9, 'num_le_booms': 1,
         't_skin': 0.0005, 't_web': 0.0009, 'A_spar': 3E-4, 'A_stringer': 7E-5},

        {'y_end': 5.8, 'num_box_stringers': 6, 'num_le_booms': 0,
         't_skin': 0.0005, 't_web': 0.0006, 'A_spar': 1E-4, 'A_stringer': 7E-5},

        {'y_end': 7.0, 'num_box_stringers': 2, 'num_le_booms': 0,
         't_skin': 0.0005, 't_web': 0.0005, 'A_spar': 1.0E-4, 'A_stringer': 5E-5},

        {'y_end': 8.4, 'num_box_stringers': 1, 'num_le_booms': 1,
         't_skin': 0.00025, 't_web': 0.0005, 'A_spar': 4E-5, 'A_stringer': 1E-5}
    ]


    # Flap Actuators
    flap_actuators = [
        {'y': 1.5, 'mass': 6, 'x_loc': 0.85 * chord},
        {'y': 2.5, 'mass': 6, 'x_loc': 0.85 * chord},
        {'y': 3.5, 'mass': 6, 'x_loc': 0.85 * chord}
    ]

    # Aileron Actuators
    aileron_actuators = [
        {'y': 6.0, 'mass': 6.0, 'x_loc': 0.85 * chord},
        {'y': 7.0, 'mass': 6.0, 'x_loc': 0.85 * chord}
    ]

    hld_actuators_combined = flap_actuators + aileron_actuators

    # Spoiler Actuators
    spoiler_actuators = [
        {'y': 4.2, 'mass': 5.0, 'x_loc': 0.55 * chord},
        {'y': 5.0, 'mass': 5.0, 'x_loc': 0.55 * chord}
    ]

    # Empty Fuel Tanks
    fuel_tanks = [
        {'y_start': 0.725, 'y_end': 4.0, 'mass_empty': 35.0}
    ]

    # For ground load, max fuel
    # fuel_tanks = [{'y_start': 0.725, 'y_end': 4.0, 'mass_empty': 235.0}]

    # 1. Run Aerodynamic/Inertial Load Pipeline with explicit secondary mass injection
    y_steps, shear, moment, torsion, _, fullwingmass, *_ = calculate_and_plot_vmt(
        half_wing_lift_target=10000, b=16.8, y_root=0.725, c=chord, q=571.0,
        flap_start=1.1, flap_end=3.8, delta_cl_flap=1.25,
        aileron_start=5.6, aileron_end=7.5, delta_cl_aileron=0.9,
        filepath=filepath, front_spar_pct=front_spar_pct, main_spar_pct=main_spar_pct, bays=wing_bays,
        engine_mass=122.0, engine_y_loc=2.2, engine_x_loc=-0.1 * chord,
        wtd_mass=20.0, wtd_x_loc=0.25 * chord, n_load_factor=3.8,
        plot=False,
        hld_actuators=hld_actuators_combined,
        spoiler_actuators=spoiler_actuators,
        fuel_tank_bays=fuel_tanks
    )

    props_arrays = calculate_section_properties_array(
        y_steps, wing_bays, filepath, chord, front_spar_pct, main_spar_pct
    )

    deflection, twist, max_bend_stress, max_shr_stress = calculate_deflections_and_stresses(
        y_steps, moment, torsion, props_arrays
    )

    print("\n--- RESULTS SUMMARY ---")
    print(f"Max Tip Deflection:      {deflection[-1]:.3f} m")
    print(f"Max Tip Twist:           {np.degrees(twist[-1]):.2f}°")
    print(f"Peak Bending Stress:     {np.max(max_bend_stress) / 1e6:.2f} MPa")
    print(f"Peak Shear Stress:       {np.max(max_shr_stress) / 1e6:.2f} MPa")
    print(f"TOTAL MASS (Both Halves): {fullwingmass:.2f} kg")

    plot_structural_response(y_steps, deflection, twist, max_bend_stress, max_shr_stress)
   # print(max_bend_stress)
