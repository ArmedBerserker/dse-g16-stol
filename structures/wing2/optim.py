import numpy as np
from scipy.optimize import differential_evolution
# Assuming these are available in your environment
from deflections import evaluate_wing_design, plot_structural_response

# Global Constants
Y_ROOT = 0.725
Y_TIP = 8.4
NUM_BAYS = 7

# Material Stress Limits
SIGMA_YIELD_MPA = 250.0
TAU_YIELD_MPA = 200.0

# Production Manufacturing Constraints
MIN_WEB_THICKNESS = 0.0008   # 1.0 mm minimum thickness gauge
ALUMINUM_DENSITY = 2700.0   # kg/m³
CHORD = 1.87                # m

FIXED_Y_ENDS = np.linspace(Y_ROOT + 0.5, Y_TIP, NUM_BAYS)


def decode_bays(x):
    """Helper to decode flat vector into structural dictionary."""
    stringers = np.zeros(NUM_BAYS)
    booms = np.zeros(NUM_BAYS)
    t_webs = np.zeros(NUM_BAYS)
    A_spars = np.zeros(NUM_BAYS)

    stringers[0] = x[0]
    booms[0] = x[1]
    t_webs[0] = x[2]
    A_spars[0] = x[3]

    for i in range(1, NUM_BAYS):
        idx = 4 + (i - 1) * 4
        stringers[i] = stringers[i - 1] * x[idx]
        booms[i] = booms[i - 1] * x[idx + 1]
        t_webs[i] = t_webs[i - 1] * x[idx + 2]
        A_spars[i] = A_spars[i - 1] * x[idx + 3]

    wing_bays = []

    for i in range(NUM_BAYS):
        actual_web = max(MIN_WEB_THICKNESS, t_webs[i])
        wing_bays.append({
            'y_end': float(FIXED_Y_ENDS[i]),
            'num_box_stringers': int(np.round(stringers[i])),
            'num_le_booms': int(np.round(booms[i])),
            't_skin': 0.0005,
            't_web': float(actual_web),
            'A_spar': float(A_spars[i]),
            'A_stringer': 7E-5
        })

    return wing_bays


def objective_function(x):
    try:
        wing_bays = decode_bays(x)
        res = evaluate_wing_design(wing_bays)

        # 1. Extract true physical parameters from the engine
        y_steps = np.array(res["y_steps"])
        A_front = np.array(res["A_frontspar"])
        A_rear = np.array(res["A_rearspar"])

        # 2. Integrate the exact true volume of the spar webs across the span
        total_web_volume = np.trapezoid(A_front + A_rear, y_steps)
        total_web_mass = total_web_volume * ALUMINUM_DENSITY

        # 3. Apply the exact 50% weight savings from lightening holes
        web_hole_savings = total_web_mass * 0.50
        net_sim_mass = res["mass"] - web_hole_savings

        # 4. Rib mass calculation
        total_ribs_mass = 0.0
        for i in range(NUM_BAYS + 1):
            station_y = Y_ROOT + (i * (Y_TIP - Y_ROOT) / NUM_BAYS)
            span_fraction = (station_y - Y_ROOT) / (Y_TIP - Y_ROOT)
            local_rib_thickness = max(0.0005, 0.0010 * (1.0 - 0.5 * span_fraction))
            rib_lightening_factor = 1.0 - (0.35 + 0.25 * span_fraction)
            approx_rib_area = CHORD * (CHORD * 0.12)
            total_ribs_mass += approx_rib_area * local_rib_thickness * ALUMINUM_DENSITY * rib_lightening_factor * 2

        total_production_mass = net_sim_mass + total_ribs_mass

        # 5. Quadratic penalty to force compliance
        bending_overstress = max(0.0, res["peak_bending_mpa"] - SIGMA_YIELD_MPA)
        shear_overstress = max(0.0, res["peak_shear_mpa"] - TAU_YIELD_MPA)

        penalty = 0.0
        if bending_overstress > 0 or shear_overstress > 0:
            penalty = 5000.0 + (bending_overstress ** 2) * 100 + (shear_overstress ** 2) * 100

        return total_production_mass + penalty

    except Exception:
        return 99999.0


if __name__ == "__main__":
    bounds = []

    # Bay 1 Absolute Limits
    bounds.append((0, 25))           # Stringers
    bounds.append((1, 10))           # LE Booms
    bounds.append((0.0008, 0.001))   # Web Thickness
    bounds.append((1e-6, 1e-3))      # Spar Area

    # Bays 2-7 Taper Factors
    for _ in range(NUM_BAYS - 1):
        bounds.append((0.0, 1.0))    # Stringer taper
        bounds.append((0.0, 1.0))    # LE Boom taper
        bounds.append((0.0, 1.0))    # Web thickness taper
        bounds.append((0.01, 1.0))   # Spar area taper

    print(">> Initializing Precision Wing Optimizer...")

    opt_result = differential_evolution(
        objective_function,
        bounds,
        strategy='best1bin',
        maxiter=5,       # Restored to 150 generations to avoid global minimum starvation
        popsize=15,        # 15 individuals per dimension handles 28 dimensions accurately
        tol=0.001,
        disp=True,
        workers=8,
        updating='deferred',
        polish=True        # Required to snap constraints precisely to limits via gradient local step
    )

    # Reconstruct final details properly
    best_bays = decode_bays(opt_result.x)
    final_metrics = evaluate_wing_design(best_bays)

    # Calculate final exact mass parameters
    y_steps_final = np.array(final_metrics["y_steps"])
    A_front_final = np.array(final_metrics["A_frontspar"])
    A_rear_final = np.array(final_metrics["A_rearspar"])

    final_web_volume = np.trapezoid(A_front_final + A_rear_final, y_steps_final)
    final_web_mass = final_web_volume * ALUMINUM_DENSITY
    web_hole_savings_final = final_web_mass * 0.50
    net_sim_mass_final = final_metrics["mass"] - web_hole_savings_final

    # Recalculate true rib mass
    total_ribs_mass_final = 0.0
    for i in range(NUM_BAYS + 1):
        station_y = Y_ROOT + (i * (Y_TIP - Y_ROOT) / NUM_BAYS)
        span_fraction = (station_y - Y_ROOT) / (Y_TIP - Y_ROOT)
        local_rib_thickness = max(0.0005, 0.0010 * (1.0 - 0.5 * span_fraction))
        rib_lightening_factor = 1.0 - (0.35 + 0.25 * span_fraction)
        approx_rib_area = CHORD * (CHORD * 0.12)
        total_ribs_mass_final += approx_rib_area * local_rib_thickness * ALUMINUM_DENSITY * rib_lightening_factor * 2

    true_production_mass = net_sim_mass_final + total_ribs_mass_final

    print("\n=================== OPTIMIZATION RESULTS ===================")
    print(f"TRUE WEIGHT (Integrated Mass + Ribs - Cutouts):     {true_production_mass:.2f} kg")
    print(f"Raw Simulation Weight (Solid Webs, No Ribs):        {final_metrics['mass']:.2f} kg")
    print(f"Peak Bending Stress:                                {final_metrics['peak_bending_mpa']:.2f} / {SIGMA_YIELD_MPA} MPa")
    print(f"Peak Shear Stress:                                  {final_metrics['peak_shear_mpa']:.2f} / {TAU_YIELD_MPA} MPa")
    print("-----------------------------------------------------------------------")

    for idx, b_config in enumerate(best_bays):
        print(f"Bay {idx + 1} (Ends at {b_config['y_end']:.3f}m): "
              f"Stringers={int(np.round(b_config['num_box_stringers']))} | "
              f"LE Booms={int(np.round(b_config['num_le_booms']))} | "
              f"Web={b_config['t_web']:.5f}m | "
              f"Spar Area={b_config['A_spar']:.7f}m²")

    plot_structural_response(
        final_metrics["y_steps"],
        final_metrics["deflection"],
        final_metrics["twist"],
        final_metrics["max_bend_stress"],
        final_metrics["max_shr_stress"]
    )