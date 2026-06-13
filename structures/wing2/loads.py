import numpy as np
import matplotlib.pyplot as plt

def schrenk_half_wing_loading(
        half_wing_lift_target,  # Target lift for ONE wing in N
        b,  # Total aircraft wingspan in m
        y_root,  # Spanwise position where the wing physically starts (fuselage radius) in m
        c,  # Wing chord (constant for rectangular wing) in m
        q,  # Dynamic pressure in Pa
        flap_start,  # Distance from centerline to flap start in m
        flap_end,  # Distance from centerline to flap end in m
        delta_cl_flap,  # Delta Cl for plain flap = 1.25
        aileron_start,  # Distance from centerline to aileron start in m
        aileron_end,  # Distance from centerline to aileron end in m
        delta_cl_aileron,  # Delta Cl for aileron = 0.9
        num_points=500  # Resolution along the actual wing panel
):
    """
    Computes the spanwise running load (distributed force) for a single
    half-wing panel starting at y_root and ending at the wingtip (b/2).

    Returns:
        y_stations: Array of positions along the isolated wing panel (y_root to b/2)
        running_load: Distributed force array (N/m) at each station
    """
    # 1. Setup stations tracking ONLY along the actual physical wing panel
    y_stations = np.linspace(y_root, b / 2, num_points)
    dy = y_stations[1] - y_stations[0]

    # 2. Define the 2D Lift Curve Slope (a0 approx 2*pi per radian)
    a0 = 5.6 # TODO - CHANGE THIS

    # 3. Calculate Base Shape Factor
    # NOTE: The ellipse still maps to full span 'b' to keep true aerodynamic downwash curve shape
    base_shape = 0.5 * q * c * a0 * (1.0 + (4.0 / np.pi) * np.sqrt(1.0 - (2.0 * y_stations / b) ** 2))

    # 4. Build the Control Surface Step Array (The flat additive blocks)
    control_steps = np.zeros(num_points)
    for i, y in enumerate(y_stations):
        # Inboard Flap Step
        if flap_start <= y <= flap_end:
            control_steps[i] += 0.5 * q * c * delta_cl_flap

        # Outboard Aileron Step
        if aileron_start <= y <= aileron_end:
            control_steps[i] += 0.5 * q * c * delta_cl_aileron

    # 5. Integrate control surface lift for ONE half-wing panel
    single_panel_control_lift = np.sum(control_steps) * dy

    # 6. Find the remaining lift this single panel must generate via Alpha
    lift_required_from_alpha = half_wing_lift_target - single_panel_control_lift

    # 7. Integrate the base shape factor for just this single panel
    single_panel_base_lift = np.sum(base_shape) * dy

    # 8. Solve for the operational alpha (radians) matching this isolated panel target
    alpha_wing = lift_required_from_alpha / single_panel_base_lift

    # 9. Calculate the final combined spanwise running load profile
    running_load = (alpha_wing * base_shape) + control_steps

    return y_stations, running_load

if __name__  == "__main__":
    y, load = schrenk_half_wing_loading(
        half_wing_lift_target = 10000, # Target lift for this side (N)
        b=16.8,                        # Total span (m)
        y_root=0.725,                    # Wing start position (m)
        c=1.87,                         # Chord (m)
        q=571.0,                      # Dynamic pressure (Pa)
        flap_start=1.1,                # Flap start (m)
        flap_end=3.8,                  # Flap end (m)
        delta_cl_flap=1.25,
        aileron_start=5.6,             # Aileron start (m)
        aileron_end=7.5,               # Aileron end (m)
        delta_cl_aileron=0.9
    )

    total_lift = np.trapezoid(load, y)
    print(f"Integrated half-wing lift = {total_lift:.2f} N")
    plt.figure()
    plt.plot(y, load, 'b-')
    plt.xlabel("Spanwise Position (m)")
    plt.ylabel("Force (N/m)")
    plt.grid(True)
    plt.show()