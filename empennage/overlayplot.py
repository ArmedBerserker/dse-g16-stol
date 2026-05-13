import numpy as np
import matplotlib.pyplot as plt

def compute_loading_diagram(
    wing_mass, fuselage_mass, wing_moment, fuselage_moment, seat_count, seat_pitch,
    x_first_row, c_MAC, max_takeoff_weight, avg_passenger_mass, x_cargo_front,
    x_cargo_back, x_wing, cargo_per_side, pax_per_step=2
):
    """
    Compute all the CG and mass data for a given aircraft configuration.
    Returns a dict of relevant arrays and scalars for plotting.
    """
    # Initial mass/moment
    initial_mass = wing_mass + fuselage_mass
    initial_moment = wing_moment + fuselage_moment
    x_cg_OEW = initial_moment / initial_mass  # Operating Empty Weight CG

    # LEMAC
    x_lemac = x_cg_OEW - 0.25 * c_MAC

    def calc_cg(new_mass, new_moment, current_mass, current_cg):
        total_mass = current_mass + new_mass
        total_moment = current_mass * current_cg + new_moment
        return total_moment / total_mass, total_mass

    # Create seat positions
    seat_positions = np.linspace(x_first_row, x_first_row + (seat_count - 1)*seat_pitch, seat_count)

    def calc_seat_loading(seat_order, init_cg, init_mass, pps):
        cur_cg, cur_mass = init_cg, init_mass
        mac_list, mass_list = [], []
        for pos in seat_order:
            step_mass = pps * avg_passenger_mass
            step_moment = pos * step_mass
            cur_cg, cur_mass = calc_cg(step_mass, step_moment, cur_mass, cur_cg)
            mac_list.append(((cur_cg - x_lemac)/c_MAC)*100)
            mass_list.append(cur_mass)
        return mac_list, mass_list, cur_cg, cur_mass

    def calc_full_loading(seat_order, cargoCG, cargoMass):
        w_mac, w_mass, final_cg, final_mass = calc_seat_loading(seat_order, cargoCG, cargoMass, pax_per_step)
        a_mac, a_mass, final_cg, final_mass = calc_seat_loading(seat_order, final_cg, final_mass, pax_per_step)
        return w_mac, w_mass, a_mac, a_mass, final_cg, final_mass

    # Cargo
    fcg, fmass = calc_cg(cargo_per_side, cargo_per_side*x_cargo_front, initial_mass, x_cg_OEW)
    bcg, bmass = calc_cg(cargo_per_side, cargo_per_side*x_cargo_back, initial_mass, x_cg_OEW)
    both_total = cargo_per_side * 2
    both_moment = cargo_per_side*x_cargo_front + cargo_per_side*x_cargo_back
    both_cg, both_mass = calc_cg(both_total, both_moment, initial_mass, x_cg_OEW)

    # Convert CGs to %MAC
    oew_mac  = ((x_cg_OEW - x_lemac)/c_MAC)*100
    fcg_mac  = ((fcg - x_lemac)/c_MAC)*100
    bcg_mac  = ((bcg - x_lemac)/c_MAC)*100
    both_mac = ((both_cg - x_lemac)/c_MAC)*100

    # Seat loading
    wmf, wmassf, amf, amassf, fin_cg_f, fin_mass_f = calc_full_loading(seat_positions, both_cg, both_mass)
    wmb, wmassb, amb, amassb, fin_cg_b, fin_mass_b = calc_full_loading(seat_positions[::-1], both_cg, both_mass)

    # Insert cargo CG at start
    wmf.insert(0, both_mac)
    wmb.insert(0, both_mac)
    wmassf.insert(0, both_mass)
    wmassb.insert(0, both_mass)

    # Align aisle with final window CG
    amf.insert(0, wmf[-1])
    amassf.insert(0, wmassf[-1])
    amb.insert(0, wmb[-1])
    amassb.insert(0, wmassb[-1])

    # Fuel loading
    max_zero_fuel_weight = amassf[-1]
    max_fuel = max_takeoff_weight - max_zero_fuel_weight
    fuel_list = np.linspace(0, max_fuel, 10)

    fuel_mac_list, fuel_mass_list = [], []
    for fuel in fuel_list:
        fuel_moment = fuel*x_wing
        new_cg, new_mass = calc_cg(fuel, fuel_moment, fin_mass_f, fin_cg_f)
        fuel_mac_list.append(((new_cg - x_lemac)/c_MAC)*100)
        fuel_mass_list.append(new_mass)

    # CG envelope extremes
    all_cg = wmf + wmb + amf + amb + fuel_mac_list + [both_mac, bcg_mac, fcg_mac, oew_mac]
    min_cg, max_cg = min(all_cg), max(all_cg)
    safety_forward = min_cg - 2
    safety_aft    = max_cg + 2

    return {
        "oew_mac": oew_mac,
        "initial_mass": initial_mass,
        "fcg_mac": fcg_mac,
        "fmass": fmass,
        "bcg_mac": bcg_mac,
        "bmass": bmass,
        "both_cg_mac": both_mac,
        "both_mass_val": both_mass,
        "wm_front": wmf,
        "wmass_front": wmassf,
        "am_front": amf,
        "amass_front": amassf,
        "wm_back": wmb,
        "wmass_back": wmassb,
        "am_back": amb,
        "amass_back": amassb,
        "fuel_mac": fuel_mac_list,
        "fuel_mass": fuel_mass_list,
        "max_zero_fuel_weight": max_zero_fuel_weight,
        "max_takeoff_weight": max_takeoff_weight,
        "min_cg": min_cg,
        "max_cg": max_cg,
        "safety_forward": safety_forward,
        "safety_aft": safety_aft,
    }

# =============== SCENARIO 1 (all black, alpha=0.5) ===============
scenario1 = compute_loading_diagram(
    wing_mass       = 6258.3,
    fuselage_mass   = 8278.2,
    wing_moment     = 6258.3 * 11.4618,
    fuselage_moment = 8278.2 * 12.9997,
    seat_count      = 14,
    seat_pitch      = (29 * 2.54)/100,
    x_first_row     = 6.239,
    c_MAC           = 2.303,
    max_takeoff_weight   = 23000,
    avg_passenger_mass   = 84,
    x_cargo_front   = 4.304,
    x_cargo_back    = 20.977,
    x_wing          = 12.24,
    cargo_per_side  = (7400 - 84*56)/2
)

# =============== SCENARIO 2 (colored, fully opaque) ===============
scenario2 = compute_loading_diagram(
    wing_mass       = 6601,
    fuselage_mass   = 6693,
    wing_moment     = 6601 * 11.50221254,
    fuselage_moment = 6693 * 12.30065292,
    seat_count      = 18,
    seat_pitch      = (29 * 2.54)/100,
    x_first_row     = 6.239,
    c_MAC           = 2.303,
    max_takeoff_weight   = 23000,
    avg_passenger_mass   = 84,
    x_cargo_front   = 4.304,
    x_cargo_back    = 20.977,
    x_wing          = 12.24,
    cargo_per_side  = (7400 - 84*72)/2
)

plt.figure(figsize=(12, 7))

# ----- Scenario 1: black lines, alpha=0.5 -----
plt.plot(
    [scenario2["oew_mac"], scenario2["fcg_mac"]],
    [scenario2["initial_mass"], scenario2["fmass"]],
    color='black', linestyle='--', lw=1, alpha=0.5,
    label='Cargo Loading Paths (Scen1)'
)
plt.plot(
    [scenario2["oew_mac"], scenario2["bcg_mac"]],
    [scenario2["initial_mass"], scenario2["bmass"]],
    color='black', linestyle='--', lw=1, alpha=0.5
)
plt.plot(
    [scenario2["fcg_mac"], scenario2["both_cg_mac"]],
    [scenario2["fmass"], scenario2["both_mass_val"]],
    color='black', linestyle='--', lw=1, alpha=0.5
)
plt.plot(
    [scenario2["bcg_mac"], scenario2["both_cg_mac"]],
    [scenario2["bmass"], scenario2["both_mass_val"]],
    color='black', linestyle='--', lw=1, alpha=0.5
)

plt.plot(
    scenario2["wm_front"], scenario2["wmass_front"],
    'k-o', alpha=0.5, label='Front->Back Window (Scen1)'
)
plt.plot(
    scenario2["wm_back"], scenario2["wmass_back"],
    'k--s', alpha=0.5, label='Back->Front Window (Scen1)'
)
plt.plot(
    scenario2["am_front"], scenario2["amass_front"],
    'k-o', alpha=0.5, label='Front->Back Aisle (Scen1)'
)
plt.plot(
    scenario2["am_back"], scenario2["amass_back"],
    'k--s', alpha=0.5, label='Back->Front Aisle (Scen1)'
)
plt.plot(
    scenario2["fuel_mac"], scenario2["fuel_mass"],
    'k--^', alpha=0.5, label='Fuel Loading (Scen1)'
)

# Points (you can also apply alpha here if desired)
plt.scatter(
    scenario2["oew_mac"], scenario2["initial_mass"],
    c='black', marker='*', s=100, zorder=5,
    label='OEW (Scen1)'
)
plt.scatter(
    scenario2["fcg_mac"], scenario2["fmass"],
    c='black', marker='X', s=80, zorder=5
)
plt.scatter(
    scenario2["bcg_mac"], scenario2["bmass"],
    c='black', marker='X', s=80, zorder=5
)
plt.scatter(
    scenario2["both_cg_mac"], scenario2["both_mass_val"],
    c='black', marker='X', s=80, zorder=5
)
plt.scatter(
    scenario2["am_back"][-1], scenario2["max_zero_fuel_weight"],
    c='black', marker='o', s=80, zorder=5,
    label='MZFW (Scen1)'
)
plt.scatter(
    scenario2["fuel_mac"][-1], scenario2["max_takeoff_weight"],
    c='black', marker='s', s=80, zorder=5,
    label='MTOW (Scen1)'
)

# CG lines (Scenario 1)
plt.axvline(scenario2["min_cg"], color='black', linestyle='-', lw=1, alpha=0.5)
plt.axvline(scenario2["max_cg"], color='black', linestyle='-', lw=1, alpha=0.5)
plt.axvline(scenario2["safety_forward"], color='black', linestyle='-', lw=1, alpha=0.5)
plt.axvline(scenario2["safety_aft"], color='black', linestyle='-', lw=1, alpha=0.5)


# ----- Scenario 2: colored lines (full opacity) -----
cargo_colors = ['green','magenta','orange','blue']
plt.plot(
    [scenario1["oew_mac"], scenario1["fcg_mac"]],
    [scenario1["initial_mass"], scenario1["fmass"]],
    color=cargo_colors[0], linestyle='--', lw=1.5,
    label='Cargo Paths (Scen2)'
)
plt.plot(
    [scenario1["oew_mac"], scenario1["bcg_mac"]],
    [scenario1["initial_mass"], scenario1["bmass"]],
    color=cargo_colors[1], linestyle='--', lw=1.5
)
plt.plot(
    [scenario1["fcg_mac"], scenario1["both_cg_mac"]],
    [scenario1["fmass"], scenario1["both_mass_val"]],
    color=cargo_colors[2], linestyle='--', lw=1.5
)
plt.plot(
    [scenario1["bcg_mac"], scenario1["both_cg_mac"]],
    [scenario1["bmass"], scenario1["both_mass_val"]],
    color=cargo_colors[3], linestyle='--', lw=1.5
)

plt.plot(
    scenario1["wm_front"], scenario1["wmass_front"],
    'b-o', label='Front->Back Window (Scen2)'
)
plt.plot(
    scenario1["wm_back"], scenario1["wmass_back"],
    'r--s', label='Back->Front Window (Scen2)'
)
plt.plot(
    scenario1["am_front"], scenario1["amass_front"],
    'g-o', label='Front->Back Aisle (Scen2)'
)
plt.plot(
    scenario1["am_back"], scenario1["amass_back"],
    'c--s', label='Back->Front Aisle (Scen2)'
)
plt.plot(
    scenario1["fuel_mac"], scenario1["fuel_mass"],
    'k--s', label='Fuel Loading (Scen2)'
)

# Points
plt.scatter(
    scenario1["oew_mac"], scenario1["initial_mass"],
    c='black', marker='*', s=150, zorder=5,
    label='OEW (Scen2)'
)
plt.scatter(
    scenario1["fcg_mac"], scenario1["fmass"],
    c='orange', marker='X', s=150, zorder=5,
    label='Front Cargo (Scen2)'
)
plt.scatter(
    scenario1["bcg_mac"], scenario1["bmass"],
    c='purple', marker='X', s=150, zorder=5,
    label='Back Cargo (Scen2)'
)
plt.scatter(
    scenario1["both_cg_mac"], scenario1["both_mass_val"],
    c='green', marker='X', s=150, zorder=5,
    label='Both Cargo (Scen2)'
)
plt.scatter(
    scenario1["am_back"][-1], scenario1["max_zero_fuel_weight"],
    c='blue', marker='*', s=200, zorder=5,
    label='MZFW (Scen2)'
)
plt.scatter(
    scenario1["fuel_mac"][-1], scenario1["max_takeoff_weight"],
    c='red', marker='*', s=200, zorder=5,
    label='MTOW (Scen2)'
)

# CG lines (Scenario 2) in a different color if desired (here: gray)
plt.axvline(scenario1["min_cg"], color='gray', linestyle='-', lw=1)
plt.axvline(scenario1["max_cg"], color='gray', linestyle='-', lw=1)
plt.axvline(scenario1["safety_forward"], color='gray', linestyle='-', lw=1)
plt.axvline(scenario1["safety_aft"], color='gray', linestyle='-', lw=1)

# A reference line at 25% MAC
plt.axvline(25, color='g', linestyle='-.', label='25% MAC')

plt.title("Overlay: Two ATR72-600 Loading Diagrams")
plt.xlabel("CG Position (% MAC)")
plt.ylabel("Total Mass (kg)")
plt.grid(True)
plt.legend(ncol=2, loc='upper left')
plt.tight_layout()
plt.show()