import numpy as np
import matplotlib.pyplot as plt

# Aircraft and loading parameters
seat_pitch = (29 * 2.54) / 100   # meters
x_first_row = 6.239              # first row position (m)
c_MAC = 2.303                    # Mean Aerodynamic Chord (m)

max_takeoff_weight = 23000
wing_mass = 6258.3
fuselage_mass = 8278.2
wing_moment = wing_mass * 11.4618
fuselage_moment = fuselage_mass * 12.9997
avg_passenger_mass = 84

mass_per_cargo = (7400 - avg_passenger_mass * 56) / 2
x_cargo_front = 4.304
x_cargo_back = 20.977
x_wing = 12.24

# Initial state and LEMAC position
initial_mass = wing_mass + fuselage_mass
initial_moment = wing_moment + fuselage_moment
x_cg_OEW = initial_moment / initial_mass
x_lemac = x_cg_OEW - 0.25 * c_MAC

# Create seat positions (14 rows)
seat_positions = np.linspace(x_first_row, x_first_row + 13 * seat_pitch, 14)

def calculate_cg(new_mass, new_moment, current_mass, current_cg):
    total_mass = current_mass + new_mass
    total_moment = current_mass * current_cg + new_moment
    return total_moment / total_mass, total_mass

def calculate_seat_loading(seat_order, init_cg, init_mass, pax_per_step):
    current_cg, current_mass = init_cg, init_mass
    mac_list, mass_list = [], []
    for pos in seat_order:
        step_mass = pax_per_step * avg_passenger_mass
        step_moment = pos * step_mass
        current_cg, current_mass = calculate_cg(step_mass, step_moment, current_mass, current_cg)
        mac_list.append(((current_cg - x_lemac) / c_MAC) * 100)
        mass_list.append(current_mass)
    return mac_list, mass_list, current_cg, current_mass

def calculate_full_loading(seat_order, cargo_cg, cargo_mass):
    window_mac, window_mass, final_cg, final_mass = calculate_seat_loading(seat_order, cargo_cg, cargo_mass, 2)
    aisle_mac, aisle_mass, final_cg, final_mass = calculate_seat_loading(seat_order, final_cg, final_mass, 2)
    return window_mac, window_mass, aisle_mac, aisle_mass, final_cg, final_mass

# Cargo configurations
fcg, fmass = calculate_cg(mass_per_cargo, mass_per_cargo * x_cargo_front, initial_mass, x_cg_OEW)
bcg, bmass = calculate_cg(mass_per_cargo, mass_per_cargo * x_cargo_back, initial_mass, x_cg_OEW)
both_mass_val = mass_per_cargo * 2
both_moment = mass_per_cargo * x_cargo_front + mass_per_cargo * x_cargo_back
both_cg, both_mass_val = calculate_cg(both_mass_val, both_moment, initial_mass, x_cg_OEW)

# Convert CG positions to %MAC
oew_mac = ((x_cg_OEW - x_lemac) / c_MAC) * 100
fcg_mac = ((fcg - x_lemac) / c_MAC) * 100
bcg_mac = ((bcg - x_lemac) / c_MAC) * 100
both_cg_mac = ((both_cg - x_lemac) / c_MAC) * 100

# Calculate seat loading progressions for two loading scenarios
window_mac_front, window_mass_front, aisle_mac_front, aisle_mass_front, front_final_cg, front_final_mass = \
    calculate_full_loading(seat_positions, both_cg, both_mass_val)
window_mac_back, window_mass_back, aisle_mac_back, aisle_mass_back, back_final_cg, back_final_mass = \
    calculate_full_loading(seat_positions[::-1], both_cg, both_mass_val)

# Connect cargo to seat loading by inserting cargo points into the progressions
window_mac_front.insert(0, both_cg_mac)
window_mac_back.insert(0, both_cg_mac)
window_mass_front.insert(0, both_mass_val)
window_mass_back.insert(0, both_mass_val)

# THESE LINES should align each aisle progression to the correct final window CG:
aisle_mac_front.insert(0, window_mac_front[-1])
aisle_mac_back.insert(0, window_mac_back[-1])
aisle_mass_front.insert(0, window_mass_front[-1])
aisle_mass_back.insert(0, window_mass_back[-1])

# Fuel loading parameters and progression
max_zero_fuel_weight = aisle_mass_front[-1]
max_fuel = max_takeoff_weight - max_zero_fuel_weight
fuel_list = np.linspace(0, max_fuel, 10)

def add_fuel_loading(final_cg, final_mass):
    cg_list, mass_list = [], []
    for fuel in fuel_list:
        fuel_moment = fuel * x_wing
        new_cg, new_mass = calculate_cg(fuel, fuel_moment, final_mass, final_cg)
        cg_list.append(((new_cg - x_lemac) / c_MAC) * 100)
        mass_list.append(new_mass)
    return cg_list, mass_list

fuel_mac, fuel_mass = add_fuel_loading(front_final_cg, front_final_mass)

# Determine overall CG envelope from the loading progressions and fuel loading
all_cg = window_mac_front + window_mac_back + aisle_mac_front + aisle_mac_back + fuel_mac + [both_cg_mac,bcg_mac,fcg_mac,oew_mac]
min_cg = min(all_cg)
max_cg = max(all_cg)
# Define safety margin lines: add 2% at the forward (minimum) side and subtract 2% at the aft (maximum) side
safety_forward = min_cg - 2
safety_aft = max_cg + 2


# Plotting
plt.figure(figsize=(12, 7))

# Cargo connection lines
cargo_colors = ['green', 'magenta', 'orange', 'blue']
plt.plot([oew_mac, fcg_mac], [initial_mass, fmass],
         color=cargo_colors[0], linestyle='--', lw=1.5, label='Cargo Loading Paths')
plt.plot([oew_mac, bcg_mac], [initial_mass, bmass],
         color=cargo_colors[1], linestyle='--', lw=1.5)
plt.plot([fcg_mac, both_cg_mac], [fmass, both_mass_val],
         color=cargo_colors[2], linestyle='--', lw=1.5)
plt.plot([bcg_mac, both_cg_mac], [bmass, both_mass_val],
         color=cargo_colors[3], linestyle='--', lw=1.5)

# Seat and fuel loading progressions with distinct colors and markers
plt.plot(window_mac_front, window_mass_front, 'b-o', label='Front to Back Loading Window')
plt.plot(window_mac_back, window_mass_back, 'r--s', label='Back to Front Loading Window')
plt.plot(aisle_mac_front, aisle_mass_front, 'g-o', label='Front to Back Loading Aisle')
plt.plot(aisle_mac_back, aisle_mass_back, 'c--s', label='Back to Front Loading Aisle')
plt.plot(fuel_mac, fuel_mass, 'k--s', label='Fuel Loading')

# Cargo points
plt.scatter(fcg_mac, fmass, c='orange', marker='X', s=150, zorder=5, label='Front Cargo')
plt.scatter(bcg_mac, bmass, c='purple', marker='X', s=150, zorder=5, label='Back Cargo')
plt.scatter(both_cg_mac, both_mass_val, c='green', marker='X', s=150, zorder=5, label='Both Cargo')

# Reference markers and lines
plt.axvline(25, color='g', linestyle='-.', label='25% MAC (OEW)')
# Safety margin lines (2% margin from the envelope extremes)
plt.axvline(safety_forward, color='black', linestyle='-', lw=1, label='Safety Margin 2%')
plt.axvline(safety_aft, color='black', linestyle='-', lw=1)  # No label here
plt.axvline(max_cg, color='black', linestyle='-', lw=1)
plt.axvline(min_cg, color='black', linestyle='-', lw=1)


plt.scatter(oew_mac, initial_mass, c='black', marker='*', s=200, label='OEW Position', zorder=5)
plt.scatter(aisle_mac_back[-1], max_zero_fuel_weight, c='blue', marker='*', s=200, 
            label='Max Zero Fuel Weight', zorder=5)
plt.scatter(fuel_mac[-1], max_takeoff_weight, c='red', marker='*', s=200, 
            label='Max Take-off Weight', zorder=5)



plt.title("Loading Diagram ATR72-600")
plt.xlabel("CG Position (% MAC)")
plt.ylabel("Total Mass (kg)")
plt.grid(True)
plt.legend(ncol=1, loc='upper left')
plt.xlim(-20, 70)
plt.tight_layout()
plt.show()