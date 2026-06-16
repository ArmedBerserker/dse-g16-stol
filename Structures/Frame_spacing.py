import numpy as np
import pandas as pd
from datetime import datetime
import math
from scipy.interpolate import interp1d

''' TODO:
    - Add curved plate functions
    - Set up stringer spacing and finsing points they are at along sections (make it symmetric, start at top or bottom center)
    - Check correct material properties and thicknesses are used everywhere
    - Stringer pitch, area, and skin thickness plot with contours areas colors mass or if viable'''

material_database = {  # https://www.aerospacemetals.com/wp-content/uploads/2023/06/Aluminum-2024-T3.pdf
    'AA2024-T3': {
        'E': 73100000000,
        'G': 28000000000,
        'density': 2780,
        'TYS': 345000000,
        'UTS': 483000000,
        'CYS': 345000000,
        'UCS': 483000000,
        'tau max': 207000000,
        'nu': 0.33,
        'alpha': 0.8,
        'n': 0.6
    },
    'AA7075-T6': {  # https://www.matweb.com/search/datasheet.aspx?MatGUID=4f19a42be94546b686bbf43f79c51b7d&ckck=1
        'E': 71700000000,
        'G': 26900000000,
        'density': 2810,
        'TYS': 503000000,
        'UTS': 572000000,
        'CYS': 503000000,
        'UCS': 572000000,
        'tau max': 331000000,
        'nu': 0.33,
        'alpha': 0.8,
        'n': 0.6
    },
    'AA7050-T7451': {  # https://www.matweb.com/search/DataSheet.aspx?MatGUID=142262cf7fbc4c83917ca5c3d17df1ed 
        'E': 71700000000,
        'G': 26900000000,
        'density': 2830,
        'TYS': 469000000,
        'UTS': 524000000,
        'CYS': 469000000,
        'UCS': 524000000,
        'tau max': 303000000,
        'nu': 0.33,
        'alpha': 0.8,
        'n': 0.6
    },
    'AA7050-T7651': {  # https://www.matweb.com/search/DataSheet.aspx?MatGUID=142262cf7fbc4c83917ca5c3d17df1ed 
        'E': 71700000000,
        'G': 26900000000,
        'density': 2830,
        'TYS': 490000000,
        'UTS': 552000000,
        'CYS': 490000000,
        'UCS': 552000000,
        'tau max': 324000000,
        'nu': 0.33,
        'alpha': 0.8,
        'n': 0.6
    }
}

# Add L, U
stringers_database = {
    'j': {
        'boundary conditions': ['SSFS', 'SSSS', 'SSFS', 'SSFS'],
        'rel length': np.array([0.863, 2.46, 0.94, 0.863])
    },
    'z': {
        'boundary conditions': ['SSFS', 'SSSS', 'SSFS'],
        'rel length': np.array([0.9, 0.6, 0.9])},
    'formed z': {
        'boundary conditions': ['SSFS', 'SSSS', 'SSSS', 'SSFS'],
        'rel length': np.array([0.33, 0.9, 2.44, 0.9])
    },
}

load_cases = {
    'take-off': {
        'ax': 0.5*9.81,
        'mass frac': 1
    },
    'landing': {
        'ax': 0.6*9.81,
        'mass frac': 0.9
    }
}

def thin_plate_buckling_stress(C, material_properties: dict, t, b):
    return C * np.pi**2 * material_properties['E'] / (12 * (1 * material_properties['nu']**2)) * (t / b)**2

def we2(C, t, material_properties: dict, cripling_stress_stiffener):
    return t * np.sqrt(C * np.pi**2 / (12 * (1 * material_properties['nu']**2))) * np.sqrt(material_properties['E'] / cripling_stress_stiffener)

def stiffener_crippling_stress(stringer_area, material_properties: dict, stiffener_type: str, l_stringer, t: float = 0.0015):
    stringer = stringers_database[stiffener_type]
    bounds = stringer['boundary conditions']
    b_rel = stringer['rel length']
    B_scale = (stringer_area - t**2 * (len(b_rel)-1)) / (t * np.sum(b_rel))
    b = B_scale * b_rel
    crip_stress_elements = np.zeros_like(b)
    A_elem = np.zeros_like(b)
    for i, b_part in enumerate(b):
        C = C_flat_plate(bounds[i], a_b=l_stringer/b_part)
        # print(f'C={C}')
        crip_stress_frac = material_properties['alpha'] * (C / material_properties['CYS'] * np.pi**2 * material_properties['E'] / (12 * (1 * material_properties['nu']**2)) * (t / b_part)**2)**(1 - material_properties['n'])
        if crip_stress_frac < 1:
            crip_stress_elements[i] = crip_stress_frac * material_properties['CYS']
        else:
            crip_stress_elements[i] = material_properties['CYS']
        A_elem[i] = t * b[i]
    return (crip_stress_elements@A_elem) / np.sum(A_elem)

def C_flat_plate(boundary_condition: str, a_b):
    C_30 = {
        'SSFF': 0.026,
        'CCFF': 0.107
    }
    if boundary_condition == 'SSFF' or boundary_condition == 'CCFF':
        return C_30.get(boundary_condition)
    else:
        df = pd.read_csv("Structures/buckling_coefficients.csv")
        x = df['a/b'].values
        y = df[boundary_condition].values

        f = interp1d(x, y, kind='linear', fill_value='extrapolate')
        return float(f(a_b))

def C(b_t):
    if b_t<40:
        return 4.0
    elif b_t>110:
        return 6.98
    else:
        return 70/2.98 * (b_t - 40) + 4.0

def pannel_buckling_stress(stringer_area, stringer_type: str, material_properties_skin: dict, material_properties_stringer: dict, b_pitch, a_pitch, t_skin: float = 0.032*0.0254, t_stringer: float = 0.0015):
    # skin_buckling_stress = thin_plate_buckling_stress(C=4.0, material_properties=material_properties_skin, t=t_skin, b=b_pitch)
    crip_stress_stiff = stiffener_crippling_stress(stringer_area, material_properties_stringer, stiffener_type=stringer_type,l_stringer=a_pitch, t=t_stringer)
    b_t = b_pitch / t_skin
    we2_current = we2(C(b_t), t_skin, material_properties_skin, crip_stress_stiff)
    b_t_new = (b_pitch - we2_current) / t_skin
    C_stiff_plate = C(b_t_new)
    new_skin_buck_stress = C_stiff_plate * np.pi**2 * material_properties_skin['E'] / (12 * (1 * material_properties_skin['nu']**2)) * (1 / b_t_new)**2
    num1 = (stringer_area + we2_current * t_skin) * crip_stress_stiff
    num2 = (b_pitch - we2_current) * t_skin * new_skin_buck_stress
    return (num1 + num2) / (stringer_area + b_pitch * t_skin)

def mass_est(n_stringers, stringer_area, n_frames, frame_area, t_skin, length_section, perimeter_skin, material_properties_skin: dict, material_properties_stringer: dict):
    m_stringers = n_stringers * stringer_area * length_section * material_properties_stringer['density']
    m_skin = perimeter_skin * t_skin * length_section * material_properties_skin['density']
    m_frames = n_frames * frame_area * perimeter * material_properties_stringer['density']
    return m_frames + m_skin + m_stringers

if __name__ == '__main__':
    # Variables
    MTOM = 1840  # kg
    safety_factor = 1.5
    l_fus = 11
    l_nc = 3.4650000000000007
    l_tc = 4.429687500000001
    l_section = l_fus - l_nc - l_tc

    # Define cross section (8 plates)
    h = 1.7
    w = 1.45
    r_upper = 0.55
    r_lower = 0.4
    s_sections = [w-2*r_upper, np.pi/2*r_upper, h-r_upper-r_lower, 
                  np.pi/2*r_lower, w-2*r_lower, np.pi/2*r_lower, 
                  h-r_upper-r_lower, np.pi/2*r_upper]
    perimeter = sum(s_sections)
    print(perimeter, l_section)

    # Define loads from accelerations
    F = []
    for load_case in load_cases.keys():
        F.append(load_cases[load_case]['mass frac'] * MTOM * safety_factor * load_cases[load_case]['ax'])
    F_crit = max(F)

    # Material and stringer selection
    stringer_type = 'z'
    material_names = list(material_database.keys())
    material_skin = 'AA2024-T3'
    material_stringer = 'AA2024-T3'
    material_properties_skin = material_database[material_skin]
    material_properties_stringer = material_database[material_stringer]

    # Define geometry search range
    # t_skin_vals = np.arange(0.0005, 0.007, 0.0001)  # 0.5-1.5mm, steps of 0.1mm
    # t_stringer_vals = np.arange(0.0008, 0.001, 0.0002)  # 0.8-3.4mm, steps of 0.2mm
    # A_stringer = (1e-6) * np.arange(30, 45, 10)  # 30-250mm^2, steps of 10mm^2
    # A_frame = 60*1e-6  # Set to constant 60mm^2 for now (not included in buckling calculations)
    # stringer_pitch = np.arange(0.1, 0.21, 0.1)  # 100-200mm, steps of 10mm
    # frame_pitch = np.arange(0.3, 0.51, 0.1)  # 300-500mm, steps of 10mm

    # t_skin_vals = np.arange(0.0003, 0.0009, 0.0001)  # 0.5-1.5mm, steps of 0.1mm
    # t_stringer_vals = np.arange(0.0008, 0.001, 0.0001)  # 0.8-3.4mm, steps of 0.2mm
    # A_stringer = (1e-6) * np.arange(35, 75, 5)  # 30-250mm^2, steps of 10mm^2
    # A_frame = 20*1e-6  # Set to constant 60mm^2 for now (not included in buckling calculations)
    # stringer_pitch = np.arange(0.1, 0.22, 0.02)  # 100-200mm, steps of 10mm
    # frame_pitch = np.arange(0.3, 0.5, 0.05)  # 300-500mm, steps of 10mm
    t_skin_vals = np.arange(0.0003, 0.00039, 0.0001)  # 0.5-1.5mm, steps of 0.1mm
    t_stringer_vals = np.arange(0.0012954, 0.002, 0.001)  # 0.8-3.4mm, steps of 0.2mm
    A_stringer = (1e-6) * np.arange(70, 75, 5)  # 30-250mm^2, steps of 10mm^2
    A_frame = 70*1e-6  # Set to constant 60mm^2 for now (not included in buckling calculations)
    stringer_pitch = np.arange(0.2, 0.22, 0.02)  # 100-200mm, steps of 10mm
    frame_pitch = np.arange(0.45, 0.5, 0.05)  # 300-500mm, steps of 10mm

    # Check how to size each section so they have the same stress they can take
    viable_solutions = []

    for h, stringer_type in enumerate(list(stringers_database.keys())):
        for i, t_skin in enumerate(t_skin_vals):
            for j, t_stringer in enumerate(t_stringer_vals):
                for k, A_string in enumerate(A_stringer):
                    for l, b_pitch in enumerate(stringer_pitch):
                        n_stringers = math.ceil(perimeter / b_pitch)

                        # Load applied
                        A_cross_section = n_stringers * A_string + perimeter * t_skin
                        applied_stress = F_crit / A_cross_section
                        for m, a_pitch in enumerate(frame_pitch):
                            n_frames = math.ceil(l_section / a_pitch)

                            # Load it can withstand
                            buck_stress_pan = pannel_buckling_stress(A_string, stringer_type, material_properties_skin, material_properties_stringer, b_pitch, a_pitch, t_skin, t_stringer)
                            viable_solution: bool = applied_stress <= buck_stress_pan
                            if viable_solution:
                                print(f'Buckling stress panel: {buck_stress_pan}, applied: {applied_stress}')
                                mass = mass_est(n_stringers, A_string, n_frames, A_frame, t_skin, l_section, perimeter, material_properties_skin, material_properties_stringer)
                                margin = (buck_stress_pan - applied_stress) / applied_stress * 100
                                # Change e for rivet type: 4 = flathead, 3 = brazier head, 1 = countersunk
                                e = 1
                                p = t_skin / np.sqrt(12 * (1 - material_database[material_skin]['nu']**2) * applied_stress / (e * np.pi**2 * material_database[material_skin]['E'])) # * 0.000145038
                                viable_solutions.append({"mass [kg]": mass,
                                                         "rivet spacing [mm]": p*1000,
                                                         "stringer type": stringer_type,
                                                         "skin thickness [mm]": t_skin*1e3, 
                                                         "stringer thickness [mm]": t_stringer*1e3, 
                                                         "stringer cross-section area [mm2]": A_string*1e6,
                                                         "frame cross-section area [mm2]": A_frame*1e6,
                                                         "stringer pitch [mm]": b_pitch*1e3,
                                                         "frame pitch [mm]": a_pitch*1e3,
                                                         "stringer/frame material": material_stringer,
                                                         "skin material": material_skin,
                                                         "stress margin [%]": margin})
    if len(viable_solutions) == 0:
        print(f'No viable solutions found')
    else:
        df = pd.DataFrame(viable_solutions)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # Lowest mass first
        df = df.sort_values("mass [kg]")
        path_all = f"Structures/all_viable_solutions_{timestamp}.csv"
        df.to_csv(path_all, index=False)
        print(f'Viable solutions saved to {path_all}')

        best_10 = df.nsmallest(10, "mass [kg]")
        path_best = f"Structures/lowest_mass_solutions_{timestamp}.csv"
        best_10.to_csv(path_best, index=False)

    # Figure out how to space frames
