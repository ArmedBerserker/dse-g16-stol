from classes.aircraft_2 import Aircraft, loader


def fuselage_cross_section(ac: Aircraft):
    fuselage_width= ac.fuselage.seat_width*ac.fuselage.num_seats + ac.fuselage.aisle_width + ac.fuselage.clearance + 2*ac.fuselage.wall_thickness
    fuselage_height = ac.fuselage.ailse_height + 2*ac.fuselage.wall_thickness + ac.fuselage.under_compartment_storage_height + ac.fuselage.floor_thickness
    fuselage_diameter = (fuselage_width + fuselage_height)/2
    return fuselage_width,fuselage_height, fuselage_diameter

def fuselage_length_components(ac: Aircraft, fuselage_diameter):
    fuselage_tot_length = 0.86 * ac.weights.m_takeoff**0.42
    fuselage_cone_length = ac.fuselage.tail_cone_fuselage_ratio*fuselage_tot_length
    tail_cone_fineness = fuselage_cone_length/fuselage_diameter
    V_50 = 1.3*ac.requirements.general['stall_speed']
    over_nose_angle = ac.fuselage.approach_angle + 0.07 * V_50







if __name__ == '__main__':
    file_path = 'yamls/aircraft.yaml'
    target_class = Aircraft
    fuselage = loader.load(file_path, target_class)


