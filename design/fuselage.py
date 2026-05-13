from classes.aircraft_2 import Aircraft, loader
from lookups.consts import *

def fuselage_cross_section(ac: Aircraft):
    fuselage_width= ac.fuselage.seat_width * ac.fuselage.num_seats + ac.fuselage.aisle_width + ac.fuselage.clearance + 2 * ac.fuselage.wall_thickness
    fuselage_height = ac.fuselage.ailse_height + 2 * ac.fuselage.wall_thickness + ac.fuselage.under_compartment_storage_height * M_TO_INCH + ac.fuselage.floor_thickness * M_TO_INCH
    fuselage_width *= INCH_TO_M
    fuselage_height *= INCH_TO_M
    fuselage_diameter = (fuselage_width + fuselage_height) / 2

    return fuselage_width,fuselage_height, fuselage_diameter


def fuselage_length_components(ac: Aircraft, fuselage_diameter):
    fuselage_tot_length = 0.86 * (ac.weights.m_takeoff / LBS_TO_KG) ** 0.42
    fuselage_tot_length *= FT_TO_M
    fuselage_cone_length = ac.fuselage.tail_cone_fuselage_ratio * fuselage_tot_length
    tail_cone_fineness = fuselage_cone_length/fuselage_diameter
    V_50 = 1.3 * ac.requirements.general['stall_speed']
    over_nose_angle = ac.fuselage.approach_angle + 0.07 * V_50
    return fuselage_tot_length, tail_cone_fineness, over_nose_angle


if __name__ == '__main__':
    file_path = '../yamls/aircraft.yaml'
    target_class = Aircraft
    aircraft = loader.load(file_path, target_class)

    w, h, d = fuselage_cross_section(aircraft)

    l, fr, na = fuselage_length_components(aircraft, d)

    print(w, h, d)
    print(l, fr, na)


