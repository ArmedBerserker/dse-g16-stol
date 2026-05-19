from classes.aircraft_2 import Aircraft, loader
from lookups.consts import *

def fuselage_cross_section(ac: Aircraft):
    fus = ac.fuselage
    fuselage_width = (2 * (fus.wall_thickness + fus.clearance + fus.seat_width) +
                      fus.aisle_width)
    fuselage_height = (fus.wall_thickness + fus.floor_thickness * M_TO_INCH + fus.aisle_height +
                       fus.top_compartment_height * M_TO_INCH)
    fuselage_width *= INCH_TO_M
    fuselage_height *= INCH_TO_M
    fuselage_diameter = (fuselage_width + fuselage_height) / 2

    return fuselage_width, fuselage_height, fuselage_diameter


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


