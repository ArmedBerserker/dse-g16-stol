# Fix path FIRST, before any local imports
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from classes.aircraft_2 import Aircraft, loader
from lookups.consts import *

def fuselage_cross_section(ac: Aircraft, update_ac: False):
    fus = ac.fuselage
    fuselage_width = (2 * (fus.wall_thickness*M_TO_INCH + fus.clearance*M_TO_INCH + fus.seat_width) +
                      fus.aisle_width*M_TO_INCH)
    fuselage_height = (fus.wall_thickness*M_TO_INCH + fus.floor_thickness * M_TO_INCH + fus.aisle_height +
                       fus.top_compartment_height * M_TO_INCH)
    fuselage_width *= INCH_TO_M
    fuselage_height *= INCH_TO_M
    fuselage_diameter = (fuselage_width + fuselage_height) / 2
    if update_ac:
        f = ac.fuselage
        f.width = fuselage_width
        f.height = fuselage_height
        f.eq_diameter = fuselage_diameter

    return fuselage_width, fuselage_height, fuselage_diameter


def fuselage_length_components(ac: Aircraft, fuselage_diameter, update_ac: False):
    fuselage_tot_length = 0.86 * (ac.weights.m_takeoff / LBS_TO_KG) ** 0.42
    fuselage_tot_length *= FT_TO_M
    fuselage_cone_length = ac.fuselage.tail_cone_fuselage_ratio * fuselage_tot_length

    tail_cone_fineness = fuselage_cone_length/fuselage_diameter

    V_50 = 1.3 * ac.requirements.general['stall_speed']
    over_nose_angle = ac.fuselage.approach_angle + 0.07 * V_50
    if update_ac:
        f = ac.fuselage
        f.length = fuselage_tot_length
        f.tail_cone_fuselage_ratio = tail_cone_fineness
        f.over_nose_angle = over_nose_angle
        f.tail_cone_length = fuselage_cone_length

    return fuselage_tot_length, tail_cone_fineness, over_nose_angle,fuselage_cone_length

def size_fuselage(ac: Aircraft) -> Aircraft:
    fuselage_cross_section(ac, update_ac=True)
    fuselage_length_components(ac, ac.fuselage.eq_diameter, update_ac=True)

if __name__ == '__main__':
    file_path = 'yamls/aircraft.yaml'
    target_class = Aircraft
    aircraft = loader.load(file_path, target_class)

    w, h, d = fuselage_cross_section(aircraft)

    l, fr, na, fc = fuselage_length_components(aircraft, d)

    print(w, h, d)
    print(l, fr, na,fc)