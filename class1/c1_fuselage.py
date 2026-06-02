# Fix path FIRST, before any local imports
import sys
import os
import math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from classes.aircraft_2 import Aircraft, loader
from lookups.consts import *

def fuselage_cross_section(ac: Aircraft):
    fus = ac.fuselage
    fus.width = 2 * (fus.wall_thickness + fus.clearance + fus.seat_width) + fus.aisle_width
    fus.height = 2 * fus.wall_thickness + fus.floor_thickness + fus.aisle_height
    fus.effective_diameter = (fus.width + fus.height) / 2


def fuselage_external_dimensions(ac: Aircraft):
    fus = ac.fuselage
    fus.nose_cone_length = fus.nose_cone_fineness_ratio * fus.effective_diameter
    fus.tail_cone_length = fus.tail_cone_fineness_ratio * fus.effective_diameter
    fus.over_nose_angle = fus.approach_angle + 0.07 * 1.3 * ac.requirements.general['stall_speed']

    fus.length = 9
    fus.base_area = 0.0875
    fus.start_cabin = 1.92
    fus.l_cabin = 4.5-1.92
    fus.n_pax = 6
    fus.n_window_seats = 2
    fus.n_middle_seats = 0
    fus.n_aisle_seats = 0
    fus.max_perimeter = math.pi * 0.55 + math.pi * 0.4 + 0.65 + 0.35 + 2 * 0.75
    fus.max_cross_section_area = math.pi * 0.55 ** 2 / 2 + math.pi * 0.4 ** 2 / 2 + 1.35 * 0.75 + 0.65 * 0.4 + 0.35 * 0.55
    internal_area = math.pi * 0.5 ** 2 / 2 + math.pi * 0.35 ** 2 / 2 + 1.35 * 0.75 + 0.65 * 0.4 + 0.35 * 0.55
    fus.vol_cabin_and_cargo = fus.max_cross_section_area * fus.l_cabin
    fus.x_pos_seats = [1.92, 3.5, 4.5]  
    fus.x_cargo_holds = (5.5*1 + 6.25*0.25)/1.25 # 5


def calculate_fuselage_parameters(ac:Aircraft):
    fuselage_cross_section(ac)
    fuselage_external_dimensions(ac)


if __name__ == '__main__':
    file_path = '../yamls/aircraft.yaml'
    target_class = Aircraft
    aircraft = loader.load(file_path, target_class)

    calculate_fuselage_parameters(aircraft)

    # print(aircraft.fuselage.length)
    # print(aircraft.fuselage.height)
    # print(aircraft.fuselage.width)
    # print(aircraft.fuselage.base_area)
    # print(aircraft.fuselage.max_cross_section_area)
    # print(aircraft.fuselage.max_perimeter)
    # print(aircraft.fuselage.start_cabin)
    # print(aircraft.fuselage.l_cabin)
    # print(aircraft.fuselage.vol_cabin_and_cargo)
    # print(aircraft.fuselage.x_pos_seats)
    # print(aircraft.fuselage.x_cargo_holds)
    # print(aircraft.fuselage.n_pax)
    # print(aircraft.fuselage.n_window_seats)
    # print(aircraft.fuselage.n_middle_seats)
    # print(aircraft.fuselage.n_aisle_seats)





# def fuselage_cross_section(ac: Aircraft, update_ac: False):
#     fus = ac.fuselage
#     fuselage_width = (2 * (fus.wall_thickness*M_TO_INCH + fus.clearance*M_TO_INCH + fus.seat_width) +
#                       fus.aisle_width*M_TO_INCH)
#     fuselage_height = (fus.wall_thickness*M_TO_INCH + fus.floor_thickness * M_TO_INCH + fus.aisle_height +
#                        fus.top_compartment_height * M_TO_INCH)
#     fuselage_width *= INCH_TO_M
#     fuselage_height *= INCH_TO_M
#     fuselage_diameter = (fuselage_width + fuselage_height) / 2
#     if update_ac:
#         f = ac.fuselage
#         f.width = fuselage_width
#         f.height = fuselage_height
#         f.eq_diameter = fuselage_diameter

#     return fuselage_width, fuselage_height, fuselage_diameter


# def fuselage_length_components(ac: Aircraft, fuselage_diameter, update_ac: False):
#     fuselage_tot_length = 0.86 * (ac.weights.m_takeoff / LBS_TO_KG) ** 0.42
#     fuselage_tot_length *= FT_TO_M
#     fuselage_cone_length = ac.fuselage.tail_cone_fuselage_ratio * fuselage_tot_length

#     tail_cone_fineness = fuselage_cone_length/fuselage_diameter

#     V_50 = 1.3 * ac.requirements.general['stall_speed']
#     over_nose_angle = ac.fuselage.approach_angle + 0.07 * V_50
#     if update_ac:
#         f = ac.fuselage
#         f.length = fuselage_tot_length
#         f.tail_cone_fuselage_ratio = tail_cone_fineness
#         f.over_nose_angle = over_nose_angle
#         f.tail_cone_length = fuselage_cone_length

#     return fuselage_tot_length, tail_cone_fineness, over_nose_angle,fuselage_cone_length

# def size_fuselage(ac: Aircraft) -> Aircraft:
#     fuselage_cross_section(ac, update_ac=True)
#     fuselage_length_components(ac, ac.fuselage.eq_diameter, update_ac=True)

# if __name__ == '__main__':
#     file_path = 'yamls/aircraft.yaml'
#     target_class = Aircraft
#     aircraft = loader.load(file_path, target_class)

#     w, h, d = fuselage_cross_section(aircraft)

#     l, fr, na, fc = fuselage_length_components(aircraft, d)

#     print(w, h, d)
#     print(l, fr, na,fc)