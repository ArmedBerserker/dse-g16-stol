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


def fuselage_external_dimensions(ac: Aircraft, x_pos_seats, x_cargo_holds, nose_length, fuselage_length=0):
    fus = ac.fuselage
    fus.nose_cone_length = fus.nose_cone_fineness_ratio * fus.effective_diameter
    fus.tail_cone_length = fus.tail_cone_fineness_ratio * fus.effective_diameter
    fus.over_nose_angle = fus.approach_angle + 0.07 * 1.3 * ac.requirements.general['stall_speed']
    fus.x_pos_seats = x_pos_seats
    fus.x_cargo_holds = x_cargo_holds
    fus.nose_length = nose_length
    fus.start_cabin = fus.x_pos_seats[0] - fus.seat_length / 2
    if fuselage_length == 0:
        fus.length = (0.86 * (ac.weights.m_takeoff / LBS_TO_KG) ** 0.42) * FT_TO_M
    else:
        fus.length = fuselage_length
    fus.l_cabin = fus.length - fus.tail_cone_length - fus.nose_length - fus.cockpit_length * INCH_TO_M  # 3.8
    

    fus.base_area = 0.1
    fus.n_pax = 6
    fus.n_window_seats = 2
    fus.n_middle_seats = 0
    fus.n_aisle_seats = 0
    
    fus.max_perimeter = math.pi * fus.upper_corner_radii + (fus.width - fus.upper_corner_radii * 2) + math.pi * fus.lower_corner_radii + (fus.width - fus.lower_corner_radii * 2) + (fus.height - fus.upper_corner_radii - fus.lower_corner_radii) * 2 # math.pi * fus.upper_corner_radii + math.pi * fus.lower_corner_radii + 0.65 + 0.35 + 2 * 0.75
    fus.max_cross_section_area = math.pi * fus.upper_corner_radii ** 2 / 2 + math.pi * fus.lower_corner_radii ** 2 / 2 + (fus.width - fus.upper_corner_radii * 2) * fus.upper_corner_radii + (fus.width - fus.lower_corner_radii * 2) * fus.lower_corner_radii + (fus.height - fus.upper_corner_radii - fus.lower_corner_radii) * fus.width # math.pi * fus.upper_corner_radii ** 2 / 2 + math.pi * fus.lower_corner_radii ** 2 / 2 + 1.35 * 0.75 + 0.65 * fus.lower_corner_radii + 0.35 * fus.upper_corner_radii
    internal_area = math.pi * (fus.upper_corner_radii - fus.wall_thickness) ** 2 / 2 + math.pi * (fus.lower_corner_radii - fus.wall_thickness) ** 2 / 2 + (fus.width - 2 * fus.wall_thickness) * (fus.height - fus.upper_corner_radii - fus.lower_corner_radii) + (fus.width - 2 * fus.upper_corner_radii) * (fus.lower_corner_radii - fus.wall_thickness) + (fus.width - 2 * fus.lower_corner_radii) * fus.lower_corner_radii + (fus.width - 2 * fus.upper_corner_radii) * fus.upper_corner_radii # 0.65 * fus.lower_corner_radii + 0.35 * fus.upper_corner_radii
    fus.vol_cabin_and_cargo = fus.cargo_length * fus.cargo_width * fus.cargo_height + internal_area * fus.l_cabin    


def calculate_fuselage_parameters(ac:Aircraft, x_pos_seats, x_cargo_holds, nose_length, fuselage_length):
    fuselage_cross_section(ac)
    fuselage_external_dimensions(ac, x_pos_seats, x_cargo_holds, nose_length, fuselage_length)


if __name__ == '__main__':
    file_path = '../yamls/aircraft.yaml'
    target_class = Aircraft
    aircraft = loader.load(file_path, target_class)

    x_pos_seats = [2.3, 3.5, 4.2]
    x_cargo_holds = 5.5
    nose_length = 2
    fuselage_length = 12
    calculate_fuselage_parameters(aircraft, x_pos_seats, x_cargo_holds, nose_length, fuselage_length)

    print(aircraft.fuselage.length)
    print(aircraft.fuselage.height)
    print(aircraft.fuselage.width)
    print(aircraft.fuselage.base_area)
    print(aircraft.fuselage.max_cross_section_area)
    print(aircraft.fuselage.max_perimeter)
    print(aircraft.fuselage.start_cabin)
    print(aircraft.fuselage.l_cabin)
    print(aircraft.fuselage.vol_cabin_and_cargo)
    print(aircraft.fuselage.x_pos_seats)
    print(aircraft.fuselage.x_cargo_holds)
    print(aircraft.fuselage.n_pax)
    print(aircraft.fuselage.n_window_seats)
    print(aircraft.fuselage.n_middle_seats)
    print(aircraft.fuselage.n_aisle_seats)



