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
    fus.start_cabin = 1.7
    fus.l_cabin = 3.8
    fus.n_pax = 6
    fus.n_window_seats = 2
    fus.n_middle_seats = 0
    fus.n_aisle_seats = 0
    fus.max_perimeter = math.pi * 0.55 + math.pi * 0.4 + 0.65 + 0.35 + 2 * 0.75
    fus.max_cross_section_area = math.pi * 0.55 ** 2 / 2 + math.pi * 0.4 ** 2 / 2 + 1.35 * 0.75 + 0.65 * 0.4 + 0.35 * 0.55
    internal_area = math.pi * 0.5 ** 2 / 2 + math.pi * 0.35 ** 2 / 2 + 1.35 * 0.75 + 0.65 * 0.4 + 0.35 * 0.55
    fus.vol_cabin_and_cargo = 1.3 + internal_area * fus.l_cabin
    fus.x_pos_seats = [1.9, 3, 3.6]
    fus.x_cargo_holds = 5


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



