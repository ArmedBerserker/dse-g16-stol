import sys
import os
import numpy as np
import math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from classes.aircraft_2 import Aircraft, loader
from classes.isa import Atmosphere
from lookups.consts import *

def calculate_aircraft_parameters(mtom, h_cruise, rho, V_cruise, S, Ixx, Iyy, Izz, Ixz):

    hp0 = h_cruise * FT_TO_M  # Pressure altitude in the stationary flight condition [m]
    V0 = V_cruise * KTS_TO_MS # True airspeed in the stationary flight condition [m/sec]
    C_L_1 = mtom * g * 2 / (rho * V_cruise ** 2 * S)
    alpha0 = 0.965 * np.pi / 180 # Angle of attack in the stationary flight condition [rad]
    th0 = 0  # Pitch angle in the stationary flight condition [rad]

    # Aircraft mass
    m = mtom  # Mass [kg]

    # Aerodynamic properties
    e = 0.8  # Oswald factor [ ]
    CD0 = 0.4  # Zero-lift drag coefficient [ ]
    CLa = 5.76  # Slope of C_L-alpha curve [ ]

    # Longitudinal stability
    Cma = -1.87  # Longitudinal stability [ ]
    Cmde = 3.646  # Elevator effectiveness [ ]

    # Aircraft geometry
    S = S  # Wing area [m^2]
    Sh = 7.26 # Stabiliser area [m^2]
    Sh_S = Sh / S  # [ ]
    lh = 6.4  # Tail length [m]
    c = 1.792  # Mean aerodynamic cord [m]
    lh_c = lh / c  # [ ]
    b = 15.2  # Wing span [m]
    bh = 5.65  # Stabiliser span [m]
    A = b ** 2 / S  # Wing aspect ratio [ ]
    Ah = bh ** 2 / Sh  # Stabiliser aspect ratio [ ]
    Vh_V = 1  # [ ]
    ih = 0 # Stabiliser angle of incidence [rad]

    # Constant values concerning aircraft inertia
    muc = mtom / (rho * S * c)
    mub = mtom / (rho * S * b)
    KX2 = Ixx / (mtom * b ** 2)
    KY2 = Iyy / (mtom * c ** 2)
    KZ2 = Izz / (mtom * b ** 2)
    KXZ = Ixz / (mtom * b ** 2)

    # Aerodynamic constants
    Cmac = 0  # Moment coefficient about the aerodynamic centre [ ]
    CNwa = CLa  # Wing normal force slope [ ]
    CNha = 2 * np.pi * Ah / (Ah + 2)  # Stabiliser normal force slope [ ]
    depsda = 4 / (A + 2)  # Downwash gradient [ ]

    # Lift and drag coefficient
    CL = 2 * mtom / (rho * V0 ** 2 * S)  # Lift coefficient [ ]
    CD = CD0 + (CLa * alpha0) ** 2 / (np.pi * A * e)  # Drag coefficient [ ]

    # Stability derivatives

    CX0 = mtom * np.sin(th0) / (0.5 * rho * V0 ** 2 * S)
    CXu = 0.017142
    CXa = 0.04387  # Positive, see FD lecture notes
    CXadot = 0
    CXq = -9.81
    CXde = -0.0086

    CZ0 = - mtom * np.cos(th0) / (0.5 * rho * V0 ** 2 * S)
    CZu = 0
    CZa = -CLa
    CZadot = -2 * 0.116 * np.pi / 180 * lh * Sh_S * depsda
    CZq = 11.222
    CZde = -1.16
    
    # Cm0 = +0.0297  # wrong but doesnt matter
    # Cmu = +0.06990
    # Cmadot = +0.17800
    # Cmq = -8.79415
    # CmTc = -0.0064
    #
    # CYb = -0.7500  # right
    # CYbdot = 0
    # CYp = -0.0304  # right
    # CYr = +0.8495  # right
    # CYda = -0.0400  # wrong
    # CYdr = +0.2300  # right
    #
    # Clb = -0.10260  # right
    # Clp = -0.71085  # right
    # Clr = +0.23760  # right
    # Clda = -0.23088  # wrong
    # Cldr = +0.03440  # right
    #
    # Cnb = +0.1348  # right
    # Cnbdot = 0
    # Cnp = -0.0602  # right
    # Cnr = -0.2061  # right
    # Cnda = -0.0120  # wrong
    # Cndr = -0.0939  # right


if __name__=="__main__":
    file_path = '../yamls/aircraft.yaml'
    target_class = Aircraft
    aircraft = loader.load(file_path, target_class)

    atmos_model = Atmosphere(8500, 20)
    rho = atmos_model.density

    Ixx = 0
    Iyy = 0
    Izz = 0
    Ixz = 0

    calculate_aircraft_parameters(1870, 8500, rho[0], 132, 25.65, Ixx, Iyy, Izz, Ixz)