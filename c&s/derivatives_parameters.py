import sys
import os
import numpy as np
import math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from classes.aircraft_2 import Aircraft, loader
from classes.isa import Atmosphere
from lookups.consts import *

#if __name__=="__main__":

file_path = 'yamls/aircraft.yaml'
target_class = Aircraft
aircraft = loader.load(file_path, target_class)


Ixx = 10262.381  # kg*m2
Iyy = 5087.016 # kg*m2
Izz = 13691.332  # kg*m2
Ixy = 0.093  # kg*m2
Iyz = 0.231  # kg*m2
Ixz = 429.382  # kg*m2

#calculate_aircraft_parameters(1870, 8500, rho[0], 132, 25.65, Ixx, Iyy, Izz, Ixz)
mtom = 1849.597 #[kg]
h_cruise = 8500 #[ft]
V_cruise = 132 #[knots]
S = 31.4 #m2

atmos_model = Atmosphere(h_cruise, 0)
rho = atmos_model.density

hp0 = h_cruise * FT_TO_M  # Pressure altitude in the stationary flight condition [m]
V0 = V_cruise * KTS_TO_MS # True airspeed in the stationary flight condition [m/sec]
C_L_1 = mtom * g * 2 / (rho * V_cruise ** 2 * S)
alpha0 = 1.5 * np.pi / 180 # Angle of attack in the stationary flight condition [rad]
th0 = 0  # Pitch angle in the stationary flight condition [rad]

# Aircraft mass
m = mtom  # Mass [kg]

# Aerodynamic properties
e = 0.783  # Oswald factor [ ]
CD0 = 0.02591  # Zero-lift drag coefficient [ ]
CLa = 5.27  # Slope of C_L-alpha curve [ ]

# Longitudinal stability
Cma = -1.01  # Longitudinal stability [ ]
Cmde = 3.646  # Elevator effectiveness [ ]

# Aircraft geometry
S = S  # Wing area [m^2]
Sh = 3.92 # Stabiliser area [m^2]
Sh_S = Sh / S  # [ ]
lh = 4.96  # Tail length [m]
c = 1.87 # Mean aerodynamic cord [m]
lh_c = lh / c  # [ ]
b = 16.8  # Wing span [m]
bh = 6.34  # Stabiliser span [m]
A = b ** 2 / S  # Wing aspect ratio [ ]
Ah = bh ** 2 / Sh  # Stabiliser aspect ratio [ ]
# Vh_V = 1  # [ ]
ih = 0 # Stabiliser angle of incidence [rad]

# Constant values concerning aircraft inertia
muc = mtom / (rho * S * c)
mub = mtom / (rho * S * b)
KX2 = Ixx / (mtom * b ** 2)
KY2 = Iyy / (mtom * c ** 2)
KZ2 = Izz / (mtom * b ** 2)
KXZ = Ixz / (mtom * b ** 2)
dc = c / V0
db = b/(2*V0) 

# Aerodynamic constants
Cmac = 0  # Moment coefficient about the aerodynamic centre [ ]
CNwa = CLa  # Wing normal force slope [ ]
CNha = 2 * np.pi * Ah / (Ah + 2)  # Stabiliser normal force slope [ ]
depsda = 4 / (A + 2)  # Downwash gradient [ ]

# Lift and drag coefficient
CL = 2 * mtom / (rho * V0 ** 2 * S)  # Lift coefficient [ ]
CD = CD0 + (CLa * alpha0) ** 2 / (np.pi * A * e)  # Drag coefficient [ ]

# Stability derivatives

CX0 = 0.0037308 #mtom * np.sin(th0) / (0.5 * rho * V0 ** 2 * S)
CXu = -0.0000572
CXa = -0.261455  # Positive, see FD lecture notes
CXadot = 0
CXq = -0.455
CXde = -0.0086

CZ0 = 0.00048#- mtom * np.cos(th0) / (0.5 * rho * V0 ** 2 * S)
CZu = 0.00408
CZa = 0.002822
CZadot = 0 #-2 * 0.116 * np.pi / 180 * lh * Sh_S * depsda
CZq =0.00344 #small angle assumption -CLq = CZq
CZde = -1.16

Cm0 = -0.01007  # wrong but doesnt matter
Cmu = 0.00408
Cmadot = 0 
Cmq = -16.857
# CmTc = -0.0064 #dont know dont need?
#
CYb = 0.0715
CYbdot = 0 #assumed zero
CYp = 0.0094
CYr = -0.0093
CYda = 0.0165  # dont know
CYdr = -0.0046  # dont know
#
Clb = -0.188
Clp = -0.6199 
Clr = 0.1058
Clda = 0.02015  # dont know
Cldr = -0.159  # dont know
#
Cnb = 0.06509
Cnbdot = 0
Cnp = -0.02269
Cnr = -0.0771
Cnda = 0.00118  # dont know
Cndr = -0.0453  # dont know


# if __name__=="__main__":
#     file_path = '../yamls/aircraft.yaml'
#     target_class = Aircraft
#     aircraft = loader.load(file_path, target_class)

#     atmos_model = Atmosphere(8500, 20)
#     rho = atmos_model.density

#     Ixx = 0
#     Iyy = 0
#     Izz = 0
#     Ixz = 0

#     calculate_aircraft_parameters(1870, 8500, rho[0], 132, 25.65, Ixx, Iyy, Izz, Ixz)