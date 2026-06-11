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


Ixx = 1609.57013  # kg*m2
Iyy = 8767.774  # kg*m2
Izz = 8643.30418  # kg*m2
Ixy = 0  # kg*m2
Iyz = 0  # kg*m2
Ixz = 186.6324083  # kg*m2

#calculate_aircraft_parameters(1870, 8500, rho[0], 132, 25.65, Ixx, Iyy, Izz, Ixz)
mtom = 1870 #[kg]
h_cruise = 8500 #[ft]
V_cruise = 132 #[knots]
S = 31.4 #m2

atmos_model = Atmosphere(h_cruise, 0)
rho = atmos_model.density

hp0 = h_cruise * FT_TO_M  # Pressure altitude in the stationary flight condition [m]
V0 = V_cruise * KTS_TO_MS # True airspeed in the stationary flight condition [m/sec]
C_L_1 = mtom * g * 2 / (rho * V_cruise ** 2 * S)
alpha0 = 0.965 * np.pi / 180 # Angle of attack in the stationary flight condition [rad]
th0 = 0  # Pitch angle in the stationary flight condition [rad]

# Aircraft mass
m = mtom  # Mass [kg]

# Aerodynamic properties
e = 0.783  # Oswald factor [ ]
CD0 = 0.02591  # Zero-lift drag coefficient [ ]
CLa = 5.76  # Slope of C_L-alpha curve [ ]

# Longitudinal stability
Cma = -1.87  # Longitudinal stability [ ]
Cmde = 3.646  # Elevator effectiveness [ ]

# Aircraft geometry
S = S  # Wing area [m^2]
Sh = 3.92 # Stabiliser area [m^2]
Sh_S = Sh / S  # [ ]
lh = 4.96  # Tail length [m]
c = 1.87 # Mean aerodynamic cord [m]
lh_c = lh / c  # [ ]
b = 16.8  # Wing span [m]
bh = 4.85  # Stabiliser span [m]
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

CX0 = mtom * np.sin(th0) / (0.5 * rho * V0 ** 2 * S)
CXu = 0.017142
CXa = 0.04387  # Positive, see FD lecture notes
CXadot = 0
CXq = -9.81
CXde = -0.0086

CZ0 = - mtom * np.cos(th0) / (0.5 * rho * V0 ** 2 * S)
CZu = 0
CZa = -CLa
CZadot = 0 #-2 * 0.116 * np.pi / 180 * lh * Sh_S * depsda
CZq = -11.222 #small angle assumption -CLq = CZq
CZde = -1.16

Cm0 = +0.0297  # wrong but doesnt matter
Cmu = 0.0073
Cmadot = 0 
Cmq = -25.426
# CmTc = -0.0064 #dont know dont need?
#
CYb = 0.19133 
CYbdot = 0 #assumed zero
CYp = -0.000609 
CYr = -0.1547 
CYda = -0.0400  # dont know
CYdr = +0.2300  # dont know
#
Clb = -0.023067 
Clp = -0.53728  
Clr = 0.053059  
Clda = -0.23088  # dont know
Cldr = +0.03440  # dont know
#
Cnb = 0.065995
Cnbdot = 0
Cnp = -0.017503
Cnr = -0.056
Cnda = -0.0120  # dont know
Cndr = -0.0939  # dont know


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