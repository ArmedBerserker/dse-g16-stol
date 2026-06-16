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


Ixx = 10078.193  # kg*m2
Iyy = 5034.751 # kg*m2
Izz = 13474.965  # kg*m2
Ixy = -0.105  # kg*m2
Iyz = 0.219  # kg*m2
Ixz = 430.507  # kg*m2

#calculate_aircraft_parameters(1870, 8500, rho[0], 132, 25.65, Ixx, Iyy, Izz, Ixz)
mtom = 1837.806 #[kg]
h_cruise = 2590.8 #[m]
V_cruise = 132 #[knots]
S = 31.4 #m2

atmos_model = Atmosphere(h_cruise, 0)
rho = float(atmos_model.density)
hp0 = h_cruise  # Pressure altitude in the stationary flight condition [m]
V0 = V_cruise * KTS_TO_MS # True airspeed in the stationary flight condition [m/sec]
CL = mtom * g * 2 / (rho * V0 ** 2 * S)
print(CL)
alpha0 = 1.5 * np.pi / 180 # Angle of attack in the stationary flight condition [rad]
th0 = 0  # Pitch angle in the stationary flight condition [rad]

# Aircraft mass
m = mtom  # Mass [kg]

# Aerodynamic properties
e = 0.783  # Oswald factor [ ]
CD0 = 0.02591  # Zero-lift drag coefficient [ ]
CLa = 5.32  # Slope of C_L-alpha curve [ ]

# Longitudinal stability
Cma = -1.057  # Longitudinal stability [ ]
Cmde = -1.7215 # Elevator effectiveness [ ]

# Aircraft geometry
S = S  # Wing area [m^2]
Sh = 6.7 # Stabiliser area [m^2]
Sh_S = Sh / S  # [ ]
#lh = 4.96  # Tail length [m]
c = 1.87 # Mean aerodynamic cord [m]
#lh_c = lh / c  # [ ]
b = 16.8  # Wing span [m]
#bh = 6.34  # Stabiliser span [m]
A = b ** 2 / S  # Wing aspect ratio [ ]
#Ah = bh ** 2 / Sh  # Stabiliser aspect ratio [ ]
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
#Cmac = 0  # Moment coefficient about the aerodynamic centre [ ]
CNwa = CLa  # Wing normal force slope [ ]
#CNha = 2 * np.pi * Ah / (Ah + 2)  # Stabiliser normal force slope [ ]
depsda = 4 / (A + 2)  # Downwash gradient [ ]

# Lift and drag coefficient
# CL = 2 * mtom / (rho * V0 ** 2 * S)  # Lift coefficient [ ]
CD = CD0 + (CLa * alpha0) ** 2 / (np.pi * A * e)  # Drag coefficient [ ]
print("CD",CD)
# Stability derivatives

CX0 = 0.0037133 #mtom * np.sin(th0) / (0.5 * rho * V0 ** 2 * S)
CXu = -0.0534
CXa = -0.26165  # Positive, see FD lecture notes
CXadot = 0
CXq = -0.4569102
CXde = -0.0114175

CZ0 = -0.2271869#- mtom * np.cos(th0) / (0.5 * rho * V0 ** 2 * S)
CZu = -CL
CZa = -5.325
CZadot = 0 #-2 * 0.116 * np.pi / 180 * lh * Sh_S * depsda
CZq =-9.3594 #small angle assumption -CLq = CZq
CZde = -0.0027

Cm0 = -0.00776
Cmu = -0.0007061
Cmadot = 0 
Cmq = -17.034
# CmTc = -0.0064 #dont know dont need?
#
CYb = -0.375
CYbdot = 0 #assumed zero
CYp = -0.0299074
CYr = 0.233864
CYda = 0.0061034 
CYdr = 0.13313
#
Clb = -0.171
Clp = -0.6153418
Clr = 0.0977866
Clda = -0.1595  
Cldr = 0.01874
#
Cnb = 0.04723
Cnbdot = 0
Cnp = -0.02735
Cnr = -0.0681
Cnda = 0.0018233
Cndr = -0.0424711


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
print("rho =", rho)
print("muc =", muc)
print("CL =", CL)
print("CX0 =", CX0)
print("CZ0 =", CZ0)