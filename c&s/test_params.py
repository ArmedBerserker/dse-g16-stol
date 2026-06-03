import sys
import os
import numpy as np
import math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from classes.aircraft_2 import Aircraft, loader
from classes.isa import Atmosphere
from lookups.consts import *
#used table citation values from svv pleun
atmos_model = Atmosphere(8500, 20)
rho = atmos_model.density[0]

# Ixx = 0
# Iyy = 0
# Izz = 0
# Ixz = 0

#calculate_aircraft_parameters(1870, 8500, rho[0], 132, 25.65, Ixx, Iyy, Izz, Ixz)
mtom = ((9197 + 2750) * 0.453592) + 101 + 80 + 80 + 85 + 94 + 73 + 65 + 82 + 100
h_cruise = 8500 #[ft]
V_cruise = 66.75 #[m/s]
S = 30 #m2

hp0 = 2100 # Pressure altitude in the stationary flight condition [m]
V0 = 102 # True airspeed in the stationary flight condition [m/sec]
alpha0 = -1.129
th0 = 0.0608

# Aircraft mass
m = mtom  # Mass [kg]

# Aerodynamic properties 
e = 0.654  # Oswald factor [ ]
CD0 = 0.0267  # Zero-lift drag coefficient [ ]
CLa = 0.0752 * (180/np.pi)  # Slope of C_L-alpha curve [rad ]

# Longitudinal stability
Cma = -0.50505 # Longitudinal stability [ ]
Cmde = -1.2  # Elevator effectiveness [ ]

# Aircraft geometry
S = S  # Wing area [m^2]
Sh = 0.2 *S # Stabiliser area [m^2]
Sh_S = Sh / S  # [ ]
lh = 0.71* 5.968 # Tail length [m]
c = 2.0569  # Mean aerodynamic cord [m]
lh_c = lh / c  # [ ]
b = 15.911  # Wing span [m]
bh = 5.791 # Stabiliser span [m]
A = b ** 2 / S  # Wing aspect ratio [ ]
Ah = bh ** 2 / Sh  # Stabiliser aspect ratio [ ]
Vh_V = 1  # [ ]
ih = -2*np.pi/180 # Stabiliser angle of incidence [rad]

# Constant values concerning aircraft inertia
mub = m/(rho*S*b)
muc = m/(rho*S*c)
KX2 = 0.019 #Kxsquared
KY2 = 1.25*1.114 #Kysquared
KZ2 = 0.042 #Kzsquared
KXZ = 0.002
dc = c / V0
db = b/(2*V0) 

# Aerodynamic constants
Cmac = 0  # Moment coefficient about the aerodynamic centre [ ]
CNwa = CLa  # Wing normal force slope [ ]
CNha = 2 * np.pi * Ah / (Ah + 2)  # Stabiliser normal force slope [ ]
depsda = 4 / (A + 2)  # Downwash gradient [ ]

# Lift and drag coefficient
CL = 2 * mtom *9.81/ (rho * V0 ** 2 * S)  # Lift coefficient [ ]
CD = CD0 + (CL) ** 2 / (np.pi * A * e)  # Drag coefficient [ ]
print(CL)

# Stability derivatives

CX0 = mtom*9.81 * np.sin(th0) / (0.5 * rho * V0 ** 2 * S)
CXu = -2 *CD
CXa = CL * (1-2*CLa/ (np.pi *A*e))  # Positive, see FD lecture notes
CXadot = 0.0833
CXq = -0.28170
CXde = -0.03728

CZ0 = - mtom *9.81* np.cos(th0) / (0.5 * rho * V0 ** 2 * S)
CZu = -2 * CL
CZa = -(CLa + CD)
CZadot = -0.00350
CZq = -5.66290 #small angle assumption -CLq = CZq
CZde = -0.69612

Cm0 = +0.0297  # wrong but doesnt matter
Cmu = 0.0400
Cmadot = +0.17800  #dont know
Cmq = -8.79415
CmTc = -0.0064 #dont know dont need?
#
CYb = -0.75
CYbdot = 0 #assumed zero
CYp = -0.0304
CYr = -0.8495
CYda = -0.0400  # dont know
CYdr = +0.2300  # dont know
#
Clb = -0.10260
Clp = -0.71085  
Clr = 0.23760  
Clda = -0.23088  # dont know
Cldr = +0.03440  # dont know
#
Cnb = 0.1348
Cnbdot = 0
Cnp = -0.0602
Cnr = -0.2061
Cnda = -0.0120  # dont know
Cndr = -0.0939  # dont know

