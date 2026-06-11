import numpy as np
from scipy.optimize import fsolve

# Figure 12.12 - x (crcv) / y (tau)
x_data = np.array([
    0.0055, 0.0124, 0.0193, 0.0262, 0.0332, 0.0401, 0.0470, 0.0539,
    0.0608, 0.0677, 0.0746, 0.0815, 0.0884, 0.0954, 0.1023, 0.1092,
    0.1161, 0.1230, 0.1299, 0.1368, 0.1437, 0.1507, 0.1576, 0.1645,
    0.1714, 0.1783, 0.1852, 0.1921, 0.1990, 0.2059, 0.2129, 0.2198,
    0.2267, 0.2336, 0.2405, 0.2474, 0.2543, 0.2612, 0.2682, 0.2751,
    0.2820, 0.2889, 0.2958, 0.3027, 0.3096, 0.3165, 0.3235, 0.3304,
    0.3373, 0.3442, 0.3511, 0.3580, 0.3649, 0.3718, 0.3787, 0.3857,
    0.3926, 0.3995, 0.4064, 0.4133, 0.4202, 0.4271, 0.4340, 0.4410,
    0.4479, 0.4548, 0.4617, 0.4686, 0.4755, 0.4824, 0.4893, 0.4963,
    0.5032, 0.5101, 0.5170, 0.5239, 0.5308, 0.5377, 0.5446, 0.5515,
    0.5585, 0.5654, 0.5723, 0.5792, 0.5861, 0.5930, 0.5999, 0.6068,
    0.6138, 0.6207, 0.6276, 0.6345, 0.6414, 0.6483, 0.6552, 0.6621,
    0.6691, 0.6760, 0.6829, 0.6898, 0.6967
])
y_data = np.array([
    0.0144, 0.0378, 0.0609, 0.0835, 0.1054, 0.1263, 0.1460, 0.1645,
    0.1818, 0.1980, 0.2133, 0.2278, 0.2414, 0.2544, 0.2667, 0.2784,
    0.2896, 0.3003, 0.3106, 0.3206, 0.3302, 0.3395, 0.3486, 0.3576,
    0.3663, 0.3749, 0.3833, 0.3916, 0.3998, 0.4079, 0.4159, 0.4239,
    0.4317, 0.4395, 0.4471, 0.4546, 0.4620, 0.4692, 0.4763, 0.4832,
    0.4899, 0.4965, 0.5029, 0.5091, 0.5152, 0.5212, 0.5271, 0.5329,
    0.5386, 0.5442, 0.5498, 0.5553, 0.5608, 0.5662, 0.5716, 0.5770,
    0.5823, 0.5876, 0.5929, 0.5982, 0.6035, 0.6087, 0.6140, 0.6192,
    0.6244, 0.6296, 0.6348, 0.6400, 0.6452, 0.6503, 0.6554, 0.6604,
    0.6655, 0.6704, 0.6754, 0.6803, 0.6851, 0.6900, 0.6948, 0.6995,
    0.7043, 0.7090, 0.7137, 0.7184, 0.7231, 0.7278, 0.7325, 0.7373,
    0.7420, 0.7467, 0.7514, 0.7562, 0.7609, 0.7657, 0.7704, 0.7752,
    0.7800, 0.7847, 0.7895, 0.7943, 0.7990
])

def get_tau(crcv):
    return np.interp(crcv, x_data, y_data)
def get_crcv(tau):
    return np.interp(tau, y_data, x_data)

# both guesstimated from flight dynamics reader
Cnb = 0.18
Cyb = -0.8

# assumed
Cdy = 0.5
dc = 0.89  # centroid distance to cg

S = 31.4
MTOW = 1840
rho = 1.079 # 2000ft at ISA +20
CLmax = 1.94
Vs = np.sqrt(MTOW * 9.81/(0.5*rho*CLmax*S))
b = 16.8
yt = 1.42 + (1.45/2)
Vmc = 0.95 * Vs
Vapp = 1.3 * Vs
Vw = 0.2 * Vs
q = 0.5 * rho * Vmc * Vmc
T = 5475 

brbv = 0.8 # 0.7-1.0
delta_r_max = np.deg2rad(30)
clav = 6
eta = 1.0 # dynamic pressure ratio at the tail
volume = 0.056 # vertical tail volume coefficient
lv = 9.9 # actually xv but whatever, adsee is fucked
Sv = volume * S * b / lv
Ss = 16.88 # side area of aircraft profile
svs = Sv/S

# Condition 1 - Asymmetric Thrust
cndeltar = (T*yt)/(q*S*b*delta_r_max)
tau = (cndeltar) / (clav*volume*eta*brbv)
if tau >= 0.8:
    raise Exception("Tau larger than 0.8")
crcv = get_crcv(tau)
if crcv >= 0.5:
    raise Exception("Cr/Cv larger than 0.5")
print(cndeltar, tau,  crcv)  # tau = control surface effectiveness vs chord ratio, cr/cv = like cf/c

# Condition 2 - Crosswind Landing
Fw = 0.5 * rho * Vw * Vw * Ss * Cdy
cydeltar = clav * eta * tau * brbv * svs
beta = np.arctan(Vw/Vapp)
Vtot = np.sqrt( Vapp*2 + Vw*2 )
print(cydeltar, beta, Vtot) # 

def func(x): # x = [sigma, deltaR]
    beta_term = beta - x[0]
    return [
        (0.5 * rho * Vtot ** 2 * S * b) *
        (Cnb * beta_term + cndeltar * x[1])
        + (Fw * dc * np.cos(x[0])),

        Fw - (0.5 * rho * Vtot ** 2 * S) *
        (Cyb * beta_term + cydeltar * x[1])
    ]
root = fsolve(func, [1, 1])

print(np.rad2deg(root))  # crab angle (<10, def <15) and rudder defl. (>1)