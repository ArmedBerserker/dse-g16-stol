import sys
import os
import numpy as np
import math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from classes.aircraft_2 import Aircraft, loader
from lookups.consts import *

# C_L_alpha, C_L_
# Values for later: C_L_alpha_dot, C_M_alpha_dot

def C_D_a(mtow, rho_cruise, V_cruise, S, CL_CD, C_L_a):
    C_L_1 = mtow * 2 / (rho_cruise * V_cruise ** 2 * S)

    C_D = CL_CD[0]
    C_L = CL_CD[1]
    closest_index = np.abs(C_L - C_L_1).argmin()
    gradient = (C_L[closest_index + 1] - C_L[closest_index - 1]) / (C_D[closest_index + 1] - C_D[closest_index - 1])
    C_D_alpha = gradient * C_L_a
    return C_D_alpha


def C_m_a(x_cg_fw, x_cg_aft, x_ac_A, C_L_a):
    x_ref = (x_cg_fw + x_cg_aft) / 2
    dcm_dcl = x_ref - x_ac_A
    C_m_alpha = dcm_dcl * C_L_a

    return C_m_alpha


def angle_of_attack_derivatives(mtow, rho_cruise, V_cruise, S, CL_CD, C_L_a):
    C_D_alpha = C_D_a(mtow, rho_cruise, V_cruise, S, CL_CD, C_L_a)
    C_m_alpha = C_m_a()


def K_i_line():
    x1 = -1.0
    y1 = 1.85
    x2 = 0
    y2 = 1
    m = (y2 - y1) / (x2 - x1)
    b = y1 - m * x1
    return lambda x: m * x + b


def C_Y_beta(dihedral_angle, zw_df2, S_O, S, C_L_a_v, V_cruise):
    C_Y_beta_w = -0.00573 * dihedral_angle

    K_i_func = K_i_line()
    K_i = K_i_func(zw_df2)

    C_Y_beta_f = -2 * K_i * (S_O / S)

    k_v = 1  # b_v/r12 - ReLU

    # A_v = b_v ** 2 / S_v
    # a_cruise = math.sqrt(1.4 * 2.87 * T_cruise)
    # M = V_cruise / a_cruise
    # bet = math.sqrt(1 - M ** 2)
    # C_L_a_v = 2 * math.pi * A_v / (2 + math.sqrt(A_v ** 2 * bet ** 2 / ))
    # C_Y_beta_v = C_L_a_v *


def angle_of_sideslip_derivatives():
    pass


if __name__ == '__main__':
    file_path = '../yamls/aircraft.yaml'
    target_class = Aircraft
    aircraft = loader.load(file_path, target_class)

    # angle_of_attack_derivatives(aircraft.weights.m_takeoff)

    # zw_df2 = -0.88
    X_1 = 7.1875
    l_f = 9
    X_O = l_f * (0.378 + 0.527 * (X_1 / l_f))
    S_O = ...

    r12 = 0.75
    b_v = 2.67

    # print(X_O)
    # S_0_S =
    # K_i_func = K_i_line()
    # K_i = K_i_func(zw_df2)
    # print(K_i)
