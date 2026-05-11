import numpy as np


def calculate_section_inertia(tskin, tspar, Astr, xfs, xrs, hfs, hrs):
    # Width of the box
    w = xrs - xfs

    # Ixx Calculation
    I_spars = (1 / 12) * tspar * (hfs ** 3 + hrs ** 3)
    int_factor = (w / 12) * (hfs ** 2 + hfs * hrs + hrs ** 2)
    I_skins = 2 * (tskin * int_factor)
    I_stringers = 2 * (Astr / 12) * (hfs ** 2 + hfs * hrs + hrs ** 2)
    Ixx = I_spars + I_skins + I_stringers

    # J Calculation
    # 1. Enclosed Area (Trapezoid)
    Ae = 0.5 * (hfs + hrs) * w

    # 2. Length of the sloping skins (Pyhtagoras)
    slope_length = np.sqrt(w ** 2 + (0.5 * (hfs - hrs)) ** 2)

    # 3. Perimeter Integral (sum of length/thickness for each wall)
    integral = (hfs / tspar) + (hrs / tspar) + (2 * slope_length / tskin)

    J = (4 * Ae ** 2) / integral #Bredt-Batho

    return Ixx, J


