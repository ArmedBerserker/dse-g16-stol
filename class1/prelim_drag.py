from classes.aircraft_2 import Aircraft
from lookups.consts import *
import pandas as pd
import numpy as np


def cd0(ac : Aircraft,
        type_to_use : str = "Single Engine Propeller Driven",
        friction_source : str = 'lookups/skin_fric.csv',
        s_wet_source : str = 'lookups/s_wets.csv') -> float:
    
    '''Applies the Roskam method for cd0 estimation'''

    # hence we have an array with first column being skin frics, the middle being a and the right being b
    # log(f) = a + b * log(S_w) in lbs and ft2

    ######### WE FIND CD0 HERE ####################################
    cf_tab = pd.read_csv(friction_source)
    cf_tab = cf_tab.sort_values('c_f').to_numpy() 

    a = np.interp(ac.wing.c_f, cf_tab[:,0], cf_tab[:, 1])
    b = np.interp(ac.wing.c_f, cf_tab[:,0], cf_tab[:, 2])
    
    swet_dict = pd.read_csv(s_wet_source)
    swet_dict = swet_dict.set_index(swet_dict.columns[0]).to_dict('index')

    c = swet_dict[type_to_use]['c']
    d = swet_dict[type_to_use]['d']
    
    # log(Sw) = c + d * log(m_to) in lbs and ft2
    
    # log(f) = a + b * c + b * d * log(m_to))

    logf = a + b * c + b * d * np.log10(ac.weights.m_takeoff * LBS_TO_KG)

    f = 10 ** logf

    cd0 = f / (ac.wing.area * M2_TO_F2)
    
    return cd0

def k(ac : Aircraft) -> float:
    '''Apply the Vos method to estimate Cdi'''
    wing = ac.wing
    assert wing.aspect_ratio != None and wing.psi != None and wing.phi != None, 'Aspect Ratio, Psi, Phi must be defined!'
    e = 1 / (np.pi * ac.wing.aspect_ratio * wing.psi + (1 / wing.phi))
    k = 1 / (np.pi * ac.wing.aspect_ratio * e)
    return k

    
def prelim_drag(ac : Aircraft,
                type_to_use : str = "Single Engine Propeller Driven",
                friction_source : str = 'lookups/skin_fric.csv',
                s_wet_source : str = 'lookups/s_wets.csv') -> float:
    
    '''Returns the max L/D ratio using the Vos+Roskam method'''

    return 0.5 * np.sqrt(1 / k(ac) / cd0(ac, type_to_use, friction_source, s_wet_source))
