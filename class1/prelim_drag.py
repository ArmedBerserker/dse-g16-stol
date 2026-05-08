"""
Preliminary aircraft drag estimation utilities.

Uses Roskam correlations to estimate zero-lift drag, CD0, and the Vos method to
estimate the induced drag factor, k. These are combined to estimate the maximum
lift-to-drag ratio, (L/D)_max, for early aircraft sizing.
"""

from classes.aircraft_2 import Aircraft
from lookups.consts import *
import pandas as pd
import numpy as np


def cd0(ac : Aircraft,
        type_to_use : str = "Single Engine Propeller Driven",
        friction_source : str = 'lookups/skin_fric.csv',
        s_wet_source : str = 'lookups/s_wets.csv') -> float:
    
    """
    Estimate the zero-lift drag coefficient, CD0, using the Roskam preliminary
    drag build-up method.

    This function estimates the equivalent parasite drag area, f, from empirical
    Roskam relationships of the form:

        log10(f) = a + b * log10(S_wet)
        log10(S_wet) = c + d * log10(W_TO)

    where f and S_wet are in imperial units. The coefficients a and b are
    interpolated from a skin-friction lookup table using the aircraft wing
    skin-friction coefficient, while c and d are selected from a wetted-area
    lookup table based on the aircraft type.

    Parameters
    ----------
    ac : Aircraft
        Aircraft object containing the required weight and wing data.
        Required fields are:
        - ac.wing.c_f : wing skin-friction coefficient
        - ac.wing.area : wing reference area in m^2
        - ac.weights.m_takeoff : takeoff mass in kg

    type_to_use : str, optional
        Aircraft type used to select the wetted-area correlation coefficients
        from `s_wet_source`. Default is "Single Engine Propeller Driven".

    friction_source : str, optional
        Path to the CSV file containing the skin-friction correlation data.
        The file is expected to contain columns for c_f, a, and b.
        Default is "lookups/skin_fric.csv".

    s_wet_source : str, optional
        Path to the CSV file containing wetted-area correlation coefficients.
        The first column is used as the aircraft type key, and the table must
        contain columns `c` and `d`. Default is "lookups/s_wets.csv".

    Returns
    -------
    float
        Estimated zero-lift drag coefficient, CD0, nondimensional.

    Notes
    -----
    The Roskam correlations are applied using imperial units. Therefore, the
    takeoff mass and wing area are converted before evaluating CD0.
    """

    # hence we have an array with first column being skin frics, the middle being a and the right being b
    # log(f) = a + b * log(S_w) in lbs and ft2

    ######### WE FIND CD0 HERE ####################################
    cf_tab = pd.read_csv(str(friction_source))
    cf_tab = cf_tab.sort_values('c_f').to_numpy() 

    a = np.interp(ac.wing.c_f, cf_tab[:,0], cf_tab[:, 1])
    b = np.interp(ac.wing.c_f, cf_tab[:,0], cf_tab[:, 2])
    
    swet_dict = pd.read_csv(s_wet_source)
    swet_dict = swet_dict.set_index(swet_dict.columns[0]).to_dict('index')

    c = swet_dict[type_to_use]['c']
    d = swet_dict[type_to_use]['d']
    
    # log(Sw) = c + d * log(m_to) in lbs and ft2
    
    # log(f) = a + b * c + b * d * log(m_to))

    logf = a + b * c + b * d * np.log10(ac.weights.m_takeoff / LBS_TO_KG)

    f = 10 ** logf

    cd0 = f / (ac.wing.area * M2_TO_F2)
    
    return cd0

def k(ac : Aircraft) -> tuple[float, float]:
    """
    Estimate the induced drag factor, k, using the Vos method.

    This function calculates the Oswald efficiency factor, e, from the Vos
    empirical relationship:

        e = 1 / (pi * AR * psi + 1 / phi)

    The induced drag factor is then calculated as:

        k = 1 / (pi * AR * e)

    where AR is the wing aspect ratio, and phi and psi are empirical correction
    factors for the aircraft configuration.

    Parameters
    ----------
    ac : Aircraft
        Aircraft object containing the required wing aerodynamic parameters.
        Required fields are:
        - ac.wing.aspect_ratio : wing aspect ratio
        - ac.wing.psi : Parasite drag parameter
        - ac.wing.phi : Spanwise efficiency coefficient 

    Returns
    -------
    float
        Induced drag factor, k, used in the drag polar:

            CD = CD0 + k * CL^2
    """

    wing = ac.wing
    assert wing.aspect_ratio != None and wing.psi != None and wing.phi != None, 'Aspect Ratio, Psi, Phi must be defined!'
    e = 1 / (np.pi * ac.wing.aspect_ratio * wing.psi + (1 / wing.phi))
    k = 1 / (np.pi * ac.wing.aspect_ratio * e)
    return k, e

    
def prelim_drag(ac : Aircraft,
                type_to_use : str = "Single Engine Propeller Driven",
                friction_source : str = 'lookups/skin_fric.csv',
                s_wet_source : str = 'lookups/s_wets.csv') -> float:
    
    """
    Estimate the maximum lift-to-drag ratio, (L/D)_max, using Roskam CD0 and
    Vos induced drag factor estimates.

    Parameters
    ----------
    ac : Aircraft
        Aircraft object used by the CD0 and induced drag calculations.
    type_to_use : str, optional
        Aircraft type for the Roskam wetted-area correlation.
    friction_source : str, optional
        Path to the Roskam skin-friction lookup CSV.
    s_wet_source : str, optional
        Path to the Roskam wetted-area lookup CSV.

    Returns
    -------
    float
        Estimated maximum lift-to-drag ratio.
    """

    return 0.5 * np.sqrt(1 / k(ac)[0] / cd0(ac, type_to_use, friction_source, s_wet_source))
