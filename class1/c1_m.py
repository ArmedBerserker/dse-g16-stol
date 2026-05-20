"""
Class I aircraft sizing utilities.

Provides preliminary estimates for mission energy fractions and operating empty
mass fraction for propeller, battery-electric, and hybrid-electric aircraft.

Main functions
--------------
energy_frac_needed
    Selects the appropriate energy fraction method based on engine type.

breguet_prop
    Estimates fuel fraction for a conventional propeller aircraft.

breguet_bat
    Estimates battery fraction for a fully electric aircraft.

breguet_hyb
    Estimates fuel and battery fractions for a hybrid-electric aircraft.

operating_empty_frac
    Estimates operating empty mass fraction from either reference-aircraft data
    or statistical mass fraction coefficients.

Notes
-----
Inputs are expected to come from an `Aircraft` object. Masses are in kg, mission
range is in km, energy densities are in J/kg, and efficiencies are dimensionless.
These methods are intended for preliminary Class I sizing only.
"""

from classes.aircraft_2 import Aircraft
from lookups.consts import *
import pandas as pd
import numpy as np
from typing import Literal


def energy_frac_needed(ac : Aircraft,
                       Phi : float = None, 
                       frac_source : str = 'lookups/fuel_fracs1.csv'):
    """
    Compute the required energy fraction for an aircraft mission.

    This function selects the appropriate Breguet-based energy fraction method
    depending on the aircraft engine type.

    Parameters
    ----------
    ac : Aircraft
        Aircraft object containing engine, wing, mission, requirement, and weight data.
    frac_source : str, optional
        Path to the CSV file containing mission fuel fraction factors for conventional
        fuel-powered aircraft. Default is 'lookups/fuel_fracs1.csv'.

    Returns
    -------
    float or tuple[float, float]
        Required energy fraction. For propeller and battery aircraft, returns a single
        fraction. For hybrid aircraft, returns a tuple containing fuel fraction and
        battery fraction.

    Raises
    ------
    AssertionError
        If the engine type and power split parameter Phi are inconsistent.
    ValueError
        If the engine type is not 'prop', 'bat', or 'hyb'.
    """
    # Check the engine type of the aircraft
    eng_type = ac.engine.engine_type
    if eng_type == 'prop':
        assert ac.engine.Phi == 0
        return breguet_prop(ac, frac_source)
    elif eng_type == 'bat':
        assert ac.engine.Phi == 1
        return breguet_bat(ac)
    elif eng_type == 'hyb':
        assert ac.engine.Phi > 0 and ac.engine.Phi < 1
        return breguet_hyb(ac, Phi)
    else:
        return ValueError('Engine type not defined correctly. Must by "prop", "bat" or "hyb"')
    
def breguet_prop(ac : Aircraft, frac_source : str) -> float | tuple[float]:
    """
    Compute the fuel fraction required for a conventional propeller aircraft.

    This method applies the Vos and Roskam combined method using Breguet range relation.\n

    Parameters
    ----------
    ac : Aircraft
        Aircraft object containing engine efficiency, fuel energy density, lift-to-drag
        ratio, range, and aircraft type data.
    frac_source : str
        Path to the CSV file containing non-cruise fuel fraction factors for different
        aircraft categories.

    Returns
    -------
    float
        Fraction of takeoff mass required as fuel for the mission.

    Notes
    -----
    If no standard aircraft type is defined in the aircraft requirements, the function
    defaults to using the 'Homebuilt' aircraft type.
    """

    if ac.requirements.general['standard_type'] == None:
        print('No standard aircraft type found. Using Homebuilt.')
        ac_type = 'Homebuilt'
    else:
        ac_type = ac.requirements.general['standard_type']


    data = pd.read_csv(frac_source)
    rel_data = data[data['Airplane Type'] == ac_type].iloc[0]

    eta_fuel = ac.engine.eta_1
    eta_prop = ac.engine.eta_prop
    e_f = ac.engine.e_1
    ld = ac.wing.ld
    R = ac.mission.range * 1000 # convert to meters
    efg = e_f / g
    lnfrac = R / (eta_prop * eta_fuel * ld * efg)
    cruise_frac = np.exp(-lnfrac)
    fuel_frac = 1
    # for keys, vals in rel_data.items():
    #     if keys == 'Airplane Type':
    #         continue
    #     fuel_frac *= float(vals)
    fuel_frac *= cruise_frac

    return (1 - fuel_frac, )

def breguet_bat(ac : Aircraft) -> float:
    """
    Compute the battery mass fraction required for a fully electric aircraft.

    This method applies the Vos method for electric form of the Breguet range equation using
    battery efficiency, propulsive efficiency, battery specific energy, lift-to-drag
    ratio, and mission range.

    Parameters
    ----------
    ac : Aircraft
        Aircraft object containing battery energy density, efficiencies, lift-to-drag
        ratio, and mission range.

    Returns
    -------
    float
        Required battery mass fraction relative to takeoff mass.
    """
    eta_bat = ac.engine.eta_2
    eta_prop = ac.engine.eta_prop
    e_b = ac.engine.e_2
    ld = ac.wing.ld
    R = ac.mission.range * 1000 # convert to meters
    ebg = e_b / g
    battery_fraction = R / (eta_prop * eta_bat * ld * ebg)
    return (battery_fraction, ) 

def breguet_hyb(ac : Aircraft, Phi : float = None) -> tuple[float, float]:
    """
    Compute fuel and battery fractions for a hybrid-electric aircraft.

    This function applies the constant power split hybrid-electric range relation
    from de Vries, Hoogreef, and Vos, corresponding to Equation 17 from
    "Range Equation for Hybrid-Electric Aircraft with Constant Power Split".

    Parameters
    ----------
    ac : Aircraft
        Aircraft object containing hybrid engine efficiencies, fuel and battery
        energy densities, power split ratio, lift-to-drag ratio, and mission range.

    Returns
    -------
    tuple[float, float]
        Tuple containing:

        - fuel_frac : float
            Required fuel mass fraction.
        - battery_frac : float
            Required battery mass fraction.

    Notes
    -----
    The power split parameter Phi should satisfy 0 < Phi < 1 for hybrid aircraft.
    """
    eta_1 = ac.engine.eta_1
    eta_2 = ac.engine.eta_2
    eta_3 = ac.engine.eta_3

    e_1 = ac.engine.e_1
    e_2 = ac.engine.e_2
    if Phi == None:
        phi = ac.engine.Phi
    else:
        phi = Phi

    ld = ac.wing.ld
    R = ac.mission.range * 1000

    lnfrac = R / (eta_3 * (e_1 / g) * ld * (eta_1 + eta_2 * (phi / (1 - phi))))
    fuel_frac = 1 - np.exp(-lnfrac)

    battery_frac = (phi / (1 - phi)) * fuel_frac * (e_1 / e_2)
    return float(fuel_frac), float(battery_frac)

def operating_empty_frac(ac : Aircraft, 
                         correction : float = 1,
                         source_for_fracs : str = 'general',
                         engine_type : str = Literal['turboprop', 'piston'],
                         gear_type : str = Literal['tricycle', 'tail dragger']):
    
    if source_for_fracs == 'general':
        # These coefficients come from a linear fit of our reference aircraft
        a, b = 0.5873688639193038, -2.683435417799367
        OEW = (a * ac.weights.m_takeoff + b) * correction # correction is for our innovations
        return OEW / ac.weights.m_takeoff
    elif source_for_fracs == 'specific':
        if engine_type == 'turboprop':
            # Turboprop regression:
            a1 = 0.659427
            b1 = -426.832379 
            # R^2 = 0.980083
            #m_oe = 0.659427 m_to + -426.832379
        elif engine_type == 'piston':
            # Piston regression:
            a1 = 0.641828
            b1 = -52.572671 
            # R^2 = 0.921436
            # m_oe = 0.641828 m_to + -52.572671
        else:
            pass
        
        if gear_type == 'tail dragger':
            # Tail dragger regression:
            a2 = 0.498280
            b2 = 121.013331 
            # R^2 = 0.936882
            # m_oe = 0.498280 m_to + 121.013331
        elif gear_type == 'tricycle':
            #Tricycle regression:
            a2 = 0.585847
            b2 = 41.994669  
            # R^2 = 0.947668
            #m_oe = 0.585847 m_to + 41.994669
        else:
            pass

        OEW1 = (a1 * ac.weights.m_takeoff + b1) * correction # correction is for our innovations
        OEW2 = (a2 * ac.weights.m_takeoff + b2) * correction # correction is for our innovations
        return 0.5 * (OEW1+OEW2) / ac.weights.m_takeoff
    else:
        return -1
