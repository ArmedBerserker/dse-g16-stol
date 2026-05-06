import class1.prelim_drag as drag1
from classes.aircraft_2 import loader, Aircraft
from lookups.consts import *
import pandas as pd
import numpy as np

class Class1:
    def __init__(self, ac : Aircraft, frac_source : str = 'lookups/fuel_fracs1.csv', mass_source : str = 'lookups/mass_relations1.csv'):
        self.ac = ac
        if ac.requirements.general['standard_type'] == None:
            print('No standard aircraft type found. Using Homebuilt.')
            self.ac_type = 'Homebuilt'
        else:
            self.ac_type = ac.requirements.general['standard_type']

        self.data = pd.read_csv(frac_source)
        self.rel_data = self.data[self.data['Airplane Type'] == self.ac_type].iloc[0]
        self.data2 = pd.read_csv(mass_source)    
    @property
    def _warm_up_frac(self) -> float:
        return self.rel_data['Warm-up']
    
    @property
    def _taxi_frac(self) -> float:
        return self.rel_data['Taxi']
    
    @property
    def _take_off_frac(self) -> float:
        return self.rel_data['Take-off']
    
    @property
    def _climb_frac(self) -> float:
        return self.rel_data['Climb']
    
    @property
    def _descent_frac(self) -> float:
        return self.rel_data['Descent']
    
    @property
    def _landing_frac(self) -> float:
        return self.rel_data['Landing Taxi']
    
    @property
    def _fuel_frac_cruise(self) -> float:
        engine_type = self.ac.engine.engine_type
        if engine_type == 'prop':
            return self.breguet_prop()
        elif engine_type == 'bat':
            return self.breguet_bat()
        elif engine_type == 'hyb':
            return self.breguet_hyb()
    
    @property
    def _fuel_frac_endurance(self) -> float:
        engine_type = self.ac.engine.engine_type
        if engine_type == 'prop':
            ...

    def breguet_prop(self):
        R = self.ac.mission.range / MIL_TO_KM
        eta_prop = self.ac.engine.eta_prop
        c_p = self.ac.engine.c_p # assumed in lbs/hr/hp
        lnfrac = R / 375 * c_p / eta_prop / self.ac.wing.ld
        fracinv = np.exp(lnfrac)
        return 1 / fracinv
    
    def breguet_(self):
        R = self.ac.mission.range / NMIL_TO_KM
        V = self.ac.mission.cruise_speed  # as long as v_cr in kts
        c_j = self.ac.engine.c_j
        lnfrac = R * c_j / V / self.ac.wing.ld
        fracinv = np.exp(lnfrac)
        return 1/fracinv
    
    def fuel_frac_used(self):
        mulitple = self._warm_up_frac * self._taxi_frac * self._take_off_frac * self._climb_frac * self._fuel_frac * self._descent_frac * self._landing_frac
        return ((1 - mulitple) * self.ac.weights.m_takeoff + self.ac.requirements.general['reserve']) / self.ac.weights.m_takeoff
     
    def struct_mass_frac(self, subtype : str = 'Propeller Driven'):
        row = self.data2[
        (self.data2["Airplane Type"] == self.ac_type) &
        (self.data2["Subtype"] == subtype)
        ]
        A = row['A'].iloc[0]
        B = row['B'].iloc[0]

        W_e = 10 ** ((np.log10(self.ac.weights.m_takeoff / LBS_TO_KG) - A) / B)
        return W_e / (self.ac.weights.m_takeoff / LBS_TO_KG)
    

def energy_frac_needed(ac : Aircraft, frac_source : str = 'lookups/fuel_fracs1.csv'):
     
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
        return breguet_hyb(ac)
    else:
        return ValueError('Engine type not defined correctly. Must by "prop", "bat" or "hyb"')
    
def breguet_prop(ac : Aircraft, frac_source : str) -> float | tuple[float]:
    '''
    Applies the Vos + Roskam methods
    '''

    if ac.requirements.general['standard_type'] == None:
        print('No standard aircraft type found. Using Homebuilt.')
        ac_type = 'Homebuilt'
    else:
        ac_type = ac.requirements.general['standard_type']


    data = pd.read_csv(frac_source)
    rel_data = data[data['Airplane Type'] == ac_type].iloc[0]

    eta_fuel = ac.engine.eta_1
    eta_prop = ac.engine.eta_2
    e_f = ac.engine.e_1
    ld = ac.wing.ld
    R = ac.mission.range * 1000 # convert to meters
    efg = e_f / g
    lnfrac = R / (eta_prop * eta_fuel * ld * efg)
    cruise_frac = np.exp(-lnfrac)
    fuel_frac = 1.
    for keys, vals in rel_data.items():
        if keys == 'Airplane Type':
            continue
        fuel_frac *= float(vals)
    fuel_frac *= cruise_frac

    return 1 - fuel_frac

def breguet_bat(ac : Aircraft) -> float:
    '''
    Applies the Vos method for energy fraction
    '''
    eta_bat = ac.engine.eta_1
    eta_prop = ac.engine.eta_2
    e_f = ac.engine.e_f
    ld = ac.wing.ld
    R = ac.mission.range * 1000 # convert to meters
    efg = e_f / g
    battery_fraction = R / (eta_prop * eta_bat * ld * efg)
    return battery_fraction 

def breguet_hyb(ac : Aircraft) -> tuple[float, float]:
    '''
    Applies the de Vries, Hoogreef, Vos method for hybrid fractions \n
    corresponding to eq 17 from *"Range Equation for Hybrid-Electric Aircraft with Constant Power Split"*
    '''

    eta_1 = ac.engine.eta_1
    eta_2 = ac.engine.eta_2
    eta_3 = ac.engine.eta_3

    e_1 = ac.engine.e_1
    e_2 = ac.engine.e_2

    Phi = ac.engine.Phi

    ld = ac.wing.ld
    R = ac.mission.range * 1000

    lnfrac = R / (eta_3 * (e_1 / g) * ld * (eta_1 + eta_2 * (Phi / (1 - Phi))))
    fuel_frac = 1 - np.exp(-lnfrac)

    battery_frac = (Phi / (1 - Phi)) * fuel_frac * (e_1 / e_2)
    return float(fuel_frac), float(battery_frac)


ac = loader.load('yamls/aircraft.yaml', Aircraft)
print(energy_frac_needed(ac))