import class1.prelim_drag as drag1
from classes.aircraft_2 import loader, Aircraft
from lookups.consts import *
import pandas as pd
import numpy as np

class Class1:
    def __init__(self, ac : Aircraft, frac_source : str = 'lookups/fuel_fracs1.csv'):
        self.ac = ac
        if ac.requirements.general['standard_type'] == None:
            print('No standard aircraft type found. Using Homebuilt.')
            self.ac_type = 'Homebuilt'
        else:
            self.ac_type = ac.general.requirements['standard_type']

        self.data = pd.read_csv(frac_source)
        self.rel_data = self.data[self.data['Airplane Type'] == self.ac_type].iloc[0]
    
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
    def _fuel_frac(self) -> float:
        engine_type = self.ac.engine.engine_type
        if engine_type == 'prop':
            return self.breguet_prop(self)
        else:
            return self.breguet_jet(self)
        
    def breguet_prop(self):
        R = self.ac.mission.range / MIL_TO_KM
        eta_prop = self.ac.engine.eta_prop
        c_p = self.ac.engine.c_p # assumed in lbs/hr/hp
        lnfrac = R / 375 * c_p / eta_prop / self.ac.wing.ld
        fracinv = np.exp(lnfrac)
        return 1 / fracinv
    
    def breguet_jet(self):
        R = self.ac.mission.range / NMIL_TO_KM
        V = self.ac.mission.cruise_speed  # as long as v_cr in kts
        c_j = self.ac.engine.c_j
        lnfrac = R * c_j / V / self.ac.wing.ld
        fracinv = np.exp(lnfrac)
        return 1/fracinv
    
    @property
    def fuel_frac_used(self):
        mulitple = self._warm_up_frac * self._taxi_frac * self._take_off_frac * self._climb_frac * self._fuel_frac * self._descent_frac * self._landing_frac
        return (1 - mulitple) * self.ac.weights.m_takeoff + self.ac.requirements.general['reserve']
    


    

    

        
        