import json
from dataclasses import dataclass, field
import pandas as pd
from typing import Literal

# some general variables
engine_types = Literal['jet', 'prop']

# classes start
@dataclass
class Wing:
    '''This class should define the basics of a wing'''
    span : float                            # span of the wing [m]
    wing_area : float                       # area of the wing [m2]
    airfoils : str | list[str]              # airfoil(s), list order must be root, tip
    aspect_ratio : float | None = None      # aspect ratio (can be computed too)
    taper_ratio : float | None = None       # c_tip / c_root (can be left if unknown)
    def __post_init__(self):
        if self.aspect_ratio is None: self.aspect_ratio = self.span ** 2 / self.wing_area

        # Ensure the aspect ratio is not correct

    @classmethod
    def from_json(cls, file : str):
        '''Alternate method of defining a wing from a json file (see wing.json)'''
        with open(file, 'r') as f:
            data = json.load(f)
        return cls(**data)

@dataclass
class Engine:
    '''This class defines the basics of the engine'''
    engine_type : engine_types      # jet or prop
    engine_name : str               # engine name
    fuel_energy : float             # specific energy of the fuel in J/kg
    
    @classmethod
    def from_json(cls, path : str):
        with open(path, 'r') as f:
            data = json.load(f)
        if data['engine_type'] == 'jet':
            return Jet(**data)
        elif data['engine_type'] == 'prop':
            return Prop(**data)
        else:
            return cls(**data)
    
@dataclass
class Jet(Engine):
    engine_type : engine_types = field(default='jet', init=False) # ensures that the Engine class gets a jet type

    thrust : float              # The thrust of the engine in N
    TSFC : float | None         # Thrust Specific Fuel Consumption in g/kN/s
    B : float                   # Bypass ratio

    def __post_init__(self):
        if self.TSFC == None:
            self.TSFC = 22 * self.B ** -0.19
    


@dataclass
class Prop(Engine):
    engine_type : engine_types = field(default='prop', init=False) 

    power : float               # Power in hp
    PSFC : float                # Power Specific Fuel Consumption in lb/hp/hr


@dataclass
class Aircraft:
    '''This class should define the basics of an aircraft, either from a csv/datafile or manually'''
    ac_type : str               # aircraft name
    ac_config : str             # high-wing, low-wing, canard, etc
    wing : Wing                 # the main lift making surface
    engine : Engine             # the engine
    n_engine : int              # the number of these engines
    range : float | None        # range in km
    payload_mass : float | None # mass in kg


    @classmethod
    def from_json(cls, file : str, wing_from_path = True, engine_from_path = True):
        '''Alternate method of defining an aircraft from a json file (see ac.json)'''
        with open(file, 'r') as f:
            ac = json.load(f)


        if wing_from_path:
            wing_path = ac['wing_path']
            with open(wing_path, 'r') as t:
                wing_params = json.load(t)
        else:
            wing_params = ac['wing']
        
        if engine_from_path:
            engine_path = ac['engine_path']
            with open(engine_path, 'r') as e:
                engine_params = json.load(e)
        else:
            engine_params = ac['engine']

        return cls(
            ac['ac_type'],
            ac['ac_config'],
            Wing(**wing_params),
            Engine(**engine_params),
            ac['n_engine'],
            ac['range'],
            ac['payload_mass']
        )



ac = Aircraft.from_json('jsons/docs/ac.json', True, True)