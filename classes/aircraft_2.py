"""
Aircraft data structures and YAML loading utilities.

Defines the main dataclasses used to store aircraft requirements, mission,
weights, wing geometry, fuselage geometry, and the complete Aircraft object.
Also provides a simple loader class for reading YAML files into these
dataclasses.
"""

from dataclasses import dataclass, is_dataclass, fields, field
from typing import Type, TypeVar, Any
import yaml

T = TypeVar('T')

class loader:
    '''Enables loading any kind of file into any class easily'''
    def __init__(self, filepath : str):
        self.filepath = filepath
    
    def load(self, target_class : Type[T]) -> T:
        '''wrapper to call easily'''
        data = self._read_file()

        if hasattr(target_class, 'from_dict'):
            return target_class.from_dict(data)
        else:
            return self._build_dataclass(target_class, data)
    
    def _read_file(self) -> dict:
        '''reads a yaml file and returns a dictionary'''
        with open(self.filepath, 'r') as f:
            d = yaml.safe_load(f)

        for k, v in d.items():
            if isinstance(v, str) and v.endswith(('yaml', 'yml')):
                print(f'{k} is being loaded from {v}')
                with open(v, 'r') as w:
                    d[k] = yaml.safe_load(w)
            else:
                pass

        return d
    
    def _build_dataclass(self, target_class : Type[T], data : dict) -> T:
        '''builds the target class from a dictionary input'''
        if not is_dataclass(target_class): 
            raise TypeError(f'{target_class.__name__} must be a dataclass')
        
        init_values = {}

        for field_info in fields(target_class):
            field_name = field_info.name
            #field_type = field_info.type

            if field_name not in data:
                print(f'Input file {self.filepath} is missing {field_name}')
            
            init_values[field_name] = data[field_name] # THIS CURRENTLY DOES NOT ALLOW CLASSES THAT HAVE CLASS INPUTS

        return target_class(**init_values)


@dataclass
class Requirements:
    take_off : dict
    climb : dict
    cruise : dict
    landing : dict
    approach : dict
    climb_gradient : dict

    def __str__(self):
        text = "The requirements are:\n"
        for field_info in fields(self):
            field_name = field_info.name
            field_value = getattr(self, field_name)

            text += f'{field_name}: {field_value} \n'
        return text

@dataclass
class Mission:
    range : float | None
    cruise_altitude : float | None
    cruise_speed : float | None
    endurance : float | None

    def __str__(self):
        text = "The mission is:\n"
        for field_info in fields(self):
            field_name = field_info.name
            field_value = getattr(self, field_name)

            text += f'{field_name}: {field_value} \n'
        return text
    

@dataclass
class Weights:
    m_takeoff : float | None
    m_empty : float | None
    m_payload : float | None
    m_energy : float | dict[float]

    def __str__(self):
        text = "The weights are:\n"
        for field_info in fields(self):
            field_name = field_info.name
            field_value = getattr(self, field_name)

            text += f'{field_name}: {field_value} \n'
        return text

@dataclass
class Wing:
    area : float | None = None
    span : float | None = None
    aspect_ratio : float | None = None
    taper_ratio : float | None = None
    sweep : float | None = None
    c_f : float | None = None
    phi : float | None = None
    psi : float | None = None
    airfoils : list[str] = None
    # ADD WHATEVER IS NEEDED

    def __str__(self):
        text = "The wing is:\n"
        for field_info in fields(self):
            field_name = field_info.name
            field_value = getattr(self, field_name)

            text += f'{field_name}: {field_value} \n'
        return text

@dataclass
class Fuselage:
    length : float | None
    span : float
    height : float
    wetted_area : float

    def __str__(self):
        text = "The fuselage is:\n"
        for field_info in fields(self):
            field_name = field_info.name
            field_value = getattr(self, field_name)

            text += f'{field_name}: {field_value} \n'
        return text

@dataclass
class Aircraft:
    requirements : Requirements
    mission : Mission
    weights : Weights
    wing : Wing
    fuselage : Fuselage

    @classmethod
    def from_dict(cls, data : dict):
        return cls(requirements = Requirements(**data['requirements']),
                   mission = Mission(**data['mission']),
                   weights = Weights(**data['weights']),
                   wing = Wing(**data['wing']),
                   fuselage = Fuselage(**data['fuselage']))
    
    def __str__(self):
        text = "The aircraft is:\n"
        for field_info in fields(self):
            field_name = field_info.name
            field_value = getattr(self, field_name)
            stripped_f_val = str(field_value).split('\n', 1)[1]
            text += f'{field_name}: {stripped_f_val} \n'
        return text

