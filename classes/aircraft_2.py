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
import os
import math as m

T = TypeVar('T')

class loader:
    '''Enables loading any kind of file into any class easily'''
    def __init__(self, filepath : str):
        self.filepath = filepath
    
    @classmethod
    def load(cls, filepath : str, target_class : Type[T]) -> T:
        return cls(filepath).instload(target_class)
    
    def instload(self, target_class : Type[T]) -> T:
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

        base_dir = os.path.dirname(self.filepath)
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(self.filepath)))

        for k, v in d.items():
            if isinstance(v, str) and v.endswith(('yaml', 'yml')):
                full_path = os.path.join(base_dir, v)
                # print(f'{k} is being loaded from {v}')
                with open(full_path, 'r') as w:
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
    general: dict
    take_off : dict
    climb : dict
    cruise : dict
    landing : dict
    approach : dict

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
    m_cargo : float | None
    m_pax : float | None
    m_energy : float | dict[float]
    m_fuel: float | None
    m_battery: float | None
    m_piston: float | None 
    m_supercap: float | None 
    m_turboprop: float | None
    m_propeller: float | None
    oew_frac: float | None
    x_cg_aft: float | None
    x_cg_fwd: float | None
    z_cg: float | None

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
    e: float | None = None
    CD0: float | None = None
    k: float | None = None
    ld : float | None = None
    wing_type: float | None = None
    c_root: float | None = None
    c_tip: float | None = None
    sweep_LE_deg: float | None = None
    sweep_TE_deg: float | None = None
    sweep_c_2_deg: float | None = None
    MAC: float | None = None
    y_MAC: float | None = None
    x_c_front_spar: float | None = None
    x_c_rear_spar: float | None = None
    x_le: float | None = None
    tip_twist: float | None = None
    dihedral: float | None = None
    airfoil_name: str | None = None
    x_c_t_c_max: float | None = None
    incidence_deg: float | None = None

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
    over_nose_angle: float | None = None 
    wall_thickness: float | None = None

    floor_thickness: float | None = None
    aisle_height: float | None = None
    top_compartment_height: float | None = None

    clearance: float | None = None
    arm_rest_width: float | None = None
    seat_width: float | None = None
    aisle_width: float | None = None
    top_compartment_width: float | None = None

    cargo_width: float | None = None
    cargo_height: float | None = None
    cargo_length: float | None = None

    tail_cone_fuselage_ratio: float | None = None
    approach_angle: float | None = None

    seat_pitch: float | None = None
    cockpit_length: float | None = None

    tail_upsweep_taildragger: float | None = None
    tail_upsweep_tricycle: float | None = None

    nose_cone_length: float | None = None
    window_angle: float | None = None

    emergency_exit_width: float | None = None
    emergency_exit_height: float | None = None

    door_width: float | None = None
    door_height: float | None = None

    length: float | None = None
    height: float | None = None
    width: float | None = None
    eq_diameter: float | None = None
    base_area: float | None = None
    tail_cone_length: float | None = None
    upsweep_angle: float | None = None
    max_cross_section_area: float | None = None
    max_perimeter: float | None = None
    start_cabin: float | None = None
    l_cabin: float | None = None
    vol_cabin_and_cargo: float | None = None
    x_pos_seats: list | None = None
    x_cargo_holds: float | None = None # list if multiple, if multiple uncomment the line below for lift of mass fractions (sum to 1)
    #mass_frac_cargo_holds: null
    n_pax: int | None = None
    n_window_seats: int | None = None
    n_middle_seats: int | None = None
    n_aisle_seats: int | None = None

    def __str__(self):
        text = "The fuselage is:\n"
        for field_info in fields(self):
            field_name = field_info.name
            field_value = getattr(self, field_name)

            text += f'{field_name}: {field_value} \n'
        return text

@dataclass
class Empennage:
    cg_and_positioning: dict
    horizontal_tail : dict
    vertical_tail : dict
    t_tail_condition : bool

    def __str__(self):
        text = "The empennage is:\n"
        for field_info in fields(self):
            field_name = field_info.name
            field_value = getattr(self, field_name)

            text += f'{field_name}: {field_value} \n'
        return text

@dataclass
class HLD_and_AIL:
    flaps: dict
    ailerons: dict
    slats: dict
    landing_lift: dict
    take_off_lift: dict
    clean_lift: dict

    def __str__(self):
        text = "The HLD and ailerons are:\n"
        for field_info in fields(self):
            field_name = field_info.name
            field_value = getattr(self, field_name)

            text += f'{field_name}: {field_value} \n'
        return text
    
@dataclass
class Engine:
    engine_type : str | None
    alpha_p_id : str | None
    count : int | None
    eta_1 : list[float] | None        # this corresponds to the fuel always
    eta_2 : list[float] | None        # this corresponds to the battery always
    eta_3 : list[float] | None        # this corresponds to the prop always. Prop is always last
    e_1 : float | None
    e_2 : float | None
    Phi : float | None
    eng_vdist_from_wing_y_c : float | None
    eng_above_wing : bool
    eng_y_pos_fuselage: float | None
    eng_x_pos: str | None
    n_fuel_tanks: int | None
    x_cg_fuel_tanks_c_r: float | None
    fuel_type: str | None
    engine_power_cruise: float | None
    super_cap_power: float | None
    power_cr: float | None
    power_to: float | None
    prop_diameter: float | None
    nac_diameter: float | None
    length_nac: float | None
    nac_t_c_max: float | None
    nac_x_c_t_c_max: float | None
    i_n: float | None
    SHP_max: float | None
    n_prop_blades: int | None

    def __post_init__(self):
        self.eta_1 = m.prod(self.eta_1)
        self.eta_2 = m.prod(self.eta_2)
        self.eta_prop = self.eta_3[-1]
        self.eta_3 = m.prod(self.eta_3)
    def __str__(self):
        text = "The engine is:\n"
        for field_info in fields(self):
            field_name = field_info.name
            field_value = getattr(self, field_name)

            text += f'{field_name}: {field_value} \n'
        return text

@dataclass
class Landing_Gear:
    gear_type: str | None
    n_struts_mlg: int | None
    n_wheels_mlg: int | None
    n_wheels_nlg: int | None
    pt: float | None
    n_nlg_min_as: float | None # Minimum nose(tail) landing gear load fraction. Value from Vos
    n_nlg_max_as: float | None # Maximum nose(tail) landing gear load fraction. Value from Vos
    tipover: float | None # Tip-over angle requirement [degrees]. Value from Vos
    scrape: float | None # Scrape angle requirement [degrees]. Value from Vos
    prop_clear: float | None # propeller clearance requirement with tail wheel [m]. From CS25.925
    fus_pitch: float | None # Fuselage inclination relative to ground plane [degrees]
    turnover: float | None # Turnover angle requirement [degrees]. Value from Roskam for rough field
    bank: float | None # Bank angle clearance requirement [degrees]. Valye from Vos
    a: float | None # Fraction of total stroke between 1-g compression and full compression. Value from Vos
    s:  float | None # Total stroke of the suspension system [m]
    fus_ground_clear:  float | None # Required distance from ground to tail cone (CHOSEN VALUE, NOT CS23 f)[m]

    selected_mlg_tire: dict | None
    selected_nlg_tire: dict | None

    longitudinal_nlg: float | None
    lateral_nlg: float | None
    height_nlg: float | None

    longitudinal_mlg: float | None
    lateral_mlg: float | None
    height_mlg: float | None

    def __str__(self):
        text = "The landing gear is:\n"
        for field_info in fields(self):
            field_name = field_info.name
            field_value = getattr(self, field_name)

            text += f'{field_name}: {field_value} \n'
        return text

@dataclass
class Aircraft:
    name : str
    requirements : Requirements
    mission : Mission
    weights : Weights
    wing : Wing
    fuselage : Fuselage
    engine : Engine
    empennage: Empennage
    hld_and_ailerons: HLD_and_AIL
    landing_gear: Landing_Gear

    @classmethod
    def from_dict(cls, data : dict):
        return cls(name = data['name'],
                   requirements = Requirements(**data['requirements']),
                   mission = Mission(**data['mission']),
                   weights = Weights(**data['weights']),
                   wing = Wing(**data['wing']),
                   fuselage = Fuselage(**data['fuselage']),
                   engine = Engine(**data['engine']))
    # def __str__(self):
    #     text = "The aircraft is:\n"
    #     for field_info in fields(self):
    #         field_name = field_info.name
    #         field_value = getattr(self, field_name)
    #         stripped_f_val = str(field_value).split('\n', 1)[1]
    #         text += f'{field_name}: {stripped_f_val} \n'
    #     return text
    def __str__(self):
        text = "The aircraft is:\n"
        for field_info in fields(self):
            field_name = field_info.name
            field_value = getattr(self, field_name)
            parts = str(field_value).split('\n', 1)
            stripped_f_val = parts[1] if len(parts) > 1 else parts[0]
            text += f'{field_name}: {stripped_f_val} \n'
        return text