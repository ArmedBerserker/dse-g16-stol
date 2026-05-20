import os, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from classes import aircraft_2 as ac
from class1 import c1_m as mass
from lookups import consts as _CONST

def explore_concept_masses(parent_file : str) -> tuple[float, float]:
    
    concept = ac.loader.load(parent_file, ac.Aircraft)

    oem_eng_type = concept.engine.alpha_p_id
    config = concept.wing.type
    if config == 'A':
        gear = 'tricycle'
    elif config == 'B':
        gear = 'taildragger'
    else:
        raise ValueError('wing config is incorrect')
    oem_frac = mass.operating_empty_frac(ac = concept,
                                         source_for_fracs = 'general',
                                         engine_type = oem_eng_type,
                                         gear_type = gear)

    print(oem_frac)