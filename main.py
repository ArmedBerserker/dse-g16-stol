from class1 import c1_m, c1_loading_and_empennage, matching_diagram, prelim_drag, c1_wing_planform, c2_drag, c1_landing_gear
import c2_m
import sys
import os

# Fix path FIRST, before any local imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from classes.aircraft_2 import Aircraft, loader, Requirements, Mission, Fuselage, Wing, Engine, Weights, Empennage, HLD_and_AIL
import numpy as np
import matplotlib.pyplot as plt
from classes.isa import Atmosphere
from lookups.consts import *
from pathlib import Path
import pandas as pd

# sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
BASE_DIR = Path(__file__).resolve().parent

ac1 = Aircraft('Boosted_piston_taildragger',
                loader.load('concepts/reqs_nturb.yaml', Requirements),
                loader.load('yamls/mission.yaml', Mission),
                loader.load('yamls/weights.yaml', Weights),
                loader.load('concepts/wing_courier.yaml', Wing),
                loader.load('concepts/fuselage_taildragger.yaml', Fuselage),
                loader.load('concepts/engine_piston_b.yaml', Engine),
                loader.load('yamls/empennage_config.yaml', Empennage),
                loader.load('yamls/HLD_and_ailerons.yaml', HLD_and_AIL))
# ac2 = Aircraft('Piston_hybrid_taildragger',
#                 loader.load('concepts/reqs_nturb.yaml', Requirements),
#                 loader.load('yamls/mission.yaml', Mission),
#                 loader.load('yamls/weights.yaml', Weights),
#                 loader.load('concepts/wing_courier.yaml', Wing),
#                 loader.load('concepts/fuselage_taildragger.yaml', Fuselage),
#                 loader.load('concepts/engine_piston_e.yaml', Engine))
ac3 = Aircraft('Boosted_turboprop_taildragger',
                loader.load('concepts/reqs_turb.yaml', Requirements),
                loader.load('yamls/mission.yaml', Mission),
                loader.load('yamls/weights.yaml', Weights),
                loader.load('concepts/wing_courier.yaml', Wing),
                loader.load('concepts/fuselage_taildragger.yaml', Fuselage),
                loader.load('concepts/engine_tprop_b.yaml', Engine),
                loader.load('yamls/empennage_config.yaml', Empennage),
                loader.load('yamls/HLD_and_ailerons.yaml', HLD_and_AIL))
# ac4 = Aircraft('Turbine_hybrid_taildragger',
#                 loader.load('concepts/reqs_turb.yaml', Requirements),
#                 loader.load('yamls/mission.yaml', Mission),
#                 loader.load('yamls/weights.yaml', Weights),
#                 loader.load('concepts/wing_courier.yaml', Wing),
#                 loader.load('concepts/fuselage_taildragger.yaml', Fuselage),
#                 loader.load('concepts/engine_turb_e.yaml', Engine))
# ac5 = Aircraft('H2_taildragger',
#                 loader.load('concepts/reqs_nturb.yaml', Requirements),
#                 loader.load('yamls/mission.yaml', Mission),
#                 loader.load('yamls/weights.yaml', Weights),
#                 loader.load('concepts/wing_courier.yaml', Wing),
#                 loader.load('concepts/fuselage_taildragger.yaml', Fuselage),
#                 loader.load('concepts/engine_h2.yaml', Engine))
ac6 = Aircraft('Boosted_piston_tricycle',
                loader.load('concepts/reqs_nturb.yaml', Requirements),
                loader.load('yamls/mission.yaml', Mission),
                loader.load('yamls/weights.yaml', Weights),
                loader.load('concepts/wing_electra.yaml', Wing),
                loader.load('concepts/fuselage_tricycle.yaml', Fuselage),
                loader.load('concepts/engine_piston_b.yaml', Engine),
                loader.load('yamls/empennage_config.yaml', Empennage),
                loader.load('yamls/HLD_and_ailerons.yaml', HLD_and_AIL))
# ac7 = Aircraft('Piston_hybrid_tricycle',
#                 loader.load('concepts/reqs_nturb.yaml', Requirements),
#                 loader.load('yamls/mission.yaml', Mission),
#                 loader.load('yamls/weights.yaml', Weights),
#                 loader.load('concepts/wing_electra.yaml', Wing),
#                 loader.load('concepts/fuselage_tricycle.yaml', Fuselage),
#                 loader.load('concepts/engine_piston_e.yaml', Engine))
ac8 = Aircraft('Boosted_turboprop_tricycle',
                loader.load('concepts/reqs_turb.yaml', Requirements),
                loader.load('yamls/mission.yaml', Mission),
                loader.load('yamls/weights.yaml', Weights),
                loader.load('concepts/wing_electra.yaml', Wing),
                loader.load('concepts/fuselage_tricycle.yaml', Fuselage),
                loader.load('concepts/engine_tprop_b.yaml', Engine),
                loader.load('yamls/empennage_config.yaml', Empennage),
                loader.load('yamls/HLD_and_ailerons.yaml', HLD_and_AIL))
# ac9 = Aircraft('Turbine_hybrid_tricycle',
#                 loader.load('concepts/reqs_turb.yaml', Requirements),
#                 loader.load('yamls/mission.yaml', Mission),
#                 loader.load('yamls/weights.yaml', Weights),
#                 loader.load('concepts/wing_electra.yaml', Wing),
#                 loader.load('concepts/fuselage_tricycle.yaml', Fuselage),
#                 loader.load('concepts/engine_turb_e.yaml', Engine))
# ac10 = Aircraft('H2_tricycle',
#                 loader.load('concepts/reqs_nturb.yaml', Requirements),
#                 loader.load('yamls/mission.yaml', Mission),
#                 loader.load('yamls/weights.yaml', Weights),
#                 loader.load('concepts/wing_electra.yaml', Wing),
#                 loader.load('concepts/fuselage_tricycle.yaml', Fuselage),
#                 loader.load('concepts/engine_h2.yaml', Engine))
for ac in [ac1, ac3, ac6, ac8]:
    ''' STEPS:
        - Preliminary drag estimation
        - Class I mass
        - Matching diagram
        - Wing planform
        - HLD and ailerons
        - Fuselage
        - Empennage sizing
        - Landing gear sizing
        - Class II weight
        - Class II drag'''
    
    # 1. Preliminary drag
    DRAG_KWARGS = {
        'type_to_use': 'Twin Engine Propeller Driven',
        'friction_source': 'lookups/skin_fric.csv',
        's_wet_source': 'lookups/s_wets.csv',
    }
    ac.wing.CD0 = prelim_drag.cd0(ac, **DRAG_KWARGS)
    ac.wing.k, ac.wing.e = prelim_drag.k(ac)
    ac.wing.ld = prelim_drag(ac, **DRAG_KWARGS)

    # 2. Class I weight estimation


    # 3. Matching Diagram
    data_cr = matching_diagram.plot_matching_and_select_design_point(ac,W_P_plot=np.arange(0.00000001,0.15,0.0001), W_S_plot=np.arange(1,1250), show_plot=False, requirement_to_meet='cruise')
    data_to = matching_diagram.plot_matching_and_select_design_point(ac,W_P_plot=np.arange(0.00000001,0.15,0.0001), W_S_plot=np.arange(1,1250), show_plot=False, requirement_to_meet='to')
    ac.engine.power_to = ac.weights.m_takeoff * g / data_to['W/P']
    ac.wing.area = ac.weights.m_takeoff * g / data_to['W/S']
    data_cr = matching_diagram.plot_matching_and_select_design_point(ac, type_to_use='Twin Engine Propeller Driven', W_P_plot=np.arange(0.00000001,0.15,0.0001), W_S_plot=np.arange(1,1250), output_filepath='outputs/Iteration_matching_plot.png', requirement_to_meet='cruise')
    ac.engine.power_cr = ac.weights.m_takeoff * g / data_cr['W/P']
    # print(f" \n Aircraft: {ac.name}:")
    # print(f" \n Cruise data: \n {data_cr}")
    # print(f" \n Take-off data: \n {data_to}")

    # 4. Wing planform NOTE: edit airfoil name!!!
    c1_wing_planform.size_wing_planform(ac)

    # 5. HLD and ailerons

    # 6. Fuselage
