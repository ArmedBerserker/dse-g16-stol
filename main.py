from class1 import c1_m, c1_loading_and_empennage, matching_diagram, prelim_drag, c1_wing_planform, c2_drag, c1_landing_gear, c1_fuselage
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
        - Preliminary drag estimation       DONE
        - Class I mass
        - Matching diagram                  DONE
        - Wing planform                     ADD FILES
        - HLD and ailerons
        - Fuselage                          DONE
        - Empennage sizing                  DONE
        - Landing gear sizing               ADD FILES
        - Class II drag
        - Class II weight                   DONE
        - Loading diagram                   DONE
        - Scissor plot????
        '''
    
    # 1. Preliminary drag
    DRAG_KWARGS = {
        'type_to_use': 'Twin Engine Propeller Driven',
        'friction_source': 'lookups/skin_fric.csv',
        's_wet_source': 'lookups/s_wets.csv',
    }
    ac.wing.CD0 = prelim_drag.cd0(ac, **DRAG_KWARGS)
    ac.wing.k, ac.wing.e = prelim_drag.k(ac)
    ac.wing.ld = prelim_drag(ac, **DRAG_KWARGS)
    print(f'\n {ac.name}: \t Preliminary drag estimation complete')

    # 2. Class I weight estimation

    print(f'\n {ac.name}: \t Class I mass estimation complete')

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
    print(f'\n {ac.name}: \t Matching complete')

    # 4. Wing planform NOTE: edit airfoil name!!!
    c1_wing_planform.size_wing_planform(ac)
    print(f'\n {ac.name}: \t Planform sizing complete')

    # 5. HLD and ailerons

    print(f'\n {ac.name}: \t HLD and aileron sizing complete')

    # 6. Fuselage
    c1_fuselage.size_fuselage(ac)
    print(f'\n {ac.name}: \t Fuselage sizing complete')

    # 7. Empennage
    tricycle_condition = ac.name.endswith('tricycle')
    c1_loading_and_empennage.class_I_loading_cgs(ac, tricycle_condition, update_ac=True)
    c1_loading_and_empennage.size_empennage_planform(ac)
    print(f'\n {ac.name}: \t Initial loading and empennage sizing complete')

    # 8. Landing gear NOTE: add .py and .yaml files
    c1_landing_gear.size_tires(ac)
    c1_landing_gear.tire_location(ac)
    print(f'\n {ac.name}: \t Landing gear sizing complete')

    # 9. Class II drag

    # 10. Class II weight
    pie_chart_folder = Path('outputs/Class_2_mass')
    mtow_pie_chart = pie_chart_folder / 'mtow_pie_chart_init.png'
    oew_pie_chart = pie_chart_folder / 'oew_pie_chart_init.png'
    struc_pie_chart = pie_chart_folder / 'struc_pie_chart_init.png'
    x_le_ht = 
    x_le_vt = 
    m_ff = 
    m_res =
    c2_m.W_oe_and_cg_from_nose(ac, ac.wing.x_le, x_le_ht, x_le_vt, update_ac=True, pie_chart_output_path=oew_pie_chart, show_pie_chart=True, struc_pie_chart_output_path=struc_pie_chart, struc_show_pie_chart=True)
    c2_m.W_to_new(ac, ac.wing.x_le, x_le_ht, x_le_vt, m_ff, m_res, update_ac=True, pie_chart_output_path=mtow_pie_chart, show_pie_chart=True)

    # 11. Loading diagram
    loading_diagram_path = pie_chart_folder / 'Initial_loading_diagram.png'
    c2_m.loading_diagram(ac.wing.x_le, ac, show_plot=True, output_filepath=loading_diagram_path, update_ac_cgs=True)

