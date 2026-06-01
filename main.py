from class1 import c1_m, c1_loading_and_empennage, matching_diagram, prelim_drag, c1_planform_sizing, c2_drag, c1_gear_sizing, c1_fuselage, c2_drag_new, c1_lift_and_ailerons
import c2_m
import sys
import os

# Fix path FIRST, before any local imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from classes.aircraft_2 import Aircraft, loader, Requirements, Mission, Fuselage, Wing, Engine, Weights, Empennage, HLD_and_AIL, Landing_Gear
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
                loader.load('concepts/fus_td.yaml', Fuselage),
                loader.load('concepts/engine_piston_b.yaml', Engine),
                loader.load('concepts/taildragger_emp.yaml', Empennage),
                loader.load('yamls/HLD_and_ailerons.yaml', HLD_and_AIL),
                loader.load('concepts/taildragger_gear.yaml', Landing_Gear))
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
                loader.load('concepts/fus_td.yaml', Fuselage),
                loader.load('concepts/engine_tprop_b.yaml', Engine),
                loader.load('concepts/taildragger_emp.yaml', Empennage),
                loader.load('yamls/HLD_and_ailerons.yaml', HLD_and_AIL),
                loader.load('concepts/taildragger_gear.yaml', Landing_Gear))
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
                loader.load('concepts/fus_tri.yaml', Fuselage),
                loader.load('concepts/engine_piston_b.yaml', Engine),
                loader.load('concepts/tricycle_empennage.yaml', Empennage),
                loader.load('yamls/HLD_and_ailerons.yaml', HLD_and_AIL),
                loader.load('concepts/tricycle_gear.yaml', Landing_Gear))
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
                loader.load('concepts/fus_tri.yaml', Fuselage),
                loader.load('concepts/engine_tprop_b.yaml', Engine),
                loader.load('concepts/tricycle_empennage.yaml', Empennage),
                loader.load('yamls/HLD_and_ailerons.yaml', HLD_and_AIL),
                loader.load('concepts/tricycle_gear.yaml', Landing_Gear))
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
acs = [ac1, ac3, ac6, ac8]
for ac in [ac8]:
    ''' STEPS:
        - Preliminary drag estimation       DONE
        - Class I mass                      DONE
        - Matching diagram                  DONE
        - Wing planform                     DONE
        - HLD and ailerons                  DONE
        - Fuselage                          DONE
        - Empennage sizing                  DONE
        - Landing gear sizing               DONE
        - Class II drag                     DONE
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
    ac.wing.ld = prelim_drag.prelim_drag(ac, **DRAG_KWARGS)
    print(f'\n {ac.name}: \t Preliminary drag estimation complete')

    # 2. Class I weight estimation
    c1_m.energy_frac_needed(ac, update_ac=True)
    c1_m.operating_empty_frac(ac, correction=1, source_for_fracs='specific', engine_type=ac.engine.alpha_p_id, gear_type=ac.landing_gear.gear_type, update_ac=True)
    print(f'Class I estimates: \n empty: {ac.weights.m_empty} \n frac: {ac.weights.oew_frac}')
    print(f'\n {ac.name}: \t Class I mass estimation complete')

    # 3. Matching Diagram
    output_matching = f'outputs/matching_init_{ac.name}.png'
    data_cr = matching_diagram.plot_matching_and_select_design_point(ac, output_filepath=output_matching, W_P_plot=np.arange(0.00000001,0.15,0.0001), W_S_plot=np.arange(1,1250), show_plot=False, requirement_to_meet='cruise')
    data_to = matching_diagram.plot_matching_and_select_design_point(ac, W_P_plot=np.arange(0.00000001,0.15,0.0001), W_S_plot=np.arange(1,1250), show_plot=False, requirement_to_meet='to')
    ac.engine.power_to = ac.weights.m_takeoff * g / data_to['W/P']
    ac.wing.area = ac.weights.m_takeoff * g / data_to['W/S']
    ac.engine.power_cr = ac.weights.m_takeoff * g / data_cr['W/P']
    # print(f" \n Aircraft: {ac.name}:")
    # print(f" \n Cruise data: \n {data_cr}")
    # print(f" \n Take-off data: \n {data_to}")
    print(f'\n {ac.name}: \t Matching complete')

    # 4. Wing planform NOTE: edit airfoil name!!!
    c1_planform_sizing.size_wing_planform(ac)
    print(f'\n {ac.name}: \t Planform sizing complete')

    # 5. HLD and ailerons
    c1_lift_and_ailerons.size_ailerons(ac, update_ac=True)
    c1_lift_and_ailerons.size_HLD(ac, update_ac=True)
    print(f'\n {ac.name}: \t HLD and aileron sizing complete')

    # # 6. Fuselage
    # c1_fuselage.size_fuselage(ac)
    # print(f'\n {ac.name}: \t Fuselage sizing complete')

    # 7. Empennage
    tricycle_condition = ac.landing_gear.gear_type == 'tricycle'
    # tricycle_condition = ac.name.endswith('tricycle')
    # c1_loading_and_empennage.class_I_loading_cgs(ac, tricycle_condition, update_ac=True)
    c1_loading_and_empennage.classI_loading_and_cgs_2(ac, update_ac=True)
    c1_loading_and_empennage.size_empennage_planform(ac)
    print(f'\n {ac.name}: \t Initial loading and empennage sizing complete')
    # print(f'Empennage: {ac.empennage}')

    # 8. Landing gear NOTE: add .py and .yaml files
    c1_gear_sizing.size_tires(ac, update_ac=True)
    c1_gear_sizing.tire_location(ac, update_ac=True)
    print(f'\n {ac.name}: \t Landing gear sizing complete')

    # 9. Class II drag
    c2_drag_new.CD0(ac, n_engine_inoperative=ac.engine.count, flight_condition='cruise', update_ac=True)
    c2_drag_new.CD0(ac, n_engine_inoperative=ac.engine.count, flight_condition='take-off', update_ac=False)
    c2_drag_new.CD0(ac, n_engine_inoperative=ac.engine.count, flight_condition='landing', update_ac=False)
    c2_drag_new.C_D_L(ac, CD0=ac.wing.CD0, flight_condition='cruise', update_ac=True, wing_tip=False)
    print(f'\n {ac.name} done (except loading diagram and Class II mass estim)')

    # 10. Class II weight
    pie_chart_folder = Path('outputs/Class_2_mass')
    mtow_pie_chart = pie_chart_folder / f'mtow_pie_chart_init_{ac.name}.png'
    oew_pie_chart = pie_chart_folder / f'oew_pie_chart_init_{ac.name}.png'
    struc_pie_chart = pie_chart_folder / f'struc_pie_chart_init_{ac.name}.png'
    x_le_ht = ac.empennage.horizontal_tail['x_le']
    x_le_vt = ac.empennage.vertical_tail['x_le']
    c2_m.W_oe_and_cg_from_nose(ac, update_ac=True, pie_chart_output_path=oew_pie_chart, show_pie_chart=False, struc_pie_chart_output_path=struc_pie_chart, struc_show_pie_chart=False)
    print(ac.weights.m_empty)
    c2_m.W_to_new(ac, W_crew=0, update_ac=True, pie_chart_output_path=mtow_pie_chart, show_pie_chart=False)
    print(ac.weights.m_takeoff)
    # print(ac)
    print(f'\n Empty mass fraction: {ac.weights.oew_frac} \n ')
    print(f'Class II estimates: \n empty: {ac.weights.m_empty} \n mtow: {ac.weights.m_takeoff}')
    # print(ac.wing.CD0)

    c1_m.energy_frac_needed(ac, update_ac=True)
    # c1_m.operating_empty_frac(ac, correction=1, source_for_fracs='specific', engine_type=ac.engine.alpha_p_id, gear_type=ac.landing_gear.gear_type, update_ac=True)
    print(f'\n {ac.name}: \t Class I mass estimation complete')
    print(f' oew: {ac.weights.m_empty}, mtow: {ac.weights.m_takeoff}, oew frac: {ac.weights.oew_frac}, fuel frac: {ac.weights.m_fuel / ac.weights.m_takeoff}')

    # print(ac.weights)
    # print(ac.landing_gear)
    # print(ac)

    # 11. Loading diagram
    # loading_diagram_path = pie_chart_folder / 'Initial_loading_diagram.png'
    # c2_m.loading_diagram(ac.wing.x_le, ac, show_plot=True, output_filepath=loading_diagram_path, update_ac_cgs=True)

