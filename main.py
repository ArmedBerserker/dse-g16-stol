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
    data_cr = matching_diagram.plot_matching_and_select_design_point(ac,W_P_plot=np.arange(0.00000001,0.15,0.0001), W_S_plot=np.arange(1,1250), show_plot=False, requirement_to_meet='cruise')
    data_to = matching_diagram.plot_matching_and_select_design_point(ac,W_P_plot=np.arange(0.00000001,0.15,0.0001), W_S_plot=np.arange(1,1250), show_plot=False, requirement_to_meet='to')

    print(f" \n Aircraft: {ac.name}:")
    print(f" \n Cruise data: \n {data_cr}")
    print(f" \n Take-off data: \n {data_to}")
