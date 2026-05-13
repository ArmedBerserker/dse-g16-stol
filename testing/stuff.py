import os, sys

ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..')
)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


from class1 import c1_matching_comparison as match
from classes import aircraft_2 as ac2
import numpy as np

aircraft = ac2.Aircraft('thing',
                        ac2.loader.load('concepts/reqs_nturb.yaml', ac2.Requirements),
                        ac2.loader.load('yamls/mission.yaml', ac2.Mission),
                        ac2.loader.load('yamls/weights.yaml', ac2.Weights),
                        ac2.loader.load('concepts/wing_courier.yaml', ac2.Wing),
                        ac2.loader.load('yamls/fuselage.yaml', ac2.Fuselage),
                        ac2.loader.load('concepts/engine_piston_b.yaml', ac2.Engine))

df = match.sensitivity_mesh(ac = aircraft, 
                            type_to_use = 'Twin Engine Propeller Driven',
                            W_S_plot = np.arange(0.1, 1200),
                            W_P_or_T_W_plot = np.arange(1e-8, 1e-1, 1e-4),
                            output_filepath_base = 'outputs',
                            A_values = np.arange(6, 14, 0.5),
                            CL_values = np.arange(1, 4, 0.2)
                            )

match.plot_sensitivity_mesh(df, 'outputs/workpls.png')