# Fix path FIRST, before any local imports
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from classes.aircraft_2 import loader, Aircraft, Requirements, Mission, Weights, Wing, Fuselage, Engine
from lookups.consts import *
from class1.prelim_drag import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines
from classes.isa import Atmosphere
from class1.matching_diagram import plot_matching_and_select_design_point
from class1 import c1_m
from pathlib import Path

''' TO DO:
    - Check the curve_labels in matching diagram point selecter and consistency with plotting function and here
    -  '''


def Resize_A(limit_W_P):
    resize = False
    if limit_W_P in ["cruise speed",
                     "aeo roc",
                     "aeo climb gradient",
                     "aeo climb gradient (turbine)",
                     "balked landing",
                     "balked landing (turbine)",
                     "oei roc/climb gradient I (turbine)",
                     "oei roc/climb gradient II (turbine)"]:
        resize = True
    return resize
    # elif limit_W_P == :  # NOTE: fill in names of lines and check it runs
    #     resize = False
    #     return resize

def Resize_CL_max_LD(limit_W_S):
    resize = False
    if limit_W_S in ["stall speed", "landing field length"]:
        resize = True
    return resize

def max_RC_and_Climb_grad(W_S, W_P, k, rho, eta_p, CD0):
    V_RC = np.sqrt(W_S / (0.5 * rho) * np.sqrt(k / (3 * CD0)))
    V_CG = 4 * W_S * k / (rho * eta_p / W_P)

    RC = eta_p / W_P - V_RC * (0.5 * rho * V_RC ** 2 / W_S * CD0 + W_S * k / (0.5 * rho * V_RC ** 2))
    CG = np.arcsin(eta_p / V_CG / W_P - 0.5 * rho * V_CG ** 2 / W_S * CD0 - W_S * k / (0.5 * rho * V_CG ** 2))
    return RC, np.rad2deg(CG)

def sensitivity_study(
        ac: Aircraft,
        type_to_use: str,
        W_S_plot: np.ndarray,
        W_P_or_T_W_plot: np.ndarray,
        output_filepath_base: str,
        param: str,                         # 'CL_max_LD' or 'A'
        step: float,
        n_steps: int,
        initial_limiting_wp_constraint: str, # from the baseline run
        initial_limiting_ws_constraint: str,
        initial_W_P: float,
        initial_W_S: float,
) -> dict:
    """
    Runs a sensitivity study on either CL_max_LD or aspect ratio (A),
    incrementing by `step` up to `n_steps` times.

    Only the target parameter is modified per call. The other is held
    fixed at its current value on ac and restored via finally.

    Returns a dict of lists, one entry per successful iteration.
    """
    # Snapshot both values so we can guarantee neither drifts
    initial_CL_LD = ac.requirements.landing['as_CL_max_la']
    initial_A     = ac.wing.aspect_ratio

    if param == 'CL_max_LD':
        should_resize = lambda wp, ws: Resize_CL_max_LD(ws)
        def apply_step(i):
            ac.requirements.landing['as_CL_max_la'] = initial_CL_LD + i * step
            ac.wing.aspect_ratio = initial_A          # explicitly hold A fixed
        def restore():
            ac.requirements.landing['as_CL_max_la'] = initial_CL_LD
            ac.wing.aspect_ratio = initial_A          # restore both
        def get_CL_LD(i): return initial_CL_LD + i * step
        def get_A(_):     return initial_A

    elif param == 'A':
        should_resize = lambda wp, ws: Resize_A(wp)
        def apply_step(i):
            ac.wing.aspect_ratio = initial_A + i * step
            ac.requirements.landing['as_CL_max_la'] = initial_CL_LD  # explicitly hold CL_LD fixed
        def restore():
            ac.wing.aspect_ratio = initial_A
            ac.requirements.landing['as_CL_max_la'] = initial_CL_LD  # restore both
        def get_CL_LD(_): return initial_CL_LD
        def get_A(i):     return initial_A + i * step

    else:
        raise ValueError(f"param must be 'CL_max_LD' or 'A', got '{param}'")

    # Output history lists 
    W_P_history                    = [initial_W_P]
    W_S_history                    = [initial_W_S]
    CL_max_LD_history              = [initial_CL_LD]
    A_history                      = [initial_A]
    limiting_wp_constraint_history = [initial_limiting_wp_constraint]
    limiting_ws_constraint_history = [initial_limiting_ws_constraint]

    # Seed the constraint check with the baseline values
    last_wp_constraint = initial_limiting_wp_constraint
    last_ws_constraint = initial_limiting_ws_constraint

    try:
        for i in range(1, n_steps + 1):
            if not should_resize(last_wp_constraint, last_ws_constraint):
                break  # Constraint no longer limiting — stop early

            apply_step(i)

            output_filepath_i = f"{output_filepath_base}_r{i}_{param}.png"
            data_i = plot_matching_and_select_design_point(ac, 
                                                           type_to_use,
                                                           W_S_plot, 
                                                           W_P_or_T_W_plot, 
                                                           output_filepath_i,
                                                           show_plot=False)
            plt.close('all')

            W_P_history.append(data_i['W/P'])
            W_S_history.append(data_i['W/S'])
            CL_max_LD_history.append(get_CL_LD(i))
            A_history.append(get_A(i))
            limiting_wp_constraint_history.append(data_i['limiting_wp_constraint'])
            limiting_ws_constraint_history.append(data_i['limiting_ws_constraint'])

            # Update for next iteration's resize check
            last_wp_constraint = data_i['limiting_wp_constraint']
            last_ws_constraint = data_i['limiting_ws_constraint']

    finally:
        restore()  # Always reset both parameters

    return {
        "W/P":                     W_P_history,
        "W/S":                     W_S_history,
        "CL_max_LD":               CL_max_LD_history,
        "A":                       A_history,
        "limiting_wp_constraint":  limiting_wp_constraint_history,
        "limiting_ws_constraint":  limiting_ws_constraint_history,
    }

def sensitivity_mesh(
        ac: Aircraft,
        type_to_use: str,
        W_S_plot: np.ndarray,
        W_P_or_T_W_plot: np.ndarray,
        output_filepath_base: str,
        A_values: np.ndarray,
        CL_values: np.ndarray,
):
    """
    Computes a full 2D sensitivity mesh over:
        - aspect ratio A
        - CL_max_LD

    Returns a dataframe with:
        A, CL_max_LD, W/S, W/P
    """

    initial_A = ac.wing.aspect_ratio
    initial_CL = ac.requirements.landing['as_CL_max_la']

    rows = []

    try:

        for A in A_values:

            ac.wing.aspect_ratio = A

            for CL in CL_values:

                ac.requirements.landing['as_CL_max_la'] = CL

                data = plot_matching_and_select_design_point(
                    ac,
                    type_to_use,
                    W_S_plot,
                    W_P_or_T_W_plot,
                    output_filepath=None,
                    show_plot=False
                )

                rows.append({
                    'A': A,
                    'CL_max_LD': CL,
                    'W/S': data['W/S'],
                    'W/P': data['W/P'],
                    'limiting_wp_constraint': data['limiting_wp_constraint'],
                    'limiting_ws_constraint': data['limiting_ws_constraint'],
                })

                plt.close('all')

    finally:
        # restore original values
        ac.wing.aspect_ratio = initial_A
        ac.requirements.landing['as_CL_max_la'] = initial_CL

    return pd.DataFrame(rows)

def plot_sensitivity_mesh(df, output_filepath):

    fig, ax = plt.subplots(figsize=(10,8))

    # -------------------------------------------------
    # CONSTANT A CURVES
    # -------------------------------------------------

    unique_A = sorted(df['A'].unique())

    for A in unique_A:

        sub = df[df['A'] == A]

        sub = sub.sort_values('CL_max_LD')

        ax.plot(
            sub['W/S'],
            sub['W/P'],
            '--o',
            linewidth=1.5,
            markersize=4,
            label=f'A={A:.1f}'
        )

        # label at end
        last = sub.iloc[-6]

        ax.text(
            last['W/S']-10,
            last['W/P'],
            f'A={A:.1f}',
            fontsize=8
        )

    # -------------------------------------------------
    # CONSTANT CL CURVES
    # -------------------------------------------------

    unique_CL = sorted(df['CL_max_LD'].unique())

    for CL in unique_CL:

        sub = df[df['CL_max_LD'] == CL]

        sub = sub.sort_values('A')

        ax.plot(
            sub['W/S'],
            sub['W/P'],
            '--',
            linewidth=1.2,
            alpha=0.8
        )

        last = sub.iloc[1]

        ax.text(
            last['W/S'],
            last['W/P'],
            f'$C_{{max_{{LD}}}}$={CL:.2f}',
            fontsize=8
        )

    # W_P_ref = [0.054147611,0.075714426,0.060130048,0.059649008,0.067577078,0.057702788,0.068512198,0.079514477,0.071580563,0.070699017,0.063456392]
    # W_S_ref = [1969.879518,972.8764492,1426.223077,1341.383441,1073.87574,1467.486818,1099.141935,921.0838509,982.8055215,1369.249615,722.3727273]

    # ax.scatter(W_S_ref, W_P_ref, label='Reference aircraft datapoints', color='black')
    ax.set_xlabel('W/S')
    ax.set_ylabel('W/P')
    # ax.legend(True)

    ax.set_title('Matching Diagram Sensitivity Mesh')

    ax.grid(True)

    plt.tight_layout()

    plt.savefig(output_filepath, dpi=200)

    plt.show()

def Weight_est_and_match_concept(ac : Aircraft,  # Change units
        type_to_use : str = "Single Engine Propeller Driven",
        W_S_plot: np.ndarray = np.arange(0,10000,5),
        W_P_or_T_W_plot: np.ndarray = np.arange(0,10000,5), 
        output_filepath_base: str = 'outputs/Matching_Diagram', 
        CL_max_step: float = 0.1, 
        A_step: float = 0.5,
        n_steps: int = 6) -> dict:

    # Initial matching diagram
    output_filepath = f"{output_filepath_base}.png"
    data = plot_matching_and_select_design_point(ac, type_to_use, W_S_plot, W_P_or_T_W_plot, output_filepath, show_plot=False)
    W_P = data['W/P']
    W_S = data['W/S']

    # Start storing history of W/P and W/S vs changes
    Initial_CL_LD = ac.requirements.landing['as_CL_max_la']
    Initial_A = ac.wing.aspect_ratio
    limiting_wp_constraint = data['limiting_wp_constraint']
    limiting_ws_constraint = data['limiting_ws_constraint']
    # W_P_history = [W_P]
    # W_S_history = [W_S]
    # CL_max_LD_history = [Initial_CL_LD]
    # A_history = [Initial_A]
    # limiting_wp_constraint_history = [limiting_wp_constraint]
    # limiting_ws_constraint_history = [limiting_ws_constraint]

    results_CL = sensitivity_study(
        ac, type_to_use,
        W_S_plot, W_P_or_T_W_plot, output_filepath_base,
        param='CL_max_LD',
        step=CL_max_step,
        n_steps=n_steps,
        initial_limiting_wp_constraint=limiting_wp_constraint,
        initial_limiting_ws_constraint=limiting_ws_constraint,
        initial_W_P=W_P,
        initial_W_S=W_S,
    )

    results_A = sensitivity_study(
        ac, type_to_use,
        W_S_plot, W_P_or_T_W_plot, output_filepath_base,
        param='A',
        step=A_step,
        n_steps=n_steps,
        initial_limiting_wp_constraint=limiting_wp_constraint,
        initial_limiting_ws_constraint=limiting_ws_constraint,
        initial_W_P=W_P,
        initial_W_S=W_S,
    )

    # Climb grad stuff
    RC, CG = max_RC_and_Climb_grad(W_S,
                                    W_P,
                                    k=k(ac)[0],
                                    rho=Atmosphere(ac.requirements.take_off['to_altitude'], ac.requirements.take_off['to_temp_shift']).density,
                                    eta_p=ac.engine.eta_prop,
                                    CD0=cd0(ac, type_to_use))

    # WEIGHT EST
    m_to: float = ac.weights.m_takeoff
    m_pl: float = ac.weights.m_payload
    # m_f_frac: float = sum(c1_m.energy_frac_needed(ac))  # Tuple with fuel_frac, battery_frac
    result = c1_m.energy_frac_needed(ac)
    if isinstance(result, tuple) and len(result) == 2:
        fuel_frac, bat_frac = result
        energy_frac = fuel_frac + bat_frac
    else:
        energy_frac = result
        fuel_frac = None
        bat_frac = result
    m_oe_frac = c1_m.operating_empty_frac(ac)
    # pl_mtow = m_pl/m_to
    mass_row = {'Fuel frac': fuel_frac, 'Battery frac': bat_frac, 'Energy frac': energy_frac, 'OEW/MTOW': m_oe_frac, 'PL/MTOW': m_pl/m_to, 'MTOW': m_to, 'Energy source mass': energy_frac*m_to, 'OEW': m_oe_frac*m_to, 'PL mass': m_pl, 'Sum mass fracs': m_oe_frac+energy_frac+m_pl/m_to, 'Max RoC [m/s]': RC, 'Max climb angle [deg at TO]': CG, 'S [m2]': m_to/W_S, 'P [kW]': m_to/W_P/1000}

    

    # # Check if need to check CL_max_landing change
    # CL_LD_resize: bool = Resize_CL_max_LD(limiting_ws_constraint)
    # if CL_LD_resize:

    #     try:
    #         # Change CL and calculate values to be tracked
    #         ac.requirements.landing['as_CL_max'] = Initial_CL_LD + CL_max_step
    #         assert (ac.requirements.landing['as_CL_max'] - (Initial_CL_LD + CL_max_step))<1e-3
            
    #         output_filepath_r1 = f"{output_filepath_base}_r1.png"
    #         data_r1 = plot_matching_and_select_design_point(ac, type_to_use, friction_source, s_wet_source,
    #                                                     W_S_plot, W_P_or_T_W_plot, output_filepath_r1)
            
    #         W_P_r1 = data_r1['W/P']
    #         W_S_r1 = data_r1['W/S']
    #         limiting_wp_constraint_r1 = data_r1['limiting_wp_constraint']
    #         limiting_ws_constraint_r1 = data_r1['limiting_ws_constraint']

    #         # Update lists
    #         W_P_history.append(W_P_r1)
    #         W_S_history.append(W_S_r1)
    #         CL_max_LD_history.append(Initial_CL_LD + CL_max_step)
    #         A_history.append(Initial_A)
    #         limiting_wp_constraint_history.append(limiting_wp_constraint_r1)
    #         limiting_ws_constraint_history.append(limiting_ws_constraint_r1)

    #         # Check second resize
    #         CL_LD_resize2: bool = Resize_CL_max_LD(limiting_ws_constraint)

    #         if CL_LD_resize2:
    #             try: 
    #                 ac.requirements.landing['as_CL_max'] = Initial_CL_LD + 2 * CL_max_step

    #                 output_filepath_r2 = f"{output_filepath_base}_r2.png"
    #                 data_r2 = plot_matching_and_select_design_point(ac, type_to_use, friction_source, s_wet_source,
    #                                                     W_S_plot, W_P_or_T_W_plot, output_filepath_r2)
                    
    #                 W_P_r2 = data_r2['W/P']
    #                 W_S_r2 = data_r2['W/S']
    #                 limiting_wp_constraint_r2 = data_r2['limiting_wp_constraint']
    #                 limiting_ws_constraint_r2 = data_r2['limiting_ws_constraint']

    #                 # Update lists
    #                 W_P_history.append(W_P_r2)
    #                 W_S_history.append(W_S_r2)
    #                 CL_max_LD_history.append(Initial_CL_LD + 2 * CL_max_step)
    #                 A_history.append(Initial_A)
    #                 limiting_wp_constraint_history.append(limiting_wp_constraint_r2)
    #                 limiting_ws_constraint_history.append(limiting_ws_constraint_r2)

    #             finally:
    #                 ac.requirements.landing['as_CL_max'] = Initial_CL_LD

    #     finally:
    #         ac.requirements.landing['as_CL_max'] = Initial_CL_LD
    
    # A_resize: bool = Resize_A(limiting_wp_constraint)
    # if A_resize:
    #     try:
    #         # Change A and calculate values to be tracked
    #         ac.wing.aspect_ratio = Initial_A + A_step
            
    #         output_filepath_r3 = f"{output_filepath_base}_r3.png"
    #         data_r3 = plot_matching_and_select_design_point(ac, type_to_use, friction_source, s_wet_source,
    #                                                     W_S_plot, W_P_or_T_W_plot, output_filepath_r3)
            
    #         W_P_r3 = data_r3['W/P']
    #         W_S_r3 = data_r3['W/S']
    #         limiting_wp_constraint_r3 = data_r3['limiting_wp_constraint']
    #         limiting_ws_constraint_r3 = data_r3['limiting_ws_constraint']

    #         # Update lists
    #         W_P_history.append(W_P_r3)
    #         W_S_history.append(W_S_r3)
    #         CL_max_LD_history.append(Initial_CL_LD)
    #         A_history.append(Initial_A + A_step)
    #         limiting_wp_constraint_history.append(limiting_wp_constraint_r3)
    #         limiting_ws_constraint_history.append(limiting_ws_constraint_r3)

    #         # Check second resize
    #         A_resize2: bool = Resize_A(limiting_wp_constraint_r3)

    #         if A_resize2:
    #             try: 
    #                 ac.wing.aspect_ratio = Initial_A + 2 * A_step

    #                 output_filepath_r4 = f"{output_filepath_base}_r4.png"
    #                 data_r4 = plot_matching_and_select_design_point(ac, type_to_use, friction_source, s_wet_source,
    #                                                     W_S_plot, W_P_or_T_W_plot, output_filepath_r4)
                    
    #                 W_P_r4 = data_r4['W/P']
    #                 W_S_r4 = data_r4['W/S']
    #                 limiting_wp_constraint_r4 = data_r4['limiting_wp_constraint']
    #                 limiting_ws_constraint_r4 = data_r4['limiting_ws_constraint']

    #                 # Update lists
    #                 W_P_history.append(W_P_r4)
    #                 W_S_history.append(W_S_r4)
    #                 CL_max_LD_history.append(Initial_CL_LD)
    #                 A_history.append(Initial_A + 2 * A_step)
    #                 limiting_wp_constraint_history.append(limiting_wp_constraint_r4)
    #                 limiting_ws_constraint_history.append(limiting_ws_constraint_r4)

    #             finally:
    #                 ac.wing.aspect_ratio = Initial_A

    #     finally:
    #         ac.wing.aspect_ratio = Initial_A
    
    # # Store results and outputs 
    # output = {
    #     "W/P": W_P_history,
    #     "W/S": W_S_history,
    #     "CL_max_LD": CL_max_LD_history,
    #     "A_history": A_history,
    #     "limiting_ws_constraint": limiting_ws_constraint_history,
    #     "limiting_wp_constraint": limiting_wp_constraint_history,
    # }

    return results_CL, results_A, mass_row, W_S, W_P


def run_sensitivity_study_save_results(aircraft_obj: list[object],
                                       concept_IDs: list[str] = ['CP_1', 'CP_2', 'CP_3'],
                                       W_S_plot: np.ndarray = np.arange(1,1250),
                                       W_P_or_T_W_plot: np.ndarray = np.arange(0.00000001,0.15,0.0001),
                                       CL_max_step: float = 0.1,
                                       A_step: float = 0.5,
                                       n_steps: int = 6,
                                       ) -> None:

    # Filepaths:
    output_dir = Path("outputs")
    folder = output_dir / 'Matching_concepts'
    folder.mkdir(parents=True, exist_ok=True)

    output_csv_path = folder / 'All_concepts_og_params_results.csv'
    output_csv_path1 = folder / 'All_concepts_mass_results.csv'

    # Original points dataframe for concepts
    rows_main = []
    rows_mass = []
    file_paths_A = []
    file_paths_CL = []

    ''' Start looping over concepts'''
    for i, ac in enumerate(aircraft_obj):
        Concept_ID: str = concept_IDs[i]
        # ac = loader.load(file, Aircraft)
        type_to_use = ac.requirements.general['standard_type']
        img_filepath_base = f"outputs/Matching_concepts/Sensitivity_study_graphs/{Concept_ID}_MD"
        output_CL, output_A, mass_row, W_S, W_P = Weight_est_and_match_concept(ac, type_to_use, W_S_plot, W_P_or_T_W_plot, img_filepath_base, CL_max_step, A_step, n_steps)

        # Add og results to main df and save its own df
        rows_main.append({'Concept_ID': i+1, 'W/S': output_CL['W/S'][0], 'W/P': output_CL['W/P'][0]})
        df1 = pd.DataFrame(output_CL)
        df2 = pd.DataFrame(output_A)
        filepath1 = folder / f'{Concept_ID}_CL_results.csv'
        filepath2 = folder / f'{Concept_ID}_A_results.csv'
        df1.to_csv(filepath1, index=False)
        df2.to_csv(filepath2, index=False)
        file_paths_CL.append(filepath1)
        file_paths_A.append(filepath2)

        # Climb grad stuff
        RC, CG = max_RC_and_Climb_grad(output_CL['W/S'][0],
                                       output_CL['W/P'][0],
                                       k=k(ac)[0],
                                       rho=Atmosphere(ac.requirements.take_off['to_altitude'], ac.requirements.take_off['to_temp_shift']).density,
                                       eta_p=ac.engine.eta_prop,
                                       CD0=cd0(ac, type_to_use))

        # WEIGHT EST
        m_to: float = ac.weights.m_takeoff
        m_pl: float = ac.weights.m_payload
        m_f_frac: float = sum(c1_m.energy_frac_needed(ac))  # Tuple with fuel_frac, battery_frac
        m_oe_frac = c1_m.operating_empty_frac(ac)
        # pl_mtow = m_pl/m_to
        rows_mass.append({'Concept_ID': i+1, 'Fuel/energy frac': m_f_frac, 'OEW/MTOW': m_oe_frac, 'PL/MTOW': m_pl/m_to, 'MTOW': m_to, 'Fuel/energy source mass': m_f_frac*m_to, 'OEW': m_oe_frac*m_to, 'PL mass': m_pl, 'Sum mass fracs': m_oe_frac+m_f_frac+m_pl/m_to, 'Max RoC [m/s]': RC, 'Max climb angle [deg at TO]': CG})


    # Save main df to csv
    df = pd.DataFrame(rows_main)
    df.to_csv(output_csv_path, index=False)
    df = pd.DataFrame(rows_mass)
    df.to_csv(output_csv_path1, index=False)

    return file_paths_A, file_paths_CL

def run_sensitivity_study_save_results1(aircraft_files: list[str] = ['yamls/aircraft.yaml','yamls/aircraft.yaml','yamls/aircraft.yaml'],
                                       concept_IDs: list[str] = ['CP_1', 'CP_2', 'CP_3'],
                                       W_S_plot: np.ndarray = np.arange(1,1250),
                                       W_P_or_T_W_plot: np.ndarray = np.arange(0.00000001,0.15,0.0001),
                                       CL_max_step: float = 0.1,
                                       A_step: float = 0.5,
                                       n_steps: int = 6,
                                       ) -> None:

    # Filepaths:
    output_dir = Path("outputs")
    folder = output_dir / 'Matching_concepts'
    folder.mkdir(parents=True, exist_ok=True)

    output_csv_path = folder / 'All_concepts_og_params_results.csv'
    output_csv_path1 = folder / 'All_concepts_mass_results.csv'

    # Original points dataframe for concepts
    rows_main = []
    rows_mass = []
    file_paths_A = []
    file_paths_CL = []

    ''' Start looping over concepts'''
    for i, file in enumerate(aircraft_files):
        Concept_ID: str = concept_IDs[i]
        ac = loader.load(file, Aircraft)
        type_to_use = ac.requirements.general['standard_type']
        img_filepath_base = f"outputs/Matching_concepts/Sensitivity_study_graphs/{Concept_ID}_MD"
        output_CL, output_A, mass_row, W_S, W_P = Weight_est_and_match_concept(ac, type_to_use, W_S_plot, W_P_or_T_W_plot, img_filepath_base, CL_max_step, A_step, n_steps)

        # Add og results to main df and save its own df
        rows_main.append({'Concept_ID': i+1, 'W/S': output_CL['W/S'][0], 'W/P': output_CL['W/P'][0]})
        df1 = pd.DataFrame(output_CL)
        df2 = pd.DataFrame(output_A)
        filepath1 = folder / f'{Concept_ID}_CL_results.csv'
        filepath2 = folder / f'{Concept_ID}_A_results.csv'
        df1.to_csv(filepath1, index=False)
        df2.to_csv(filepath2, index=False)
        file_paths_CL.append(filepath1)
        file_paths_A.append(filepath2)

        # Climb grad stuff
        RC, CG = max_RC_and_Climb_grad(output_CL['W/S'][0],
                                       output_CL['W/P'][0],
                                       k=k(ac)[0],
                                       rho=Atmosphere(ac.requirements.take_off['to_altitude'], ac.requirements.take_off['to_temp_shift']).density,
                                       eta_p=ac.engine.eta_prop,
                                       CD0=cd0(ac, type_to_use))

        # WEIGHT EST
        m_to: float = ac.weights.m_takeoff
        m_pl: float = ac.weights.m_payload
        m_f_frac: float = sum(c1_m.energy_frac_needed(ac))  # Tuple with fuel_frac, battery_frac
        m_oe_frac = c1_m.operating_empty_frac(ac)
        # pl_mtow = m_pl/m_to
        rows_mass.append({'Concept_ID': i+1, 'Fuel/energy frac': m_f_frac, 'OEW/MTOW': m_oe_frac, 'PL/MTOW': m_pl/m_to, 'MTOW': m_to, 'Fuel/energy source mass': m_f_frac*m_to, 'OEW': m_oe_frac*m_to, 'PL mass': m_pl, 'Sum mass fracs': m_oe_frac+m_f_frac+m_pl/m_to, 'Max RoC [m/s]': RC, 'Max climb angle [deg at TO]': CG})


    # Save main df to csv
    df = pd.DataFrame(rows_main)
    df.to_csv(output_csv_path, index=False)
    df = pd.DataFrame(rows_mass)
    df.to_csv(output_csv_path1, index=False)

    return file_paths_A, file_paths_CL

def plot_sensitivity_study1(
        A_csv_paths: list[str],
        CL_csv_paths: list[str],
        output_filepath: str,
        concept_names: list[str] = None) -> None:
    """
    Plots sensitivity study results for aspect ratio (A) and max lift coefficient (CL_max_LD)
    across multiple aircraft concepts.

    Args:
        A_csv_paths:      One CSV per concept for the A sensitivity study.
        CL_csv_paths:     One CSV per concept for the CL_max_LD sensitivity study.
        output_filepath:  Path to save the output plot.
        concept_names:    Optional list of concept names for the legend.
                          Defaults to 'Concept 0', 'Concept 1', etc.
    """
    assert len(A_csv_paths) == len(CL_csv_paths), \
        "Must provide the same number of CSV paths for A and CL studies."

    n_concepts = len(A_csv_paths)
    if concept_names is None:
        concept_names = [f"Concept {i+1}" for i in range(n_concepts)]

    colors = plt.cm.tab10(np.linspace(0, 1, n_concepts))
    marker_A  = 'x'
    marker_CL = 'o'

    # Collect data to scale axis limits (zoom in)
    all_W_S = []
    all_W_P = []
    datasets = []  # list of (A_df, CL_df) per concept

    for i, (a_path, cl_path) in enumerate(zip(A_csv_paths, CL_csv_paths)):
        df_A  = pd.read_csv(a_path)
        df_CL = pd.read_csv(cl_path)
        datasets.append((df_A, df_CL))
        all_W_S.extend(df_A['W/S'].tolist())
        all_W_S.extend(df_CL['W/S'].tolist())
        all_W_P.extend(df_A['W/P'].tolist())
        all_W_P.extend(df_CL['W/P'].tolist())

    # Scale axis limits
    pad = 0.05
    W_S_min, W_S_max = min(all_W_S), max(all_W_S)
    W_P_min, W_P_max = min(all_W_P), max(all_W_P)
    W_S_range = W_S_max - W_S_min or 1
    W_P_range = W_P_max - W_P_min or 1

    # Plot
    fig, ax = plt.subplots(figsize=(9, 6))

    for i, (df_A, df_CL) in enumerate(datasets):
        color = colors[i]
        label = concept_names[i]

        # Aspect ratio plotting
        ax.plot(
            df_A['W/S'], df_A['W/P'],
            marker=marker_A, linestyle='-', color=color,
            markersize=7, linewidth=1.5,
        )
        # Annotate each point with its A value
        for _, row in df_A.iterrows():
            ax.annotate(
                f"A={row['A']:.1f}",
                xy=(row['W/S'], row['W/P']),
                xytext=(4, 4), textcoords='offset points',
                fontsize=7, color=color,
            )

        # CL plotting
        ax.plot(
            df_CL['W/S'], df_CL['W/P'],
            marker=marker_CL, linestyle='--', color=color,
            markersize=7, linewidth=1.5,
        )
        # Annotate each point with its CL_max_LD value
        for _, row in df_CL.iterrows():
            ax.annotate(
                f"CL={row['CL_max_LD']:.2f}",
                xy=(row['W/S'], row['W/P']),
                xytext=(4, -10), textcoords='offset points',
                fontsize=7, color=color,
            )

    # Tight axis limits
    ax.set_xlim(W_S_min - pad * W_S_range, W_S_max + pad * W_S_range)
    ax.set_ylim(W_P_min - pad * W_P_range, W_P_max + pad * W_P_range)

    ax.set_xlabel('W/S', fontsize=12)
    ax.set_ylabel('W/P', fontsize=12)
    ax.set_title('Sensitivity Study: Design Point Shift', fontsize=13)
    ax.grid(True, linestyle='--', alpha=0.5)

    # Legend
    # Color patches for concepts
    concept_handles = [
        mlines.Line2D([], [], color=colors[i], marker='o', linestyle='-',
                      markersize=7, label=concept_names[i])
        for i in range(n_concepts)
    ]
    # Marker style handles for study type
    handle_A  = mlines.Line2D([], [], color='grey', marker=marker_A,
                               linestyle='-',  markersize=7, label='A sensitivity')
    handle_CL = mlines.Line2D([], [], color='grey', marker=marker_CL,
                               linestyle='--', markersize=7, label='$C_{L,max}$ sensitivity')

    ax.legend(
        handles=concept_handles + [handle_A, handle_CL],
        fontsize=9, loc='best', framealpha=0.8,
    )

    plt.tight_layout()
    plt.savefig(output_filepath, dpi=150)
    plt.show()

def plot_sensitivity_study(
        A_csv_paths: list[str],
        CL_csv_paths: list[str],
        output_filepath: str,
        concept_names: list[str] = None,
        param: str = 'both'  # options: 'A', 'CL_max_LD', 'both'
) -> None:
    """
    Plots sensitivity study results for aspect ratio (A) and max lift coefficient (CL_max_LD)
    across multiple aircraft concepts.

    Args:
        A_csv_paths:      One CSV per concept for the A sensitivity study.
        CL_csv_paths:     One CSV per concept for the CL_max_LD sensitivity study.
        output_filepath:  Path to save the output plot.
        concept_names:    Optional list of concept names for the legend.
                          Defaults to 'Concept 1', 'Concept 2', etc.
        param:            Which sensitivity to plot: 'A', 'CL_max_LD', or 'both' (default).
    """
    assert len(A_csv_paths) == len(CL_csv_paths), \
        "Must provide the same number of CSV paths for A and CL studies."
    assert param in ('A', 'CL_max_LD', 'both'), \
        f"param must be 'A', 'CL_max_LD', or 'both', got '{param}'"

    n_concepts = len(A_csv_paths)
    if concept_names is None:
        concept_names = [f"Concept {i+1}" for i in range(n_concepts)]

    colors = plt.cm.tab10(np.linspace(0, 1, n_concepts))
    marker_A  = 'x'
    marker_CL = 'o'

    # Collect data to scale axis limits (zoom in)
    all_W_S = []
    all_W_P = []
    datasets = []  # list of (A_df, CL_df) per concept

    for i, (a_path, cl_path) in enumerate(zip(A_csv_paths, CL_csv_paths)):
        df_A  = pd.read_csv(a_path)
        df_CL = pd.read_csv(cl_path)
        datasets.append((df_A, df_CL))
        if param in ('A', 'both'):
            all_W_S.extend(df_A['W/S'].tolist())
            all_W_P.extend(df_A['W/P'].tolist())
        if param in ('CL_max_LD', 'both'):
            all_W_S.extend(df_CL['W/S'].tolist())
            all_W_P.extend(df_CL['W/P'].tolist())

    # Scale axis limits
    pad = 0.05
    W_S_min, W_S_max = min(all_W_S), max(all_W_S)
    W_P_min, W_P_max = min(all_W_P), max(all_W_P)
    W_S_range = W_S_max - W_S_min or 1
    W_P_range = W_P_max - W_P_min or 1

    # Plot
    fig, ax = plt.subplots(figsize=(9, 6))

    for i, (df_A, df_CL) in enumerate(datasets):
        color = colors[i]

        if param in ('A', 'both'):
            ax.plot(
                df_A['W/S'], df_A['W/P'],
                marker=marker_A, linestyle='-', color=color,
                markersize=7, linewidth=1.5,
            )
            for _, row in df_A.iterrows():
                ax.annotate(
                    f"A={row['A']:.1f}",
                    xy=(row['W/S'], row['W/P']),
                    xytext=(4, 4), textcoords='offset points',
                    fontsize=7, color=color,
                )

        if param in ('CL_max_LD', 'both'):
            ax.plot(
                df_CL['W/S'], df_CL['W/P'],
                marker=marker_CL, linestyle='--', color=color,
                markersize=7, linewidth=1.5,
            )
            for _, row in df_CL.iterrows():
                ax.annotate(
                    f"CL={row['CL_max_LD']:.2f}",
                    xy=(row['W/S'], row['W/P']),
                    xytext=(4, -10), textcoords='offset points',
                    fontsize=7, color=color,
                )

    # Tight axis limits
    ax.set_xlim(W_S_min - pad * W_S_range, W_S_max + pad * W_S_range)
    ax.set_ylim(W_P_min - pad * W_P_range, W_P_max + pad * W_P_range)

    ax.set_xlabel('W/S', fontsize=12)
    ax.set_ylabel('W/P', fontsize=12)
    title_map = {'A': 'Aspect Ratio', 'CL_max_LD': '$C_{L,max}$', 'both': 'A and $C_{L,max}$'}
    ax.set_title(f'Sensitivity Study: {title_map[param]} Design Point Shift', fontsize=13)
    ax.grid(True, linestyle='--', alpha=0.5)

    # Legend
    concept_handles = [
        mlines.Line2D([], [], color=colors[i], marker='o', linestyle='-',
                      markersize=7, label=concept_names[i])
        for i in range(n_concepts)
    ]
    style_handles = []
    if param in ('A', 'both'):
        style_handles.append(mlines.Line2D([], [], color='grey', marker=marker_A,
                                           linestyle='-', markersize=7, label='A sensitivity'))
    if param in ('CL_max_LD', 'both'):
        style_handles.append(mlines.Line2D([], [], color='grey', marker=marker_CL,
                                           linestyle='--', markersize=7, label='$C_{L,max}$ sensitivity'))

    ax.legend(
        handles=concept_handles + style_handles,
        fontsize=9, loc='best', framealpha=0.8,
    )

    plt.tight_layout()
    plt.savefig(output_filepath, dpi=150)
    plt.show()

def run_concepts_c1(ac_lst: list[object],
                     output_filepath: str)->None:
    rows_lst = []
    ids_lst = []
    x_plot = []
    y_plot = []
    for i, ac in enumerate(ac_lst):
        Concept_ID = f'Final_CP_{i}'
        # type_to_use = ac.requirements.general['standard_type']
        type_to_use = ac.requirements.general['type_to_use']
        img_filepath_base = f"outputs/Matching_concepts/Sensitivity_study_graphs/{Concept_ID}_MD"
        output_CL, output_A, mass_row, W_S, W_P = Weight_est_and_match_concept(ac, type_to_use, W_S_plot=np.arange(1,1500), W_P_or_T_W_plot=np.arange(0.00000001,0.15,0.0001), output_filepath_base=img_filepath_base,CL_max_step=0.1, A_step=1, n_steps=1)
        rows_lst.append(mass_row)
        ids_lst.append(ac.name)
        x_plot.append(W_S)
        y_plot.append(W_P)
        print(f'ac name: {ac.name}')

    # Build the DataFrame
    df = pd.DataFrame(rows_lst)
    df.insert(0, 'ID', ids_lst)  # insert ID as first column

    df.to_csv(output_filepath, index=False)
    print(f' \n Done! \n Saved results to {output_filepath}')

    sample_ac = ac_lst[0]
    plot_matching_and_select_design_point(sample_ac, type_to_use=sample_ac.requirements.general['type_to_use'],  W_S_plot=np.arange(750,850), W_P_plot=np.arange(0.08,0.1,0.0001), output_filepath=f'concepts/Matching_diagram_summary.png', show_plot=True, requirement_to_meet='all', other_design_points_x=x_plot, other_design_points_y=y_plot, other_design_points_labels=ids_lst)


if __name__ == '__main__':
    file_path = "yamls/aircraft.yaml"
    target_class = Aircraft
    ac = loader.load(file_path, target_class)

    ac1 = Aircraft('Boosted_piston_taildragger',
                  loader.load('concepts/reqs_nturb.yaml', Requirements),
                  loader.load('yamls/mission.yaml', Mission),
                  loader.load('yamls/weights.yaml', Weights),
                  loader.load('concepts/wing_courier.yaml', Wing),
                  loader.load('yamls/fuselage.yaml', Fuselage),
                  loader.load('concepts/engine_piston_b.yaml', Engine))
    ac2 = Aircraft('Piston_hybrid_taildragger',
                  loader.load('concepts/reqs_nturb.yaml', Requirements),
                  loader.load('yamls/mission.yaml', Mission),
                  loader.load('yamls/weights.yaml', Weights),
                  loader.load('concepts/wing_courier.yaml', Wing),
                  loader.load('yamls/fuselage.yaml', Fuselage),
                  loader.load('concepts/engine_piston_e.yaml', Engine))
    ac3 = Aircraft('Boosted_turboprop_taildragger',
                  loader.load('concepts/reqs_turb.yaml', Requirements),
                  loader.load('yamls/mission.yaml', Mission),
                  loader.load('yamls/weights.yaml', Weights),
                  loader.load('concepts/wing_courier.yaml', Wing),
                  loader.load('yamls/fuselage.yaml', Fuselage),
                  loader.load('concepts/engine_tprop_b.yaml', Engine))
    ac4 = Aircraft('Turbine_hybrid_taildragger',
                  loader.load('concepts/reqs_turb.yaml', Requirements),
                  loader.load('yamls/mission.yaml', Mission),
                  loader.load('yamls/weights.yaml', Weights),
                  loader.load('concepts/wing_courier.yaml', Wing),
                  loader.load('yamls/fuselage.yaml', Fuselage),
                  loader.load('concepts/engine_turb_e.yaml', Engine))
    ac5 = Aircraft('H2_taildragger',
                  loader.load('concepts/reqs_nturb.yaml', Requirements),
                  loader.load('yamls/mission.yaml', Mission),
                  loader.load('yamls/weights.yaml', Weights),
                  loader.load('concepts/wing_courier.yaml', Wing),
                  loader.load('yamls/fuselage.yaml', Fuselage),
                  loader.load('concepts/engine_h2.yaml', Engine))
    ac6 = Aircraft('Boosted_piston_tricycle',
                  loader.load('concepts/reqs_nturb.yaml', Requirements),
                  loader.load('yamls/mission.yaml', Mission),
                  loader.load('yamls/weights.yaml', Weights),
                  loader.load('concepts/wing_electra.yaml', Wing),
                  loader.load('yamls/fuselage.yaml', Fuselage),
                  loader.load('concepts/engine_piston_b.yaml', Engine))
    ac7 = Aircraft('Piston_hybrid_tricycle',
                  loader.load('concepts/reqs_nturb.yaml', Requirements),
                  loader.load('yamls/mission.yaml', Mission),
                  loader.load('yamls/weights.yaml', Weights),
                  loader.load('concepts/wing_electra.yaml', Wing),
                  loader.load('yamls/fuselage.yaml', Fuselage),
                  loader.load('concepts/engine_piston_e.yaml', Engine))
    ac8 = Aircraft('Boosted_turboprop_tricycle',
                  loader.load('concepts/reqs_turb.yaml', Requirements),
                  loader.load('yamls/mission.yaml', Mission),
                  loader.load('yamls/weights.yaml', Weights),
                  loader.load('concepts/wing_electra.yaml', Wing),
                  loader.load('yamls/fuselage.yaml', Fuselage),
                  loader.load('concepts/engine_tprop_b.yaml', Engine))
    ac9 = Aircraft('Turbine_hybrid_tricycle',
                  loader.load('concepts/reqs_turb.yaml', Requirements),
                  loader.load('yamls/mission.yaml', Mission),
                  loader.load('yamls/weights.yaml', Weights),
                  loader.load('concepts/wing_electra.yaml', Wing),
                  loader.load('yamls/fuselage.yaml', Fuselage),
                  loader.load('concepts/engine_turb_e.yaml', Engine))
    ac10 = Aircraft('H2_tricycle',
                  loader.load('concepts/reqs_nturb.yaml', Requirements),
                  loader.load('yamls/mission.yaml', Mission),
                  loader.load('yamls/weights.yaml', Weights),
                  loader.load('concepts/wing_electra.yaml', Wing),
                  loader.load('yamls/fuselage.yaml', Fuselage),
                  loader.load('concepts/engine_h2.yaml', Engine))
    

    # ac_lst = [ac1, ac2, ac3, ac4, ac5, ac6, ac7, ac8, ac9, ac10]
    ac_lst = [ac1, ac2, ac3, ac4, ac6, ac7, ac8, ac9]

    run_concepts_c1(ac_lst, output_filepath=f'concepts/test_c1_results.csv')
    # # Sensitivity study matching plot output path:
    # output_dir = Path("outputs")
    # folder = output_dir / 'Matching_concepts'
    # folder.mkdir(parents=True, exist_ok=True)
    # output_path1 = folder / "Sensitivity_study_graph_A.png"
    # output_path2 = folder / "Sensitivity_study_graph_CL.png"
    # output_path3 = folder / "Mesh_sensitivity_study_graph.png"

    # folder1 = folder / 'Sensitivity_study_graphs'
    # folder1.mkdir(parents=True, exist_ok=True)

    # A_values = np.arange(6, 9, 1)
    # CL_values = np.arange(1.8, 3.3, 0.1)

    # df_mesh = sensitivity_mesh(
    #     ac,
    #     type_to_use=ac.requirements.general['standard_type'],
    #     W_S_plot=np.arange(1,1500),
    #     W_P_or_T_W_plot=np.arange(0.00000001,0.15,0.0001),
    #     output_filepath_base="outputs/test",
    #     A_values=A_values,
    #     CL_values=CL_values
    # )

    # df_mesh.to_csv("outputs/sensitivity_mesh.csv", index=False)
    # plot_sensitivity_mesh(df_mesh, output_path3)

    # file_paths_A, file_paths_CL = run_sensitivity_study_save_results(aircraft_obj=[ac], concept_IDs=['CP_test'], n_steps=2)
    # plot_sensitivity_study(file_paths_A, file_paths_CL, output_path1, param='A')
    # plot_sensitivity_study(file_paths_A, file_paths_CL, output_path2, param='CL_max_LD')
    # plot_matching_and_select_design_point(ac,W_P_plot=np.arange(0.00000001,0.15,0.0001), W_S_plot=np.arange(1,1250), show_plot=True)
