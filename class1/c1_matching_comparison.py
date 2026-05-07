from classes.aircraft_2 import loader, Aircraft
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
    if limit_W_P == "Cruise speed" or "AEO RoC" or "AEO Climb gradient" or "AEO Climb gradient (turbine)" or "Balked landing" or "Balked landing (turbine)" or "OEI RoC/Climb gradient I (turbine)" or "OEI RoC/Climb gradient II (turbine)":
        resize = True
    return resize
    # elif limit_W_P == :  # NOTE: fill in names of lines and check it runs
    #     resize = False
    #     return resize

def Resize_CL_max_LD(limit_W_S):
    resize = False
    if limit_W_S == "Stall speed" or "Landing field length":
        resize = True
    return resize

def sensitivity_study(
        ac: Aircraft,
        type_to_use: str,
        friction_source: str,
        s_wet_source: str,
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
            data_i = plot_matching_and_select_design_point(
                ac, type_to_use, friction_source, s_wet_source,
                W_S_plot, W_P_or_T_W_plot, output_filepath_i
            )

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

def Weight_est_and_match_concept(ac : Aircraft,  # Change units
        type_to_use : str = "Single Engine Propeller Driven",
        friction_source : str = 'lookups/skin_fric.csv',
        s_wet_source : str = 'lookups/s_wets.csv',
        W_S_plot: np.ndarray = np.arange(0,10000,5),
        W_P_or_T_W_plot: np.ndarray = np.arange(0,10000,5), 
        output_filepath_base: str = 'outputs/Matching_Diagram.png', 
        CL_max_step: float = 0.2, 
        A_step: float = 2,
        n_steps: int = 6) -> dict:

    # Initial matching diagram
    output_filepath = f"{output_filepath_base}.png"
    data = plot_matching_and_select_design_point(ac, type_to_use, friction_source, s_wet_source,
                                                 W_S_plot, W_P_or_T_W_plot, output_filepath)
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
        ac, type_to_use, friction_source, s_wet_source,
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
        ac, type_to_use, friction_source, s_wet_source,
        W_S_plot, W_P_or_T_W_plot, output_filepath_base,
        param='A',
        step=A_step,
        n_steps=n_steps,
        initial_limiting_wp_constraint=limiting_wp_constraint,
        initial_limiting_ws_constraint=limiting_ws_constraint,
        initial_W_P=W_P,
        initial_W_S=W_S,
    )

    

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

    return results_CL, results_A


def run_sensitivity_study_save_results(aircraft_files: list[str] = ['yamls/aircraft.yaml','yamls/aircraft.yaml','yamls/aircraft.yaml'],
                                       concept_IDs: list[str] = ['CP_1', 'CP_2', 'CP_3'],
                                       W_S_plot: np.ndarray = np.arange(0,1200,1),
                                       W_P_or_T_W_plot: np.ndarray = np.arange(0,100000,2),
                                       friction_source: str = 'lookups/skin_fric.csv',
                                       s_wet_source: str = 'lookups/s_wets.csv',
                                       CL_max_step: float = 0.2,
                                       A_step: float = 2,
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

    ''' Start looping over concepts'''
    for i, file in enumerate(aircraft_files):
        Concept_ID: str = concept_IDs[i]
        ac = loader.load(file, Aircraft)
        type_to_use = ac.requirements.general['standard_type']
        img_filepath_base = f"outputs/{Concept_ID}_MD"
        output_CL, output_A = Weight_est_and_match_concept(ac, type_to_use, friction_source, s_wet_source, W_S_plot, W_P_or_T_W_plot, img_filepath_base, CL_max_step, A_step, n_steps)

        # Add og results to main df and save its own df
        rows_main.append({'Concept_ID': i+1, 'W/S': output_CL['W/S'][0], 'W/P': output_CL['W/P'][0]})
        df1 = pd.DataFrame(output_CL)
        df2 = pd.DataFrame(output_A)
        filepath1 = folder / f'{Concept_ID}_CL_results.csv'
        filepath2 = folder / f'{Concept_ID}_A_results.csv'
        df1.to_csv(filepath1, index=False)
        df2.to_csv(filepath2, index=False)

        # INSERT WEIGHT EST
        m_to = ac.weights.m_takeoff
        m_pl = ac.weights.m_payload
        m_f_frac = c1_m.energy_frac_needed(ac)  # Tuple with fuel_frac, battery_frac
        m_oe_frac = c1_m.operating_empty_frac(ac)

        rows_mass.append({'Concept_ID': i+1, 'Fuel frac': m_f_frac, 'OEW/MTOW': m_oe_frac, 'PL/MTOW': m_pl/m_to, 'MTOW': m_to, 'Fuel mass': m_f_frac*m_to, 'OEW': m_oe_frac*m_to, 'PL mass': m_pl, 'Sum mass fracs': m_oe_frac+m_f_frac+m_pl/m_to})


    # Save main df to csv
    df = pd.DataFrame(rows_main)
    df.to_csv(output_csv_path, index=False)
    df = pd.DataFrame(rows_mass)
    df.to_csv(output_csv_path1, index=False)


def plot_sensitivity_study(
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