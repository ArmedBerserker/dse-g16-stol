"""
design_loop.py
==============
Iterative aircraft sizing / design loop.

Usage
-----
    from design_loop import run_design_loop, DesignLoopConfig
    from structures import loader, Aircraft

    ac = loader.load("path/to/aircraft.yaml", Aircraft)
    final = run_design_loop(ac)

The history file is written once per epoch and is STRICTLY a read-only log.
It is never read back at runtime; it exists only for post-processing / review.
To restart from a checkpoint see `load_epoch_as_dict()` at the bottom.
"""

from __future__ import annotations

import copy
import json
import os
import sys
from dataclasses import asdict, dataclass, field
from typing import Callable
import numpy as np

from classes.aircraft_2 import Aircraft, loader, Requirements, Mission, Fuselage, Wing, Engine, Weights, Empennage, HLD_and_AIL, Landing_Gear

from class1 import c1_m, prelim_drag, matching_diagram, c2_drag, c1_gear_sizing, c1_planform_sizing, c1_loading_and_empennage, c1_fuselage, c2_drag_new, c1_lift_and_ailerons
from c2_m import W_oe_and_cg_from_nose, W_to_new, loading_diagram, x_cg_structural_from_nose, overlay_wing_pos_and_scissor_plot
from lookups.consts import *

''' TO DO:
    - Position ht on vt right for t-tail
    - Add way to save and load aircraft so we can make the plots separately
    - Add checks in landing gear code to check sqrts are only taken of +ve values and the locations make sense
    - Add wing tip to C_D_L in drag'''

# 0. Helper function
def check_power_requirement(ac: Aircraft):
    e = ac.engine
    P_a_cr = e.eta_1 * e.eta_3 * e.engine_power_cruise
    P_r_cr = e.power_cr
    phi = (e.super_cap_power) / (e.super_cap_power + e.engine_power_cruise)
    P_a_to = e.eta_1 * e.engine_power_takeoff * e.eta_3 + e.eta_2 * e.super_cap_power * e.eta_3
    P_r_to = e.power_to
    condition_cr = P_r_cr < P_a_cr
    condition_to = P_r_to < P_a_to
    return condition_cr, condition_to

# 1. Configuration

@dataclass
class ConvergenceParam:
    """One parameter to watch during convergence checking."""
    dotpath : str          # e.g. "weights.m_takeoff" or "wing.area"
    tolerance : float      # absolute tolerance


@dataclass
class DesignLoopConfig:
    """
    All knobs for the design loop in one place.
    Modify this to change which parameters are watched and how tightly.
    """
    max_epochs   : int   = 100
    history_file : str   = "aircraft_history.json"

    # Parameters checked for convergence.
    # Paths use dot-notation into the nested dataclass tree.
    convergence_params: list[ConvergenceParam] = field(default_factory=lambda: [
        ConvergenceParam("weights.m_takeoff", 1.0),    # kg
        ConvergenceParam("wing.area",         0.01),   # m2
        ConvergenceParam("wing.ld",           0.001),
        # add / remove entries as your model grows
    ])



# 2. Individual functions for iteration steps 
#       
#       - Update params of ac object
#       - Order of execution defined in ITERATION_STEPS

def pre_loop_calculations(ac: Aircraft) -> Aircraft:
    # Calculate L/D for range eqn
    if ac.engine.count == 1:
        type = "Single Engine Propeller Driven"
    elif ac.engine.count > 1:
        type = "Twin Engine Propeller Driven"
    prelim_drag.prelim_drag(ac, type, update_ac=True)

    # Calculate Class I masses:
    c1_m.energy_frac_needed(ac, Phi=ac.engine.Phi, update_ac=True)
    c1_m.operating_empty_frac(ac, source_for_fracs='specific', engine_type=ac.engine.alpha_p_id, gear_type=ac.landing_gear.gear_type, update_ac=True)

    # print(f' Fuel weight: {ac.weights.m_fuel}')
    # Matching diagram:
    type_to_use = "Twin Engine Propeller Driven"
    if ac.engine.count == 1:
        type_to_use = "Single Engine Propeller Driven"
    data_to = matching_diagram.plot_matching_and_select_design_point(ac, type_to_use, W_P_plot=np.arange(0.00000001,0.15,0.0001), W_S_plot=np.arange(1,1250), output_filepath='outputs/Iteration_matching_plot.png', requirement_to_meet='to')
    W_P_to = data_to['W/P']
    W_S = data_to['W/S']
    ac.wing.area = ac.weights.m_takeoff * g / W_S
    data_cr = matching_diagram.plot_matching_and_select_design_point(ac, type_to_use, W_P_plot=np.arange(0.00000001,0.15,0.0001), W_S_plot=np.arange(1,1250), output_filepath='outputs/Iteration_matching_plot.png', requirement_to_meet='cruise')
    W_P_cr = data_cr['W/P']
    ac.engine.power_cr = ac.weights.m_takeoff * g / W_P_cr
    ac.engine.power_to = max(ac.weights.m_takeoff * g / W_P_to, ac.weights.m_takeoff * g / W_P_cr)

    print(f' W_e: {ac.weights.m_empty}')

    # Fuselage sizing:
    c1_fuselage.calculate_fuselage_parameters

def compute_aerodynamics(ac: Aircraft) -> Aircraft:
    """
    Step 1 — update aerodynamic coefficients and L/D from current geometry.
    Replace the body with your aero model.
    """
    # Wing planform
    c1_planform_sizing.size_wing_planform(ac)
    # HLD sizing
    c1_lift_and_ailerons.size_ailerons(ac, update_ac=True)
    c1_lift_and_ailerons.size_HLD(ac, update_ac=True)

    # Lift and drag params?

    return ac

def compute_landing_gear_positions(ac: Aircraft, epoch) -> Aircraft:
    ac = c1_gear_sizing.size_tires(ac, update_ac=True)
    ac = c1_gear_sizing.tire_location(ac, update_ac=True)
    Z1 = ac.landing_gear.height_mlg
    Z2 = ac.landing_gear.height_nlg
    ac.landing_gear.height_mlg = max(Z1, ac.landing_gear.selected_mlg_tire["Outside Diameter Max (In)"] * 2.54 / 100)
    ac.landing_gear.height_nlg = max(Z2, ac.landing_gear.selected_nlg_tire["Outside Diameter Max (In)"] * 2.54 / 100)
    return ac

def compute_empennage(ac: Aircraft, epoch) -> Aircraft:
    ac = c1_loading_and_empennage.size_empennage_planform(ac, epoch)
    return ac

def compute_class_I_mass(ac: Aircraft) -> Aircraft:
    """Initial mass estimate for all but first epoch"""
    oew_mtow = ac.weights.oew_frac
    result = c1_m.energy_frac_needed(ac)
    if ac.engine.engine_type == 'prop': # prop, bat or hyb
        fuel_frac = result[0]
        bat_frac = 0
        energy_frac = fuel_frac
    elif ac.engine.engine_type == 'bat':
        fuel_frac = 0
        bat_frac = result[0]
        energy_frac = bat_frac
    elif ac.engine.engine_type == 'hyb':
        fuel_frac, bat_frac = result
        energy_frac = fuel_frac + bat_frac
    # if isinstance(result, tuple) and len(result) == 2:
    #     fuel_frac, bat_frac = result
    #     energy_frac = fuel_frac + bat_frac
    # else:
    #     bat_frac = float(result[0])
    #     fuel_frac = 0.0
    #     energy_frac = bat_frac
    mtow = ac.weights.m_payload / (1 - fuel_frac - bat_frac - oew_mtow)

    ac.weights.m_takeoff = mtow 
    ac.weights.m_empty = oew_mtow * mtow
    ac.weights.m_energy = energy_frac * mtow
    ac.weights.m_fuel = fuel_frac * mtow
    ac.weights.m_battery = bat_frac * mtow
    ac.weights.pl_frac = ac.weights.m_payload / mtow

    type_to_use = "Twin Engine Propeller Driven"
    if ac.engine.count == 1:
        type_to_use = "Single Engine Propeller Driven"
    # NOTE: check if we need to add tw options for requirements to meet or add W/P result used by Shubhankar for weight est
    data_to = matching_diagram.plot_matching_and_select_design_point(ac, type_to_use, W_P_plot=np.arange(0.00000001,0.15,0.0001), W_S_plot=np.arange(1,1250), output_filepath='outputs/Iteration_matching_plot.png', requirement_to_meet='to', initial_est=False, show_plot=False)
    W_P_to = data_to['W/P']
    W_S = data_to['W/S']
    # NOTE: check if we need to update multiple power values and if they exist already
    ac.wing.area = ac.weights.m_takeoff * g / W_S
    data_cr = matching_diagram.plot_matching_and_select_design_point(ac, type_to_use, W_P_plot=np.arange(0.00000001,0.15,0.0001), W_S_plot=np.arange(1,1250), output_filepath='outputs/Iteration_matching_plot.png', requirement_to_meet='cruise', initial_est=False, show_plot=False)
    W_P_cr = data_cr['W/P']
    ac.engine.power_cr = ac.weights.m_takeoff * g / W_P_cr
    ac.engine.power_to = max(ac.weights.m_takeoff * g / W_P_to, ac.weights.m_takeoff * g / W_P_cr)
    return ac

def compute_class_II_mass_and_cg(ac: Aircraft, iteration: int) -> Aircraft:
    """Component build-up from actual geometry computed this epoch."""
    x_le_w = ac.wing.x_le
    m_tfo = 0.007 * ac.weights.m_takeoff
    # pie_chart_output_path = f'outputs/Class_II_weight/OEW_pie_chart_{iteration}.png'
    # struc_pie_chart_output_path = f'outputs/Class_II_weight/Structure_pie_chart_{iteration}.png'
    W_oe, x_cg_oe, ac = W_oe_and_cg_from_nose(ac, update_ac=True)
    W_to, ac = W_to_new(ac, m_res=0, update_ac=True)
    x_cg_struc, x_cg_data, ac = x_cg_structural_from_nose(ac, update_ac=True)
    fwd_cg, aft_cg, ac, xc_cg_nose_ftb, x_cg_nose_btf = loading_diagram(x_le_w, ac, update_ac_cgs=True)
    return ac

def tail_sizing_wing_positioning(ac: Aircraft, epoch: int) -> Aircraft:
    loading_diagram(ac.wing.x_le, ac, show_plot=False, output_filepath='outputs/init_loading_diagram.png', update_ac_cgs=False)
    wing_pos_arr = np.arange(0,1.01,0.01)
    # print(f'wing positions: {wing_pos_arr}')
    stability_output = overlay_wing_pos_and_scissor_plot(ac, x_le_w_fus_length_arr=wing_pos_arr, output_filepath=f'outputs/Stability_and_Control/Scissor_plot_{epoch}', show_plot=True, update_ac=False)
    if stability_output[0]>0:
        print(f' \n Sh_S: {stability_output[0]}, x_lemac/mac: {stability_output[1]}, aft_cg: {stability_output[2]}, fwd_cg: {stability_output[3]}, x_le: {stability_output[4]}')
        # NOTE: add updating Sh_S and wing_pos stored in stability output, add option to not update value if not happy with value
        update_ac = int(input('Enter 0 if you want to update these parameters into the aircraft, else 1'))
        if update_ac == 0:
            ac.weights.x_cg_aft = stability_output[2]
            ac.weights.x_cg_fwd = stability_output[3]
            ac.empennage.horizontal_tail['area div S'] = stability_output[0]
            ac.wing.x_le =  stability_output[4]
            ac.empennage.horizontal_tail['area'] = stability_output[0] * ac.wing.area
            return ac, epoch
        else: 
            epoch += 1e10
            return ac, epoch
    else: 
        return ac, epoch
    

def Class_II_drag(ac: Aircraft, epoch):
    # Cruise
    CD0 = c2_drag_new.CD0(ac, n_engine_inoperative=0, flight_condition='cruise', update_ac=True)
    c2_drag_new.C_D_L(ac, CD0, flight_condition='cruise', update_ac=True, wing_tip=False)
    # Take-off
    ac.wing.CD0_to = c2_drag_new.CD0(ac, n_engine_inoperative=0, flight_condition='take-off', update_ac=False)
    CDi, ac.wing.e_to, ac.wing.k_to, ac.wing.ld_to = c2_drag_new.C_D_L(ac, ac.wing.CD0_to, flight_condition='take-off', update_ac=True, wing_tip=False)
    # Landing
    ac.wing.CD0_ld = c2_drag_new.CD0(ac, n_engine_inoperative=0, flight_condition='landing', update_ac=False)
    CDi, ac.wing.e_ld, ac.wing.k_ld, ac.wing.ld_landing = c2_drag_new.C_D_L(ac, ac.wing.CD0_ld, flight_condition='landing', update_ac=True, wing_tip=False)
    return ac

ITERATION_STEPS1: list[Callable[[Aircraft], Aircraft]] = [
    compute_aerodynamics
]

ITERATION_STEPS2: list[Callable[[Aircraft], Aircraft]] = [
    compute_landing_gear_positions,
    compute_empennage,
    Class_II_drag
]


# 3. Define single iteration cycle
def run_iteration(ac: Aircraft,
                  epoch: int,
                  INNER_TOLERANCE: float = 0.01  # 1 % standard for Class I and II OEW convergence tolerance
                  ) -> tuple[Aircraft, bool]:
    """
    Returns the updated aircraft AND whether the inner mass loop
    (class I vs class II) converged this epoch.
    """

    # Class I mass
    if epoch > 1:
        ac = compute_class_I_mass(ac)
    mtow_class_I = ac.weights.m_takeoff      # snapshot before geometry steps
    if epoch == 1:
        tricycle_condition = (ac.landing_gear.gear_type == 'tricycle')
        c1_loading_and_empennage.class_I_loading_cgs(ac, tricycle_condition, update_ac=True)

    # Geometry / aero / propulsion etc.
    """Execute every registered step in order."""
    for step in ITERATION_STEPS1:
        ac = step(ac)

    if epoch > 1:
        # Initial loading
        c1_loading_and_empennage.classI_loading_and_cgs_2(ac, update_ac=True)

    for step in ITERATION_STEPS2:
        ac = step(ac, epoch)

    # Class II mass 
    ac = compute_class_II_mass_and_cg(ac, epoch)
    mtow_class_II = ac.weights.m_takeoff
    # print(ac)

    # ── Inner convergence check ───────────────────────────────────────────
    # inner_converged = (
    #     abs(mtow_class_II - mtow_class_I) / mtow_class_I
    # ) < INNER_TOLERANCE
    inner_converged = bool(
        abs(mtow_class_II - mtow_class_I) / mtow_class_I
        < INNER_TOLERANCE
    )

    # Stability and control and Class II drag
    ac, epoch = tail_sizing_wing_positioning(ac, epoch)

    return ac, inner_converged, epoch


# 4. Check convergence
def _get_nested(ac: Aircraft, dotpath: str) -> float:
    """
    Resolve a dot-separated attribute path against the aircraft object.

    Examples
    --------
    _get_nested(ac, "weights.m_takeoff")  ->  ac.weights.m_takeoff
    _get_nested(ac, "wing.area")          ->  ac.wing.area
    """
    obj = ac
    for attr in dotpath.split("."):
        obj = getattr(obj, attr)
    return obj

def has_converged(prev: Aircraft, curr: Aircraft,
                  params: list[ConvergenceParam]) -> tuple[bool, dict]:
    """
    Return (converged, deltas) where deltas maps each dotpath to its
    absolute change this epoch.  Useful for printing / debugging.
    """
    deltas: dict[str, float] = {}
    converged = True

    for cp in params:
        prev_val = _get_nested(prev, cp.dotpath)
        curr_val = _get_nested(curr, cp.dotpath)

        # Skip parameters that haven't been set yet (still None)
        if prev_val is None or curr_val is None:
            continue

        delta = abs(curr_val - prev_val)
        deltas[cp.dotpath] = delta

        if delta > cp.tolerance:
            converged = False

    return converged, deltas


# 5. History logging -> read-only, can look at after running, loop always uses memory not this log
def _snapshot(epoch: int, ac: Aircraft, deltas: dict, inner_converged: bool) -> dict:
    """Build a JSON-serialisable record for one epoch."""
    return {
        "epoch"  : epoch,
        "inner_converged" : inner_converged,   # class I vs II agreement
        "deltas" : deltas,          # convergence deltas for quick inspection
        "aircraft": asdict(ac),     # full recursive snapshot of every field
    }

def _append_epoch(record: dict, history_file: str) -> None:
    """
    Append one epoch record to the history file.

    The file holds a JSON array.  We load it, append, and rewrite.
    For very long runs (thousands of epochs) you may prefer a JSON-Lines
    format instead — see the note in load_epoch_as_dict() below.
    """
    if os.path.exists(history_file):
        with open(history_file, "r") as f:
            history: list = json.load(f)
    else:
        history = []

    history.append(record)

    with open(history_file, "w") as f:
        json.dump(history, f, indent=2)


# 6. Main loop for convergence
def run_design_loop(
    ac      : Aircraft,
    config  : DesignLoopConfig | None = None,
    INNER_TOLERANCE: float = 0.01,
) -> Aircraft:
    """
    Run the iterative design loop until convergence or max_epochs.

    Parameters
    ----------
    ac      : Aircraft loaded from YAML (will be mutated in-place each epoch).
    config  : DesignLoopConfig — uses defaults if not supplied.

    Returns
    -------
    The final Aircraft object after convergence (or hitting max_epochs).
    The full per-epoch history is written to config.history_file as a
    read-only log; it is never read back by this function.

    I added a check for if the tolerance entered is a fraction not 
    a percentage because I kept messing this up lol
    """

    # Add check for if tolerance was added as % instead of fraction
    if INNER_TOLERANCE > 1:
        user_input = input(f"INNER_TOLERANCE={INNER_TOLERANCE} looks like a percentage. Convert to fraction (INNER_TOLERANCE={INNER_TOLERANCE/100:.4f})? [y/n]: ")
        if user_input.lower() == 'y':
            INNER_TOLERANCE = INNER_TOLERANCE / 100
            print(f"Converted: INNER_TOLERANCE={INNER_TOLERANCE}")
        else:
            raise ValueError(f"Iteration loop stopped. Please provide INNER_TOLERANCE as a fraction (0-1), got INNER_TOLERANE={INNER_TOLERANCE}")
        
    if config is None:
        config = DesignLoopConfig()

    # Start fresh — remove any leftover history from a previous run
    if os.path.exists(config.history_file):
        os.remove(config.history_file)

    print(f"Starting design loop  (max {config.max_epochs} epochs)")
    print(f"History log → {os.path.abspath(config.history_file)}\n")
    _print_header(config.convergence_params)

    insufficient_to_power_counter = 0
    insufficient_cr_power_counter = 0
    pre_loop_calculations(ac)
    for epoch in range(1, config.max_epochs + 1):

        # Deep-copy BEFORE the iteration so we can diff afterwards
        prev = copy.deepcopy(ac)

        # Run one full design cycle
        ac, inner_converged, epoch_change = run_iteration(ac, epoch, INNER_TOLERANCE)
        if epoch_change > config.max_epochs + 1:
            break

        # Power requirement check
        enough_cr_power, enough_to_power =  check_power_requirement(ac)
        if enough_to_power:
            insufficient_to_power_counter = 0
        else:
            insufficient_to_power_counter += 1
        if enough_cr_power:
            insufficient_cr_power_counter = 0
        else:
            insufficient_cr_power_counter += 1

        # Check convergence
        converged, deltas = has_converged(prev, ac, config.convergence_params)
        print(type(inner_converged))

        # Log to file (read-only history — never read back by the loop)
        record = _snapshot(epoch, ac, deltas, inner_converged)
        _append_epoch(record, config.history_file)

        # Console output
        _print_epoch(epoch, ac, deltas, config.convergence_params, converged, inner_converged)

        # Stop if not enough power
        if insufficient_to_power_counter >5:
            print(f"\n⚠  Insufficient take-off power for 6 consecutive epochs.")
            break
        if insufficient_cr_power_counter >5:
            print(f"\n⚠  Insufficient cruise power for 6 consecutive epochs.")
            break

        if converged:
            print(f"\n✓  Converged after {epoch} epoch(s).")
            break

    else:
        print(f"\n ⚠  Reached max epochs ({config.max_epochs}) without full convergence.")

    return ac


# 7. Helper functions for printing stuff
def _print_header(params: list[ConvergenceParam]) -> None:
    col_w = 18
    header = f"{'Epoch':>6}  " + "".join(f"{p.dotpath:>{col_w}}" for p in params)
    print(header)
    print("─" * len(header))

def _print_epoch(epoch: int, ac: Aircraft, deltas: dict,
                 params: list[ConvergenceParam], converged: bool,
                 inner_converged: bool) -> None:
    col_w = 18
    row = f"{epoch:>6}  "
    for cp in params:
        val = _get_nested(ac, cp.dotpath)
        row += f"{'—':>{col_w}}" if val is None else f"{val:>{col_w}.4g}"
        # if val is None:
        #     row += f"{'—':>{col_w}}"
        # else:
        #     row += f"{val:>{col_w}.4g}"
    row += "  [I≈II ✓]" if inner_converged else "  [I≈II ✗]"
    row += "  ✓" if converged else ""
    print(row)



# 8. Stuff for after running
def load_history(history_file: str = "aircraft_history.json") -> list[dict]:
    """
    Load the full history log from a finished run.

    Returns a list of epoch records, each with keys:
        'epoch', 'deltas', 'aircraft'

    The 'aircraft' value is a plain dict — NOT an Aircraft object.
    This is intentional: history is raw data, not a live object.

    Example
    -------
        history = load_history()
        mtow_over_time = [r["aircraft"]["weights"]["m_takeoff"] for r in history]
    """
    with open(history_file, "r") as f:
        return json.load(f)


def load_epoch_as_dict(epoch: int,
                       history_file: str = "aircraft_history.json") -> dict:
    """
    Return the raw aircraft dict for one specific epoch.

    Use this to inspect a particular iteration, NOT to feed back into
    the design loop.  If you want to restart a run from a checkpoint,
    edit your YAML file with the converged values and call loader.load()
    again — that way Engine.__post_init__ runs correctly on the fresh values.

    Example
    -------
        state = load_epoch_as_dict(10)
        print(state["wing"]["area"])
    """
    history = load_history(history_file)
    for record in history:
        if record["epoch"] == epoch:
            return record["aircraft"]
    raise ValueError(f"Epoch {epoch} not found in {history_file}")


def summarise_convergence(history_file: str = "aircraft_history.json") -> None:
    """
    Print a concise convergence summary table from a finished run.

    Example
    -------
        summarise_convergence("aircraft_history.json")
    """
    history = load_history(history_file)
    print(f"{'Epoch':>6}  {'Parameter':<30}  {'Delta':>12}")
    print("─" * 54)
    for record in history:
        epoch = record["epoch"]
        for param, delta in record.get("deltas", {}).items():
            print(f"{epoch:>6}  {param:<30}  {delta:>12.6g}")
        print()


if __name__ == "__main__":
    ac = Aircraft('Boosted_turboprop_tricycle',
                loader.load('concepts/reqs_turb.yaml', Requirements),
                loader.load('yamls/mission.yaml', Mission),
                loader.load('yamls/weights.yaml', Weights),
                loader.load('concepts/wing_electra.yaml', Wing),
                loader.load('concepts/fus_tri.yaml', Fuselage),
                loader.load('concepts/engine_tprop_b.yaml', Engine),
                loader.load('concepts/tricycle_empennage.yaml', Empennage),
                loader.load('yamls/HLD_and_ailerons.yaml', HLD_and_AIL),
                loader.load('concepts/tricycle_gear.yaml', Landing_Gear))

    config = DesignLoopConfig(
        max_epochs   = 3,
        history_file = "aircraft_history.json",
        convergence_params = [
            ConvergenceParam("weights.m_empty", 1.0),
            ConvergenceParam("wing.area", 0.01),
            ConvergenceParam("wing.ld", 0.01),
        ],
    )

    final_aircraft = run_design_loop(ac, config)
    print("\nFinal aircraft state:")
    print(final_aircraft)