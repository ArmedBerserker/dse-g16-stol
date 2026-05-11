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
from dataclasses import asdict, dataclass, field
from typing import Callable

from classes.aircraft_2 import Aircraft, loader
from class1 import c1_m


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

def compute_aerodynamics(ac: Aircraft) -> Aircraft:
    """
    Step 1 — update aerodynamic coefficients and L/D from current geometry.
    Replace the body with your aero model.
    """
    # ── placeholder logic ────────────────────────────────────────────────────
    # e.g. ac.wing.ld = your_aero_module.compute_ld(ac.wing, ac.mission)
    # e.g. ac.wing.c_f = your_aero_module.compute_cf(ac.wing)
    pass
    # ─────────────────────────────────────────────────────────────────────────
    return ac


def compute_propulsion(ac: Aircraft) -> Aircraft:
    """
    Step 2 — update engine/propulsion parameters.
    Replace the body with your propulsion model.
    """
    # ── placeholder logic ────────────────────────────────────────────────────
    # e.g. required_power = your_prop_module.shaft_power(ac)
    # e.g. ac.engine.Phi   = your_prop_module.hybridisation_ratio(ac)
    pass
    # ─────────────────────────────────────────────────────────────────────────
    return ac


def compute_energy_weights(ac: Aircraft) -> Aircraft:
    """
    Step 3 — size fuel / battery from updated mission & propulsion.
    Replace the body with your energy model.
    """
    # ── placeholder logic ────────────────────────────────────────────────────
    # e.g. ac.weights.m_fuel    = your_energy_module.fuel_mass(ac)
    # e.g. ac.weights.m_battery = your_energy_module.battery_mass(ac)
    # e.g. ac.weights.m_energy  = ac.weights.m_fuel + ac.weights.m_battery
    pass
    # ─────────────────────────────────────────────────────────────────────────
    return ac


def compute_weights(ac: Aircraft) -> Aircraft:
    """
    Step 4 — update MTOW and empty weight from component contributions.
    Replace the body with your weight model.
    """
    # ── placeholder logic ────────────────────────────────────────────────────
    # e.g. ac.weights.m_empty   = your_weight_module.empty_mass(ac)
    # e.g. ac.weights.m_takeoff = (ac.weights.m_empty
    #                              + ac.weights.m_payload
    #                              + ac.weights.m_energy)
    pass
    # ─────────────────────────────────────────────────────────────────────────
    return ac


def compute_wing_geometry(ac: Aircraft) -> Aircraft:
    """
    Step 5 — re-size wing area and span from updated MTOW.
    Replace the body with your wing sizing model.
    """
    # ── placeholder logic ────────────────────────────────────────────────────
    # wing_loading = your_sizing_module.wing_loading(ac)
    # ac.wing.area = ac.weights.m_takeoff * 9.81 / wing_loading
    # ac.wing.span = (ac.wing.aspect_ratio * ac.wing.area) ** 0.5
    pass
    # ─────────────────────────────────────────────────────────────────────────
    return ac


def compute_fuselage(ac: Aircraft) -> Aircraft:
    """
    Step 6 — update fuselage geometry if driven by weights / payload.
    Replace the body with your fuselage model.
    """
    # ── placeholder logic ────────────────────────────────────────────────────
    pass
    # ─────────────────────────────────────────────────────────────────────────
    return ac

def compute_class_I_mass(ac: Aircraft) -> Aircraft:
    """Initial mass estimate from wing loading, T/W, statistical methods etc."""
    ac.weights.m_takeoff = ...
    oew_frac = c1_m.operating_empty_frac(ac)

    result = c1_m.energy_frac_needed(ac)

    if isinstance(result, tuple) and len(result) == 2:
        fuel_frac, bat_frac = result
        energy_frac = fuel_frac + bat_frac
        pl_frac = 1 - oew_frac - bat_frac - fuel_frac
    else:
        bat_frac = result
        fuel_frac = 0.0
        pl_frac = 1 - oew_frac - bat_frac
        energy_frac = bat_frac
    
    ac.weights.m_empty = 
    ac.weights.m_energy = sum(c1_m.energy_frac_needed(ac)) * ac.weights.m_takeoff
    ac.weights.m_fuel = c1_m.energy_frac_needed(ac)[0]
    ac.weights.m_battery = 
    ac.weights.m_payload = ac.weights.m_takeoff - 
    return ac

def compute_class_II_mass(ac: Aircraft) -> Aircraft:
    """Component build-up from actual geometry computed this epoch."""
    ac.weights.m_takeoff = ...   # overwrites class I estimate
    ac.weights.m_oew     = ...
    return ac

ITERATION_STEPS: list[Callable[[Aircraft], Aircraft]] = [
    compute_aerodynamics,
    compute_wing_geometry,
    compute_fuselage,
    compute_propulsion,
    compute_energy_weights
]


# 3. Define single iteration cycle

def run_iteration(ac: Aircraft,
                  INNER_TOLERANCE: float = 0.01  # 1 % standard for Class I and II OEW convergence tolerance
                  ) -> tuple[Aircraft, bool]:
    """
    Returns the updated aircraft AND whether the inner mass loop
    (class I vs class II) converged this epoch.
    """

    # Class I mass
    ac = compute_class_I_mass(ac)
    mtow_class_I = ac.weights.m_takeoff      # snapshot before geometry steps

    # Geometry / aero / propulsion etc.
    """Execute every registered step in order."""
    for step in ITERATION_STEPS:
        ac = step(ac)

    # Class II mass 
    ac = compute_class_II_mass(ac)
    mtow_class_II = ac.weights.m_takeoff

    # ── Inner convergence check ───────────────────────────────────────────
    inner_converged = (
        abs(mtow_class_II - mtow_class_I) / mtow_class_I
    ) < INNER_TOLERANCE

    return ac, inner_converged



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

    for epoch in range(1, config.max_epochs + 1):

        # Deep-copy BEFORE the iteration so we can diff afterwards
        prev = copy.deepcopy(ac)

        # Run one full design cycle
        ac, inner_converged = run_iteration(ac, INNER_TOLERANCE)

        # Check convergence
        converged, deltas = has_converged(prev, ac, config.convergence_params)

        # Log to file (read-only history — never read back by the loop)
        record = _snapshot(epoch, ac, deltas, inner_converged)
        _append_epoch(record, config.history_file)

        # Console output
        _print_epoch(epoch, ac, deltas, config.convergence_params, converged, inner_converged)

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
    import sys

    yaml_path = sys.argv[1] if len(sys.argv) > 1 else "aircraft.yaml"

    ac = loader.load(yaml_path, Aircraft)

    config = DesignLoopConfig(
        max_epochs   = 30,
        history_file = "aircraft_history.json",
        convergence_params = [
            ConvergenceParam("weights.m_empty", 1.0),
            ConvergenceParam("wing.area",         0.01),
            ConvergenceParam("wing.ld",           0.001),
        ],
    )

    final_aircraft = run_design_loop(ac, config)
    print("\nFinal aircraft state:")
    print(final_aircraft)