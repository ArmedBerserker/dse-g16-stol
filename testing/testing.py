import os, sys
ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..')
)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from class1 import c1_m as c1
from classes import aircraft_2 as ac
import matplotlib.pyplot as plt
import pandas as pd
from class1 import c1_matching_comparison as match
import numpy as np
from lookups.consts import *
from class1 import prelim_drag as drag

def explore_concepts(eng_paths : str,
                     wing_paths : str,
                     reqs_paths : str,
                     fuse_path : str,
                     mission_path : str,
                     weights_path : str,
                     do_powers : bool, 
                     do_sensitivity : bool):
    ############## THIS WHOLE SEGMENT CALCULATES THE MASS FRACTIONS FOR EACH ONE, ASSUMING CRUISE SIZING
    ac_dict = {}
    match_cr = {}
    match_to = {}
    sens_dict = {}
    gear = ("Taildragger", "Tricycle")
    engines = ("H2", "Boosted Piston", "Piston Hybrid", "Boosted Turboprop", "Turbine Hybrid")
    for i, e in enumerate(eng_paths):
        for j, w in enumerate(wing_paths):
            eng = ac.loader.load(e[0], ac.Engine)
            fuse = ac.loader.load(fuse_path, ac.Fuselage)
            mission = ac.loader.load(mission_path, ac.Mission)
            reqs = ac.loader.load(reqs_paths[1], ac.Requirements)
            weights = ac.loader.load(weights_path, ac.Weights)
            wing = ac.loader.load(w[0], ac.Wing)

            # name = str(e).split('/')[1].split('.')[0] + ' ' + str(w).split('/')[1].split('.')[0]
            name = gear[j] + " " + engines[i]
            aircraft = ac.Aircraft(name = name, 
                                requirements = reqs, 
                                mission = mission, 
                                weights = weights, 
                                wing = wing, 
                                fuselage = fuse, 
                                engine = eng)
            
            fuel_frac = c1.energy_frac_needed(aircraft)
            oem = c1.operating_empty_frac(aircraft, 
                                        source_for_fracs='specific',
                                        engine_type=e[1], 
                                        gear_type=w[1])
            ac_dict[name] = {'Energy Frac': sum(fuel_frac),
                            'Operating Empty Mass' : oem,
                            'Engine Type' : e[1],
                            'Gear Type' : w[1]}
            
            cd0 = drag.cd0(aircraft, 'Twin Engine Propeller Driven')
            k, eff = drag.k(aircraft)
            print('###############################################################################')
            print(f'Aircraft: {name}\nEnergy Fraction is: {fuel_frac}\nOperating Empty Mass is: {oem}\nZero-lift drag is: {cd0}\nOswald efficiency is: {eff}\n')
            print(f'Drag is: {drag.prelim_drag(aircraft, 'Twin Engine Propeller Driven')}')
            if ('H2' in name or 'Boosted' in name) and do_powers:
                data_cr = match.plot_matching_and_select_design_point(aircraft, 'Twin Engine Propeller Driven', 
                                                                    W_S_plot = np.arange(0.1, 1000), 
                                                                    W_P_plot = np.arange(1e-8, 1e-1, 1e-4), 
                                                                    requirement_to_meet = 'cruise', show_plot = False)
                data_to = match.plot_matching_and_select_design_point(aircraft, 'Twin Engine Propeller Driven',
                                                                    W_S_plot = np.arange(0.1, 1000),
                                                                    W_P_plot = np.arange(1e-8, 1e-1, 1e-4),
                                                                    requirement_to_meet = 'to', show_plot = False)
                match_cr[name] = data_cr
                match_to[name] = data_to
                print(f'Drag for this config is {drag.prelim_drag(aircraft, type_to_use= 'Twin Engine Propeller Driven')}')
            if do_sensitivity:
                sens_dict[name] = match.run_sensitivity_study_save_results(aircraft,
                                                          W_S_plot = np.arange(0.1, 1000),
                                                          W_P_or_T_W_plot = np.arange(1e-8, 1e-1, 1e-4))
    df = pd.DataFrame(ac_dict).T
    df_cr = pd.DataFrame(match_cr).T
    df_to = pd.DataFrame(match_to).T
        

    return df, df_cr, df_to, sens_dict

def main(masses = True, powers = True, sensitivity = True, plots = False):
    eng_paths = [('concepts/engine_h2.yaml', 'turboprop'),                 # Hydrogen Engine
                ('concepts/engine_piston_b.yaml', 'piston'),           # Piston Engine + Booster
                ('concepts/engine_piston_e.yaml', 'piston'),           # Piston Engine as generator
                ('concepts/engine_tprop_b.yaml', 'turboprop'),            # Turboprop engine + Booster
                ('concepts/engine_turb_e.yaml', 'turboprop')]             # Turbine engine as generator

    wing_paths = [('concepts/wing_courier.yaml', 'tail dragger'),             # Wing for courier-like config
                    ('concepts/wing_electra.yaml', 'tricycle')]             # Wing for electra-like config

    reqs_paths = ['concepts/reqs_nturb.yaml',
                    'concepts/reqs_turb.yaml']                        # THING

    fuse_path = 'yamls/fuselage.yaml'
    mission_path = 'yamls/mission.yaml'
    weights_path = 'yamls/weights.yaml'

    if masses:
        results = explore_concepts(
        eng_paths,
        wing_paths,
        reqs_paths,
        fuse_path,
        mission_path,
        weights_path,
        False,
        False)[0]


        mass_fraction_cols = [
        "Operating Empty Mass",
        "Energy Frac",
        ]

        plot_df = results[mass_fraction_cols].astype(float)

        fig, ax = plt.subplots(figsize=(13, 6.5))

        plot_df.plot(
            kind="bar",
            stacked=True,
            ax=ax,
            width=0.72,
            edgecolor="black",
            linewidth=0.6,
        )

        ax.set_title(
            "Mass Fraction Breakdown by Aircraft Configuration",
            fontsize=15,
            weight="bold",
            pad=15,
        )

        ax.set_xlabel("")
        ax.set_ylabel("Mass Fraction", fontsize=12)

        ax.set_ylim(0, plot_df.sum(axis=1).max() * 1.15)

        ax.tick_params(axis="x", labelrotation=35, labelsize=10)
        ax.tick_params(axis="y", labelsize=10)

        ax.set_xticklabels(
            plot_df.index,
            ha="right",
        )

        ax.grid(
            axis="y",
            linestyle="--",
            alpha=0.4,
        )

        ax.axhline(
            y=0.65,
            color="gray",
            linestyle=":",
            linewidth=1.5,
            alpha=0.8,
        )

        ax.axhline(
            y=1.0,
            color="red",
            linestyle=":",
            linewidth=1.5,
            alpha=0.8,
        )
        ax.legend(
            title="Component",
            frameon=True,
            loc="upper left",
            bbox_to_anchor=(1.01, 1),
        )

        # Add total mass-fraction labels above each stacked bar
        totals = plot_df.sum(axis=1)

        for i, total in enumerate(totals):
            ax.text(
                i,
                total + 0.01,
                f"{total:.3f}",
                ha="center",
                va="bottom",
                fontsize=9,
            )

        plt.tight_layout()
        plt.savefig('outputs/class_1configs10.png')

        eng = ac.loader.load(eng_paths[2][0], ac.Engine)
        fuse = ac.loader.load(fuse_path, ac.Fuselage)
        mission = ac.loader.load(mission_path, ac.Mission)
        reqs = ac.loader.load(reqs_paths[0], ac.Requirements)
        weights = ac.loader.load(weights_path, ac.Weights)
        wing = ac.loader.load(wing_paths[0][0], ac.Wing)

        # name = str(e).split('/')[1].split('.')[0] + ' ' + str(w).split('/')[1].split('.')[0]
        name = "Hybrid"
        aircraft = ac.Aircraft(name = name, 
                            requirements = reqs, 
                            mission = mission, 
                            weights = weights, 
                            wing = wing, 
                            fuselage = fuse, 
                            engine = eng)
        
        Phi = [0.1 * n for n in range(1, 10)]

        frac_evol = [
            sum(c1.breguet_hyb(aircraft, p))
            for p in Phi
        ]

        oem = c1.operating_empty_frac(aircraft)

        overall_evol = [
            f + oem
            for f in frac_evol
        ]

        fig, ax = plt.subplots(figsize=(8, 5))

        ax.plot(
            Phi,
            overall_evol,
            marker="o",
            linewidth=2,
            markersize=6,
        )

        ax.set_title(
            "Piston-Electric Aircraft Total Mass Fraction vs Phi",
            fontsize=14,
            weight="bold",
            pad=12,
        )

        ax.set_xlabel(r"$\Phi$", fontsize=12)
        ax.set_ylabel("Total Mass Fraction: OEM + Energy", fontsize=12)

        ax.grid(
            axis="both",
            linestyle="--",
            alpha=0.4,
        )

        for phi, total_frac in zip(Phi, overall_evol):
            ax.text(
                phi,
                total_frac + 0.002,
                f"{total_frac:.3f}",
                ha="center",
                va="bottom",
                fontsize=9,
            )

        ax.tick_params(axis="both", labelsize=10)

        plt.tight_layout()
        plt.savefig("outputs/piston_electric_total_mass_fraction_phi.svg")
        #plt.show()
    if powers:
        _, cr_results, to_results, JUNK = explore_concepts(eng_paths, 
                                                     wing_paths,
                                                     reqs_paths, 
                                                     fuse_path, 
                                                     mission_path, 
                                                     weights_path,
                                                     True,
                                                     False)

        wp_difference = cr_results["W/P"] - to_results["W/P"]
        print(wp_difference)

        to_time = 200 / (1.05 * 50 * KTS_TO_MS / np.sqrt(2))
        print(to_time)

        #### MASS IS HARD CODED RN!!! #####
        weights = ac.loader.load('yamls/weights.yaml', ac.Weights)
        engine = ac.loader.load('concepts/engine_turb_e.yaml', ac.Engine)
        power_needed = (weights.m_takeoff * g) / wp_difference
        battery_needed = power_needed * to_time / engine.e_2 
        
        # --------------------------------------------------
        # Plot battery needed
        # --------------------------------------------------

        battery_needed = battery_needed.astype(float)

        fig, ax = plt.subplots(figsize=(11, 6))

        ax.bar(
            battery_needed.index,
            battery_needed.values,
            width=0.7,
            edgecolor="black",
            linewidth=0.6,
        )

        ax.set_title(
            "Battery Mass Needed for Takeoff Power Deficit",
            fontsize=14,
            weight="bold",
            pad=12,
        )

        ax.set_xlabel("Aircraft Configuration", fontsize=12)
        ax.set_ylabel("Battery Needed", fontsize=12)

        ax.tick_params(axis="x", rotation=35, labelsize=10)
        ax.tick_params(axis="y", labelsize=10)

        ax.set_xticklabels(
            battery_needed.index,
            ha="right",
        )

        ax.grid(
            axis="y",
            linestyle="--",
            alpha=0.4,
        )

        for i, value in enumerate(battery_needed.values):
            ax.text(
                i,
                value,
                f"{value:.1f}",
                ha="center",
                va="bottom",
                fontsize=9,
            )

        plt.tight_layout()
        plt.savefig("outputs/battery_needed_takeoff.svg")
        #plt.show()

        # if sensitivity:
        #     results = explore_concepts()
    if sensitivity:
        _, _, _, thing = explore_concepts(eng_paths, wing_paths, reqs_paths, fuse_path, mission_path, weights_path, False, True)
        print(thing)
        
if __name__ == '__main__':
    main(masses = False,
         powers = True,
         sensitivity = False, 
         plots = False)
    