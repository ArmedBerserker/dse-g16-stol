import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from classes.aircraft_2 import loader, Aircraft
from lookups.consts import *
import numpy as np
import pandas as pd


def size_tires(ac: Aircraft):
    W_to = ac.weights.m_takeoff * g
    n_min_nlg = ac.landing_gear.n_nlg_min_as
    n_max_nlg = ac.landing_gear.n_nlg_max_as

    f_mlg = 1 - n_min_nlg
    f_nlg = n_max_nlg
    Ft_mlg = 0.5 * f_mlg * W_to * 1.25  # [N] Main landing gear load requirement
    Ft_nlg = f_nlg * W_to * 1.25  # [N] Nose/tail landing gear load requirement

    tire_to_speed = 1.1 * ac.requirements.take_off['to_speed'] # [Kts]
    tire_app_speed = 1.2 * ac.requirements.approach['app_speed'] # [Kts]
    V_max_tire = max(tire_to_speed, tire_app_speed)

    V_max_tire = V_max_tire * KTS_TO_MPH  # [Mph]
    pt = ac.landing_gear.pt             # [Psi]
    Load_mlg = Ft_mlg * 1/(g*LBS_TO_KG)   # [Lbs] per tire
    Load_nlg = Ft_nlg * 1/(g*LBS_TO_KG)   # [Lbs]

    tire_data = pd.read_csv('lookups/type3_tires.csv')

    # Converting speeds given in Kts to Mph
    tire_data["Rated Speed (MPH)"] = tire_data["Rated Speed (MPH)"].astype(str).apply(
        lambda x: float(x.replace("K", "")) * KTS_TO_MPH if x.endswith("K") else float(x))

    def filter_tires(db: pd.DataFrame, req_load: float) -> pd.DataFrame:
        mask = (
                (db["Rated Load (Lbs)"] >= req_load)  # load requirement
                & (db["Rated Inflation (PSI)"] >= pt)  # pressure requirement
                & (db["Rated Speed (MPH)"] >= V_max_tire)  # speed requirement
        )
        return db[mask].copy()

    tire_options_mlg = filter_tires(tire_data, Load_mlg)
    tire_options_nlg = filter_tires(tire_data, Load_nlg)

    def pick_smallest_tire(tire_options: pd.DataFrame) -> pd.Series:
        if tire_options.empty:
            raise ValueError("No Type III tire found.")
        # Tires sorted by outside diameter first, then by section width as tiebreaker
        tire_options_sorted = tire_options.sort_values(by=["Outside Diameter Max (In)", "Section Width Max (In)"])
        return tire_options_sorted.iloc[0]  # first row = smallest qualifying tire

    best_mlg_tire = pick_smallest_tire(tire_options_mlg)
    best_nlg_tire = pick_smallest_tire(tire_options_nlg)

    ac.landing_gear.selected_mlg_tire = best_mlg_tire.to_dict()
    ac.landing_gear.selected_nlg_tire = best_nlg_tire.to_dict()


def tire_location(ac: Aircraft):
    # For taildragger nlg = tail landing gear

    gear_type = ac.landing_gear.gear_type
    n_min_nlg = ac.landing_gear.n_nlg_min_as
    n_max_nlg = ac.landing_gear.n_nlg_max_as
    tipover = ac.landing_gear.tipover * np.pi/180             # [rad]
    scrape = ac.landing_gear.scrape * np.pi/180               # [rad]
    turnover = ac.landing_gear.turnover * np.pi/180           # [rad]
    bank = ac.landing_gear.bank * np.pi/180                   # [rad]
    prop_clear = ac.landing_gear.prop_clear  # [m]

    # The coordinate frame assumes [0,0] = [aircraft nose, fuselage bottom]. Z upwards
    X_cg_fwd = ac.weights.cg_fwd   # [m]
    X_cg_aft = ac.weights.cg_aft   # [m]
    Z_fus = 0                      # [m]
    fus_pitch = ac.landing_gear.fus_pitch * np.pi/180    # [rad]
    fus_ground_clear = ac.landing_gear.fus_ground_clear  # [m]
    Z_cg = ac.weights.cg_height    # [m]
    X_tcone = ac.fuselage.x_tcone  # [m]
    Z_tcone = 0                    # [m]
    Z_prop = ac.engine.z_prop      # [m]
    Y_prop = ac.engine.y_prop      # [m] for the outermost wing-mounted propeller
    X_prop = ...                   # [m] nose to prop
    bw = ac.wing.span              # [m]
    Z_tip = ...                    # [m]
    X_tip = ...                    # [m]
    a = ac.landing_gear.a
    s = ac.landing_gear.s          # [m]

    mlg_tire = ac.landing_gear.selected_mlg_tire
    nlg_tire = ac.landing_gear.selected_nlg_tire

    D_tire = mlg_tire["Outside Diameter Max (In)"] * IN_TO_M  # [m]
    D_rim = mlg_tire["Specified Rim Diameter (In)"] * IN_TO_M  # [m]

    # LONGITUDINAL POSITIONING
    if gear_type == 'tricycle':
        # Tipover constraint line
        x1 = [X_cg_aft, X_cg_aft - np.tan(tipover)]
        z1 = [Z_cg, Z_cg + 1]
        slope1, intercept1 = np.polyfit(x1, z1, 1)

        # Scrape constraint line
        x2 = [X_tcone, X_tcone + 1]
        z2 = [Z_tcone, Z_tcone + np.tan(scrape)]
        slope2, intercept2 = np.polyfit(x2, z2, 1)

        # Find intersection point
        X_mlg = (intercept2 - intercept1) / (slope1 - slope2)
        Z_mlg = slope1 * X_mlg + intercept1
        X_nlg_fwd_lim = X_cg_aft - (1/n_min_nlg - 1)*(X_mlg - X_cg_aft)
        X_nlg_aft_lim = X_cg_fwd - (1/n_max_nlg - 1)*(X_mlg - X_cg_fwd)
        # location of nlg should be chosen as the most fwd possible one if structurally allowed.
        X_nlg = X_nlg_fwd_lim

        # Need to check propeller clearance req. met
        if Z_prop - Z_mlg < prop_clear:
            extra_height_req = Z_prop - Z_mlg - prop_clear
            Z_mlg = Z_mlg + extra_height_req

        Z_nlg = Z_mlg

        # LATERAL POSITIONING
        l_mlg_fwd = abs(X_mlg - X_cg_fwd)
        l_nlg_fwd = abs(X_cg_fwd - X_nlg)
        Y_mlg_req1 = (l_mlg_fwd + l_nlg_fwd) / np.sqrt(
            l_nlg_fwd ** 2 * np.tan(turnover) ** 2 / (Z_cg - Z_mlg) ** 2 - 1)
        Y_mlg_req2 = Y_prop - (Z_prop - Z_mlg) / np.tan(bank)
        Y_mlg_req3 = bw / 2 - ((Z_tip - Z_mlg) + (X_mlg - X_tip) * np.sin(scrape)) / np.tan(bank)
        delta_Z_mlg = (a * s + (D_tire - D_rim) / 2)
        Y_mlg_req4 = delta_Z_mlg / (2 * (Z_prop - Z_mlg) - delta_Z_mlg) * Y_prop
        min_Y_mlg = max(Y_mlg_req1, Y_mlg_req2, Y_mlg_req3, Y_mlg_req4)
        # Lateral location on mlg should be chosen to be as close as possible to the primary structure
        Y_mlg = min_Y_mlg
        Y_nlg = 0


    if gear_type == 'taildragger':
        # Tipover constraint line
        x1 = [X_cg_fwd, X_cg_fwd + 1]
        z1 = [Z_cg, Z_cg + np.tan(np.pi/2 + fus_pitch - tipover)]
        slope1, intercept1 = np.polyfit(x1, z1, 1)

        # Prop clearance line
        slope5 = np.sin(fus_pitch)
        intercept5 = Z_prop*np.cos(fus_pitch) - X_prop*np.sin(fus_pitch) - prop_clear

        # Fuselage pitch constraint line
        x2 = [X_tcone + fus_ground_clear*np.sin(fus_pitch), X_tcone + fus_ground_clear*np.sin(fus_pitch) + 1]
        z2 = [Z_tcone - fus_ground_clear*np.cos(fus_pitch), Z_tcone - fus_ground_clear*np.cos(fus_pitch) + np.tan(fus_pitch)]
        slope2, intercept2 = np.polyfit(x2, z2, 1)

        # Find Intersection
        X_mlg = (intercept2 - intercept1) / (slope1 - slope2)
        Z_mlg = slope1 * X_mlg + intercept1

        X_nlg_fwd_lim = X_cg_aft - (1 / n_max_nlg - 1) * (X_mlg - X_cg_aft)
        X_nlg_aft_lim = X_cg_fwd - (1 / n_min_nlg - 1) * (X_mlg - X_cg_fwd)
        # The tail landing gear position is chosen as the most aft possible one
        X_nlg = X_nlg_aft_lim
        Z_nlg = Z_mlg + np.tan(fus_pitch) * (X_nlg - X_mlg)

        # LATERAL POSITIONING
        l_mlg_aft = abs(X_mlg - X_cg_aft)
        l_nlg_aft = abs(X_cg_aft - X_nlg)
        Y_mlg_req1 = (l_mlg_aft + l_nlg_aft) / np.sqrt(
            l_nlg_aft ** 2 * np.tan(turnover) ** 2 / (Z_cg - Z_mlg) ** 2 - 1)
        Y_mlg_req2 = Y_prop - (Z_prop - Z_mlg) / np.tan(bank)
        # no req3
        delta_Z_mlg = (a * s + (D_tire - D_rim) / 2)
        Y_mlg_req4 = delta_Z_mlg / (2 * (Z_prop - Z_mlg) - delta_Z_mlg) * Y_prop
        min_Y_mlg = max(Y_mlg_req1, Y_mlg_req2, Y_mlg_req4)
        # Lateral location on mlg should be chosen to be as close as possible to the primary structure
        Y_mlg = min_Y_mlg
        Y_nlg = 0

    # ALL DIMENSIONS WRT TO NOSE & BOTTOM OF FUSELAGE, [X, Y, Z] = CENTRE POINT OF WHEELS
    # NOSE/TAIL LANDING GEAR LOCATION
    ac.landing_gear.longitudinal_nlg = X_nlg
    ac.landing_gear.lateral_nlg = 0
    ac.landing_gear.height_nlg = abs(Z_nlg)

    # MAIN LANDING GEAR LOCATION
    ac.landing_gear.longitudinal_mlg = X_mlg
    ac.landing_gear.lateral_mlg = Y_mlg
    ac.landing_gear.height_mlg = abs(Z_mlg)
