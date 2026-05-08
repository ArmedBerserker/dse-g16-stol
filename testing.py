import class1.c1_m as c1
from classes import aircraft_2 as ac

eng_paths = ['concepts/engine_h2.yaml',                 # Hydrogen Engine
             'concepts/engine_piston_b.yaml',           # Piston Engine + Booster
             'concepts/engine_piston_e.yaml',           # Piston Engine as generator
             'concepts/engine_tprop_b.yaml',            # Turboprop engine + Booster
             'concepts/engine_turb_e.yaml']             # Turbine engine as generator

wing_paths = ['concepts/wing_courier.yaml',             # Wing for courier-like config
              'concepts/wing_electra.yaml']             # Wing for electra-like config

reqs_paths = ['yamls/reqs.yaml']                        # THING

fuse_path = 'yamls/fuselage.yaml'
mission_path = 'yamls/mission.yaml'
weights_path = 'yamls/weights.yaml'

# eng = ac.loader.load(eng_paths[3], ac.Engine)
# fuse = ac.loader.load(fuse_path, ac.Fuselage)
# mission = ac.loader.load(mission_path, ac.Mission)
# reqs = ac.loader.load(reqs_paths[0], ac.Requirements)
# weights = ac.loader.load(weights_path, ac.Weights)
# wing = ac.loader.load(wing_paths[0], ac.Wing)

# aircraft = ac.Aircraft('Config 1A', requirements = reqs, mission = mission, weights = weights, wing = wing, fuselage = fuse, engine = eng)

# fuel_frac = c1.energy_frac_needed(aircraft)
# oem = c1.operating_empty_frac(aircraft)
# print(fuel_frac, oem)

for e in eng_paths:
    for w in wing_paths:
        eng = ac.loader.load(e, ac.Engine)
        fuse = ac.loader.load(fuse_path, ac.Fuselage)
        mission = ac.loader.load(mission_path, ac.Mission)
        reqs = ac.loader.load(reqs_paths[0], ac.Requirements)
        weights = ac.loader.load(weights_path, ac.Weights)
        wing = ac.loader.load(w, ac.Wing)

        name = str(e).split('/')[1].split('.')[0] + ' ' + str(w).split('/')[1].split('.')[0]
        aircraft = ac.Aircraft(name, requirements = reqs, mission = mission, weights = weights, wing = wing, fuselage = fuse, engine = eng)

        fuel_frac = c1.energy_frac_needed(aircraft)
        oem = c1.operating_empty_frac(aircraft)
        print('###############################################################################')
        print(f'Aircraft: {name}\nEnergy Fraction is: {fuel_frac}\nOperating Empty Mass is: {oem}\n')