import class1.c1_m as c1
from classes import aircraft_2 as ac
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

eng_paths = [('concepts/engine_h2.yaml', 'turboprop'),                 # Hydrogen Engine
             ('concepts/engine_piston_b.yaml', 'piston'),           # Piston Engine + Booster
             ('concepts/engine_piston_e.yaml', 'piston'),           # Piston Engine as generator
             ('concepts/engine_tprop_b.yaml', 'turboprop'),            # Turboprop engine + Booster
             ('concepts/engine_turb_e.yaml', 'turboprop')]             # Turbine engine as generator

wing_paths = [('concepts/wing_courier.yaml', 'tail dragger'),             # Wing for courier-like config
              ('concepts/wing_electra.yaml', 'tricycle')]             # Wing for electra-like config

reqs_paths = ['yamls/reqs.yaml']                        # THING

fuse_path = 'yamls/fuselage.yaml'
mission_path = 'yamls/mission.yaml'
weights_path = 'yamls/weights.yaml'


ac_dict = {}

for e in eng_paths:
    for w in wing_paths:
        eng = ac.loader.load(e[0], ac.Engine)
        fuse = ac.loader.load(fuse_path, ac.Fuselage)
        mission = ac.loader.load(mission_path, ac.Mission)
        reqs = ac.loader.load(reqs_paths[0], ac.Requirements)
        weights = ac.loader.load(weights_path, ac.Weights)
        wing = ac.loader.load(w[0], ac.Wing)

        name = str(e).split('/')[1].split('.')[0] + ' ' + str(w).split('/')[1].split('.')[0]

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

        print('###############################################################################')
        print(f'Aircraft: {name}\nEnergy Fraction is: {fuel_frac}\nOperating Empty Mass is: {oem}\n')


df = pd.DataFrame(ac_dict).T

plt.figure(1, figsize=(12, 6))

energy = df['Energy Frac'].astype(float)
oem = df['Operating Empty Mass'].astype(float)

plt.bar(
    df.index,
    oem,
    label='Operating Empty Mass'
)

plt.bar(
    df.index,
    energy,
    bottom=oem,
    label='Energy Fraction'
)

plt.xlabel('Aircraft design')
plt.ylabel('Mass fraction')
plt.title('Stacked mass fractions for all designs')

plt.xticks(rotation=45, ha='right')
plt.legend()
plt.grid(axis='y')

plt.figure(2)

eng = ac.loader.load('concepts/engine_piston_e.yaml', ac.Engine) # choose a hybrid please!
fuse = ac.loader.load(fuse_path, ac.Fuselage)
mission = ac.loader.load(mission_path, ac.Mission)
reqs = ac.loader.load(reqs_paths[0], ac.Requirements)
weights = ac.loader.load(weights_path, ac.Weights)
wing = ac.loader.load('concepts/wing_courier.yaml', ac.Wing)

name = str(e).split('/')[1].split('.')[0] + ' ' + str(w).split('/')[1].split('.')[0]

aircraft = ac.Aircraft(name = name, 
                        requirements = reqs, 
                        mission = mission, 
                        weights = weights, 
                        wing = wing, 
                        fuselage = fuse, 
                        engine = eng)

Phi = [0.1*n for n in range(1, 10)]
frac_evol = [sum(c1.breguet_hyb(aircraft, p)) for p in Phi]

plt.scatter(Phi, frac_evol)
plt.show()

