import os, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from classes import aircraft_2 as ac
from class1 import c1_m as mass
import matplotlib.pyplot as plt
from datetime import datetime

def explore_concept_masses(parent_file : str) -> tuple[float, float]:
    
    concept = ac.loader.load(parent_file, ac.Aircraft)
    name = concept.name 
    oem_eng_type = concept.engine.alpha_p_id
    config = concept.wing.type
    if config == 'A':
        gear = 'tricycle'
    elif config == 'B':
        gear = 'tail dragger'
    else:
        raise ValueError('wing config is incorrect')
    # OUTPUT 1
    oem_frac = mass.operating_empty_frac(ac = concept,
                                         source_for_fracs = 'specific',
                                         engine_type = oem_eng_type,
                                         gear_type = gear)

    # OUTPUT 2
    energy_frac = mass.energy_frac_needed(ac = concept,)
    
    return (concept.name, oem_frac, energy_frac)

def run_loop(path_list : list[str]) -> dict:
    exploration = {}
    for p in path_list:
        name, oem, energy = explore_concept_masses(p)
        if isinstance(energy, tuple):
            energy = sum(energy)
        
        exploration[name] = {'oem_frac' : oem,
                             'energy_frac' : energy}
    return exploration

def plot_exploration(thing : dict, thing2 : dict, show = True):
    names = list(thing.keys())
    names.extend(list(thing2.keys()))
    oem_fracs = [data['oem_frac'] for data in thing.values()]
    oem_fracs.extend([data['oem_frac'] for data in thing2.values()])
    energy_fracs = [data['energy_frac'] for data in thing.values()]
    energy_fracs.extend([data['energy_frac'] for data in thing2.values()])
    fig, ax = plt.subplots(figsize=(7, 5))

    bars_oem = ax.bar(names,
           oem_fracs,
           color = 'tab:blue',
           label = r'$\frac{m_{OE}}{m_{TO}}$')
    bars_energy = ax.bar(names,
           energy_fracs,
           bottom = oem_fracs,
           color = 'tab:orange',
           label = r'$\frac{m_{e}}{m_{TO}}$')
    
    # ax.bar_label(bars_oem,
    #              fmt  = '%.3f',
    #              label_type = 'center')
    
    # ax.bar_label(bars_energy,
    #              fmt  = '%.3f',
    #              label_type = 'center',
    #              fontsize = 8)
    ax.axhline(y = 1.0,
               color = 'tab:red',
               linestyle = 'solid',
               linewidth = 2,
               label = r'$m_{TO}$ limit')
    
    ax.axhline(y = 1.0 - (704/2000),
               color = 'tab:red',
               linestyle = 'dashed',
               linewidth = 2,
               label = r'$m_{PL}$ limit (0.35)')
    ax.axhline(y = 1.0 - (704/1870),
               color = 'tab:red',
               linestyle = 'dotted',
               linewidth = 2,
               label = r'$m_{TO}$ limit with contingency (0.38)')
    
    ax.grid(axis = 'y', linestyle = 'dashed', alpha = 0.4)
    ax.set_axisbelow(True)

    ax.set_xlabel('Configuration ID')
    ax.set_ylabel('Mass Fraction ($m_{OE}+m_e$)')

    ax.legend(loc = 'upper right')
    plt.tight_layout()
    plt.savefig(f'testing/results/fractions_{datetime.now():%Y%m%d_%H%M%S}.png', dpi = 300)

    if show:
        plt.show()

def log_exploration(thing : dict, thing2 : dict):
    names = list(thing.keys())
    names.extend(list(thing2.keys()))
    oem_fracs = [data['oem_frac'] for data in thing.values()]
    oem_fracs.extend([data['oem_frac'] for data in thing2.values()])
    energy_fracs = [data['energy_frac'] for data in thing.values()]
    energy_fracs.extend([data['energy_frac'] for data in thing2.values()])
    with open(f'testing/results/fractions_{datetime.now():%Y%m%d_%H%M%S}.txt', 'w') as f:
        for i in range(len(names)):
            f.write(
            f"Configuration {names[i]}:\n"
            f"\tOperating Empty Mass fraction: {oem_fracs[i]:.4f}\n"
            f"\tEnergy Mass fraction: {energy_fracs[i]:.4f}\n"
            f"\tTotal fraction: {oem_fracs[i] + energy_fracs[i]:.4f}\n\n"
        )
if __name__ == '__main__':
    
    A_paths = ['concepts2/A1P.yaml',
               'concepts2/A1T.yaml',
               'concepts2/A2H.yaml',
               'concepts2/A3P.yaml',
               'concepts2/A3T.yaml']
    
    B_paths = ['concepts2/B1P.yaml',
               'concepts2/B1T.yaml',
               'concepts2/B2H.yaml',
               'concepts2/B3P.yaml',
               'concepts2/B3T.yaml']

    A_dict = run_loop(A_paths)
    B_dict = run_loop(B_paths)
    print(B_dict['B3T'])
    plot_exploration(A_dict, B_dict, False)
    log_exploration(A_dict, B_dict)
   
