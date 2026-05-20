import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tracemalloc

tracemalloc.start()
dataset_path = 'lookups/ref.csv'

df = pd.read_csv(dataset_path)

turbo_prop_mtow = df.loc[df['Is turboprop'] == True, 'MTOW [kg]'].to_numpy()
turbo_prop_oem = df.loc[df['Is turboprop'] == True, 'OEW [kg]'].to_numpy()
a_turbo_prop, b_turbo_prop = np.polyfit(turbo_prop_mtow, turbo_prop_oem, 1)


piston_mtow = df.loc[df['Is turboprop'] == False, 'MTOW [kg]'].to_numpy()
piston_oem = df.loc[df['Is turboprop'] == False, 'OEW [kg]'].to_numpy()
a_piston, b_piston = np.polyfit(piston_mtow, piston_oem, 1)


tail_dragger_mtow = df.loc[df['Gear config'] == 'taildragger', 'MTOW [kg]'].to_numpy()
tail_dragger_oem = df.loc[df['Gear config'] == 'taildragger', 'OEW [kg]'].to_numpy()
a_tail_dragger, b_tail_dragger = np.polyfit(tail_dragger_mtow, tail_dragger_oem, 1)


tricycle_mtow = df.loc[df['Gear config'] == 'tricycle', 'MTOW [kg]'].to_numpy()
tricycle_oem = df.loc[df['Gear config'] == 'tricycle', 'OEW [kg]'].to_numpy()
a_tricycle, b_tricycle = np.polyfit(tricycle_mtow, tricycle_oem, 1)


def plot_regression(x, y, filename, colour):
    a, b = np.polyfit(x, y, 1)

    x_vals = np.linspace(np.min(x), np.max(x), 100)
    y_pred = a * x + b

    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - (ss_res / ss_tot)

    plt.figure(figsize=(7, 5))
    plt.scatter(x, y, c=colour, marker='x')
    plt.plot(x_vals, a * x_vals + b, colour)

    plt.xlabel(r'$m_{TO}$ [kg]')
    plt.ylabel(r'$m_{OE}$ [kg]')

    plt.text(
        0.05, 0.90,
        rf'$m_{{OE}} = {a:.3f}m_{{TO}} + {b:.1f}$' + '\n' + rf'$R^2 = {r2:.3f}$',
        transform=plt.gca().transAxes,
        bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3')
    )

    plt.tight_layout()
    plt.savefig(filename)
    #plt.show()

    return a, b, r2

# Propulsion: red
a_turbo_prop, b_turbo_prop, r2_turbo_prop = plot_regression(
    turbo_prop_mtow,
    turbo_prop_oem,
    'regression/regression_turboprop.png',
    'r'
)

a_piston, b_piston, r2_piston = plot_regression(
    piston_mtow,
    piston_oem,
    'regression/regression_piston.png',
    'r'
)

# Configuration: blue
a_tail_dragger, b_tail_dragger, r2_tail_dragger = plot_regression(
    tail_dragger_mtow,
    tail_dragger_oem,
    'regression/regression_taildragger.png',
    'b'
)

a_tricycle, b_tricycle, r2_tricycle = plot_regression(
    tricycle_mtow,
    tricycle_oem,
    'regression/regression_tricycle.png',
    'b'
)