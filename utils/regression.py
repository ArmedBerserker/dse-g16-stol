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


plt.figure(1, figsize=(12, 8))

# Turboprop
plt.subplot(221)
plt.scatter(turbo_prop_mtow, turbo_prop_oem)

a_turbo_prop, b_turbo_prop = np.polyfit(turbo_prop_mtow, turbo_prop_oem, 1)
x_vals = np.linspace(np.min(turbo_prop_mtow), np.max(turbo_prop_mtow), 100)
y_pred = a_turbo_prop * turbo_prop_mtow + b_turbo_prop

ss_res = np.sum((turbo_prop_oem - y_pred) ** 2)
ss_tot = np.sum((turbo_prop_oem - np.mean(turbo_prop_oem)) ** 2)
r2_turbo_prop = 1 - (ss_res / ss_tot)

plt.plot(x_vals, a_turbo_prop * x_vals + b_turbo_prop)
plt.title(r'$m_{oe}$ vs $m_{to}$ for turboprop aircraft')
plt.xlabel(r'$m_{to}$ [kg]')
plt.ylabel(r'$m_{oe}$ [kg]')
plt.text(
    0.05, 0.90,
    rf'$m_{{oe}} = {a_turbo_prop:.3f}m_{{to}} + {b_turbo_prop:.1f}$' + '\n' + rf'$R^2 = {r2_turbo_prop:.3f}$',
    transform=plt.gca().transAxes
)


# Piston
plt.subplot(222)
plt.scatter(piston_mtow, piston_oem)

a_piston, b_piston = np.polyfit(piston_mtow, piston_oem, 1)
x_vals = np.linspace(np.min(piston_mtow), np.max(piston_mtow), 100)
y_pred = a_piston * piston_mtow + b_piston

ss_res = np.sum((piston_oem - y_pred) ** 2)
ss_tot = np.sum((piston_oem - np.mean(piston_oem)) ** 2)
r2_piston = 1 - (ss_res / ss_tot)

plt.plot(x_vals, a_piston * x_vals + b_piston)
plt.title(r'$m_{oe}$ vs $m_{to}$ for piston aircraft')
plt.xlabel(r'$m_{to}$ [kg]')
plt.ylabel(r'$m_{oe}$ [kg]')
plt.text(
    0.05, 0.90,
    rf'$m_{{oe}} = {a_piston:.3f}m_{{to}} + {b_piston:.1f}$' + '\n' + rf'$R^2 = {r2_piston:.3f}$',
    transform=plt.gca().transAxes
)


# Tail dragger
plt.subplot(223)
plt.scatter(tail_dragger_mtow, tail_dragger_oem)

a_tail_dragger, b_tail_dragger = np.polyfit(tail_dragger_mtow, tail_dragger_oem, 1)
x_vals = np.linspace(np.min(tail_dragger_mtow), np.max(tail_dragger_mtow), 100)
y_pred = a_tail_dragger * tail_dragger_mtow + b_tail_dragger

ss_res = np.sum((tail_dragger_oem - y_pred) ** 2)
ss_tot = np.sum((tail_dragger_oem - np.mean(tail_dragger_oem)) ** 2)
r2_tail_dragger = 1 - (ss_res / ss_tot)

plt.plot(x_vals, a_tail_dragger * x_vals + b_tail_dragger)
plt.title(r'$m_{oe}$ vs $m_{to}$ for tail dragger aircraft')
plt.xlabel(r'$m_{to}$ [kg]')
plt.ylabel(r'$m_{oe}$ [kg]')
plt.text(
    0.05, 0.90,
    rf'$m_{{oe}} = {a_tail_dragger:.3f}m_{{to}} + {b_tail_dragger:.1f}$' + '\n' + rf'$R^2 = {r2_tail_dragger:.3f}$',
    transform=plt.gca().transAxes
)

print('Memory used is: '+ str(tracemalloc.get_traced_memory()))

tracemalloc.stop()
# Tricycle
plt.subplot(224)
plt.scatter(tricycle_mtow, tricycle_oem)

a_tricycle, b_tricycle = np.polyfit(tricycle_mtow, tricycle_oem, 1)
x_vals = np.linspace(np.min(tricycle_mtow), np.max(tricycle_mtow), 100)
y_pred = a_tricycle * tricycle_mtow + b_tricycle

ss_res = np.sum((tricycle_oem - y_pred) ** 2)
ss_tot = np.sum((tricycle_oem - np.mean(tricycle_oem)) ** 2)
r2_tricycle = 1 - (ss_res / ss_tot)

plt.plot(x_vals, a_tricycle * x_vals + b_tricycle)
plt.title(r'$m_{oe}$ vs $m_{to}$ for tricycle aircraft')
plt.xlabel(r'$m_{to}$ [kg]')
plt.ylabel(r'$m_{oe}$ [kg]')
plt.text(
    0.05, 0.90,
    rf'$m_{{oe}} = {a_tricycle:.3f}m_{{to}} + {b_tricycle:.1f}$' + '\n' + rf'$R^2 = {r2_tricycle:.3f}$',
    transform=plt.gca().transAxes
)


print("Turboprop regression:")
print(f"a = {a_turbo_prop:.6f}, b = {b_turbo_prop:.6f}, R^2 = {r2_turbo_prop:.6f}")
print(f"m_oe = {a_turbo_prop:.6f} m_to + {b_turbo_prop:.6f}")
print()

print("Piston regression:")
print(f"a = {a_piston:.6f}, b = {b_piston:.6f}, R^2 = {r2_piston:.6f}")
print(f"m_oe = {a_piston:.6f} m_to + {b_piston:.6f}")
print()

print("Tail dragger regression:")
print(f"a = {a_tail_dragger:.6f}, b = {b_tail_dragger:.6f}, R^2 = {r2_tail_dragger:.6f}")
print(f"m_oe = {a_tail_dragger:.6f} m_to + {b_tail_dragger:.6f}")
print()

print("Tricycle regression:")
print(f"a = {a_tricycle:.6f}, b = {b_tricycle:.6f}, R^2 = {r2_tricycle:.6f}")
print(f"m_oe = {a_tricycle:.6f} m_to + {b_tricycle:.6f}")
print()

plt.tight_layout()
plt.savefig('outputs/regression.svg')
plt.show()
