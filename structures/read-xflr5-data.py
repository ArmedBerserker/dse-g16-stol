import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid


def parse_xflr5(filepath, rho=1.225, wingbox_x_pc=0.40):
    q_dyn = 0.0
    strips = []
    current_strip = None

    with open(filepath, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if "QInf" in line:
            v_inf = float(line.split('=')[1].split('m/s')[0].strip())
            q_dyn = 0.5 * rho * v_inf ** 2
            break

    parsing_cp = False
    for line in lines:
        clean = line.strip()
        if "Main Wing Cp Coefficients" in line:
            parsing_cp = True
            continue
        if not parsing_cp: continue

        if clean.startswith("Strip"):
            if current_strip: strips.append(current_strip)
            current_strip = {k: [] for k in ['x', 'y', 'z', 'nx', 'ny', 'nz', 'area', 'cp']}
        elif current_strip is not None:
            parts = clean.split()
            if len(parts) == 9:
                try:
                    for i, k in enumerate(['x', 'y', 'z', 'nx', 'ny', 'nz', 'area', 'cp'], 1):
                        current_strip[k].append(float(parts[i]))
                except ValueError:
                    continue
    if current_strip: strips.append(current_strip)

    y_stations = []
    # Force/Moment densities per unit span
    fz_p, fx_p, ty_p = [], [], []

    for s in strips:
        cp, area = np.array(s['cp']), np.array(s['area'])
        pts = np.vstack((s['x'], s['y'], s['z']))
        n = np.vstack((s['nx'], s['ny'], s['nz']))

        # 3D Force per panel (N)
        f_vecs = -cp * q_dyn * area * n

        # Calculate local Elastic Axis X-position
        x_le, x_te = np.min(s['x']), np.max(s['x'])
        ea_x = x_le + (wingbox_x_pc * (x_te - x_le))

        # 3D Moment about local EA [Mx, My, Mz]
        # Arm relative to EA (x) and Root (y=0, z=0)
        arms = np.vstack((pts[0, :] - ea_x, pts[1, :], pts[2, :]))
        m_vecs = np.cross(arms, f_vecs, axis=0)

        y_stations.append(np.mean(s['y']))
        fz_p.append(np.sum(f_vecs[2, :]))  # Vertical force
        fx_p.append(np.sum(f_vecs[0, :]))  # Drag/Backwards force
        ty_p.append(np.sum(m_vecs[1, :]))  # Torsion around Y-axis

    # Vectorization and sorting
    y_f = np.array(y_stations)
    idx = np.argsort(y_f)
    mask = y_f[idx] >= -1e-3

    y = y_f[idx][mask]
    fz = np.array(fz_p)[idx][mask]
    fx = np.array(fx_p)[idx][mask]
    ty = np.array(ty_p)[idx][mask]

    # Root interpolation (y=0) from nearest panel
    if y[0] > 1e-4:
        y = np.insert(y, 0, 0.0)
        fz = np.insert(fz, 0, fz[0])
        fx = np.insert(fx, 0, fx[0])
        ty = np.insert(ty, 0, ty[0])

    # Pre-allocate NVM arrays
    vz, vx = np.zeros_like(y), np.zeros_like(y)
    mx, mz, ty_int = np.zeros_like(y), np.zeros_like(y), np.zeros_like(y)

    for i in range(len(y)):
        y_ob = y[i:]
        if len(y_ob) > 1:
            # Shear Forces
            vz[i] = trapezoid(fz[i:], y_ob)
            vx[i] = trapezoid(fx[i:], y_ob)
            # Moments (Integral of Force * Arm)
            mx[i] = trapezoid(fz[i:] * (y_ob - y[i]), y_ob)
            mz[i] = trapezoid(fx[i:] * (y_ob - y[i]), y_ob)
            # Torsion (Direct integral of torsional density)
            ty_int[i] = trapezoid(ty[i:], y_ob)

    return {"y": y, "Vz": vz, "Vx": vx, "Mx": mx, "Mz": mz, "Ty": ty_int}


# Execution
res = parse_xflr5("MainWing_a=5.00_v=75.00ms.txt")

fig, ax = plt.subplots(5, 1, figsize=(10, 15), sharex=True)
ax[0].plot(res['y'], res['Vz'], 'r', label='Vertical Shear (Vz)')
ax[1].plot(res['y'], res['Mx'], 'g', label='Vertical Bending (Mx)')
ax[2].plot(res['y'], res['Vx'], 'orange', label='Backwards Shear (Vx)')
ax[3].plot(res['y'], res['Mz'], 'm', label='Backwards Bending (Mz)')
ax[4].plot(res['y'], res['Ty'], 'b', label='Torsion (Ty)')

for a in ax: a.grid(True); a.legend(loc='upper right')
ax[-1].set_xlabel('Spanwise Position y [m]')
plt.tight_layout()
plt.show()

print(f" Root Vertical Bending: {res['Mx'][0]:.2f} Nm")
print(f" Root Weak Bending: {res['Mz'][0]:.2f} Nm")
print(f" Root Torsion: {res['Ty'][0]:.2f} Nm")

# Exporting into separate NumPy arrays
y_coords = res['y'] # Spanwise positions
shear_v_z = res['Vz'] # Vertical Shear
shear_v_x = res['Vx'] # Backwards Shear (Drag-induced)
moment_m_x = res['Mx'] # Vertical Bending Moment
moment_m_z = res['Mz'] # Weak-axis Bending Moment
torsion_t_y = res['Ty'] # Torsional Moment

internal_loads_df = pd.DataFrame(
    {'y_m': y_coords, 'Vz_N': shear_v_z, 'Mx_Nm': moment_m_x, 'Vx_N': shear_v_x, 'Mz_Nm': moment_m_z, 'Ty_Nm': torsion_t_y})

# Save to CSV
internal_loads_df.to_csv('internal_loads.csv', index=False)




#import numpy as np

# Load the data, skipping the first row (the header)
#data = np.genfromtxt("Internal_Loads_Output.csv", delimiter=',', skip_header=1)

# Extract columns by index (0=y, 1=Vz, 2=Mx, etc.)
#y_coords = data[:, 0]
#bending_mx = data[:, 2]
