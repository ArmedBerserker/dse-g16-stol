import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, interp1d


def get_xflr5_loads(filepath, rho=1.225, wingbox_x_pc=0.25, engines=None, struct_weight=None):

    # Wing load calculator handling aero, distributed wing weight,
    # and 3D engine point loads (Thrust/Mass with offsets).

    # --- 1. DATA EXTRACTION ---
    with open(filepath, 'r') as f:
        lines = f.readlines()

    v_inf = next(float(line.split('=')[1].split('m/s')[0]) for line in lines if "QInf" in line)
    q_dyn = 0.5 * rho * v_inf ** 2
    strips_raw = "".join(lines).split("Strip")[1:]

    y_raw, fz_raw, fx_raw, ty_raw = [], [], [], []

    for block in strips_raw:
        data = np.array([line.split() for line in block.strip().split('\n') if len(line.split()) == 9], dtype=float)
        if not data.size: continue

        y_avg = np.mean(data[:, 2])
        if y_avg < 0: continue

        area, cp = data[:, 7], data[:, 8]
        normals, coords = data[:, 4:7], data[:, 1:4]

        # Aerodynamic Force (N)
        f_vecs = (-cp * q_dyn * area)[:, None] * normals

        # Moment relative to Elastic Axis
        x_le, x_te = np.min(coords[:, 0]), np.max(coords[:, 0])
        ea_x = x_le + wingbox_x_pc * (x_te - x_le)
        arms = coords.copy()
        arms[:, 0] -= ea_x
        m_vecs = np.cross(arms, f_vecs)

        y_raw.append(y_avg)
        fz_raw.append(np.sum(f_vecs[:, 2]))  # Lift
        fx_raw.append(np.sum(f_vecs[:, 0]))  # Drag
        ty_raw.append(np.sum(m_vecs[:, 1]))  # Aero Torsion

    # Establish Span
    sort_idx = np.argsort(y_raw)
    y_c = np.array(y_raw)[sort_idx]
    y_tip = y_c[-1] + (y_c[-1] - y_c[-2]) / 2 if len(y_c) > 1 else y_c[0] * 1.05

    # --- 2. CONTINUOUS DISTRIBUTIONS (N/m or Nm/m) ---
    # noinspection PyTypeChecker
    def get_spline_dist(y_centers, forces, tip_val=0.0):
        midpoints = 0.5 * (y_centers[:-1] + y_centers[1:])
        boundaries = np.concatenate(([0], midpoints, [y_tip]))
        widths = np.diff(boundaries)
        w_vals = np.array(forces) / widths
        y_fit = np.concatenate(([0], y_centers, [y_tip]))
        w_fit = np.concatenate(([w_vals[0]], w_vals, [tip_val]))
        return CubicSpline(y_fit, w_fit, bc_type=(((1, 0.0)), 'natural'))

    w_fz_aero = get_spline_dist(y_c, fz_raw)
    w_fx_aero = get_spline_dist(y_c, fx_raw)
    w_ty_aero = get_spline_dist(y_c, ty_raw)

    # Distributed Structural Weight (Linear Interpolation)
    if struct_weight:
        # Interpolate provided mass/weight grid to the wing span
        w_struct_func = interp1d(struct_weight['y'], struct_weight['w'],
                                 kind='linear', fill_value="extrapolate")
    else:
        w_struct_func = lambda y: 0

    # --- 3. NUMERICAL INTEGRATION (Tip-to-Root) ---
    y_fine = np.linspace(0, y_tip, 1000)
    Vz, Vx, Mx, Mz, Ty = [np.zeros_like(y_fine) for _ in range(5)]

    for i in range(len(y_fine) - 2, -1, -1):
        dy = y_fine[i + 1] - y_fine[i]
        ym = 0.5 * (y_fine[i] + y_fine[i + 1])

        # Distributed Force Summation (Aero + Structural)
        # Note: Structural weight is negative (acting down)
        w_z = w_fz_aero(ym) + w_struct_func(ym)
        w_x = w_fx_aero(ym)
        w_t = w_ty_aero(ym)

        # Shear Integration (Trapezoidal)
        Vz[i] = Vz[i + 1] + w_z * dy
        Vx[i] = Vx[i + 1] + w_x * dy
        Ty[i] = Ty[i + 1] + w_t * dy

        # Bending Moment Integration
        # M = integral of Shear
        Mx[i] = Mx[i + 1] + 0.5 * (Vz[i] + Vz[i + 1]) * dy
        Mz[i] = Mz[i + 1] - 0.5 * (Vx[i] + Vx[i + 1]) * dy

    # --- 4. SUPERIMPOSE ENGINE POINT LOADS & MOMENTS ---
    if engines:
        for eng in engines:
            mask = y_fine <= eng['y']

            f_weight = -eng['mass'] * 9.81
            f_thrust = eng['thrust']

            # Distance from engine station to current station
            lever_y = eng['y'] - y_fine[mask]

            # A. Point Forces -> Jumps in Shear, Kinks in Bending
            Vz[mask] += f_weight
            Mx[mask] += f_weight * lever_y

            Vx[mask] += f_thrust
            Mz[mask] -= f_thrust * lever_y  # Thrust produces -Mz bending

            # B. Point Moments -> Jumps in Bending/Torsion
            # Torsion jump from X-offset (Weight) and Z-offset (Thrust)
            # Ty = Fz * x_offset + Fx * z_offset
            t_jump = (f_weight * eng.get('x_offset', 0)) + (f_thrust * eng.get('z_offset', 0))
            Ty[mask] += t_jump

            # Optional: If the pylon produces a pitching moment (My)
            # that doesn't exist in this 1D Ty-only torsion model,
            # you could add jumps to Mx or Mz here if needed.

    return {"y": y_fine, "Vz": Vz, "Vx": Vx, "Mx": Mx, "Mz": Mz, "Ty": Ty}


# --- 5. EXECUTION ---
# x_offset: pos = forward of EA | z_offset: pos = above EA
engine_list = [{
    'y': 2.5,
    'mass': 700,
    'thrust': 5000,
    'x_offset': 0.8,  # Engine is forward
    'z_offset': -0.4  # Engine is underslung
}]

# Wing weight in N/m (Trapezoidal distribution example)
wingweight = {'y': [0, 4, 8], 'w': [-800, -400, -100]}

res = get_xflr5_loads("MainWing_a=5.00_v=75.00ms.txt", engines=None, struct_weight=None)

# --- PLOTTING ---
fig, axes = plt.subplots(5, 1, figsize=(10, 15), sharex=True)
plots = [
    ('Vz', 'crimson', 'Vert. Shear Vz [N]'),
    ('Mx', 'forestgreen', 'Bending Mx [Nm]'),
    ('Vx', 'darkorange', 'Lat. Shear Vx [N]'),
    ('Mz', 'darkmagenta', 'Bending Mz [Nm]'),
    ('Ty', 'dodgerblue', 'Torsion Ty [Nm]')
]

for ax, (key, color, label) in zip(axes, plots):
    ax.plot(res['y'], res[key], color=color, lw=2)
    ax.set_ylabel(label, fontsize=10, fontweight='bold')
    ax.grid(True, alpha=0.3)

print(res['Vz'][0])
print(res['Mx'][0])
print(res['Vx'][0])
print(res['Mz'][0])
print(res['Ty'][0])

axes[-1].set_xlabel('Spanwise Position y [m]')
plt.tight_layout()
plt.show()