import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button


def visualize_comparison(csv_path):
    df = pd.read_csv(csv_path)

    fig = plt.figure(figsize=(12, 8))
    plt.subplots_adjust(bottom=0.2)
    ax = fig.add_subplot(111, projection='3d')

    # Scaling relative to the largest force to keep arrows visible
    force_scale = 1.0 / df['Fz'].max()

    def draw_wing():
        for i, row in df.iterrows():
            # 1. Geometry Setup
            le = np.array([row['p_le_x'], row['p_le_y'], row['p_le_z']])
            chord_u = np.array([row['cx'], row['cy'], row['cz']])
            norm_u = np.array([row['nx'], row['ny'], row['nz']])
            chord_len = row['chord']
            p_c4 = le + (0.25 * chord_len * chord_u)

            # Total force vector from CSV
            f_total = np.array([row['Fx'], row['Fy'], row['Fz']])

            # --- LOCAL COMPONENTS (BLUE) ---
            # Normal force
            fn_mag = np.dot(f_total, norm_u)
            fn_vec = fn_mag * norm_u

            # Tangential force
            ft_mag = np.dot(f_total, chord_u)
            ft_vec = ft_mag * chord_u

            # --- GLOBAL COMPONENTS (RED) ---
            # Global Lift (Purely Vertical Z)
            fz_global = np.array([0, 0, row['Fz']])

            # Global Drag (Purely Longitudinal X)
            fx_global = np.array([row['Fx'], 0, 0])

            # Plotting
            ax.plot([le[0], le[0] + chord_len * chord_u[0]],
                    [le[1], le[1] + chord_len * chord_u[1]],
                    [le[2], le[2] + chord_len * chord_u[2]], color='gray', alpha=0.3)

            # Local Blue Vectors
            ax.quiver(p_c4[0], p_c4[1], p_c4[2], fn_vec[0] * force_scale, fn_vec[1] * force_scale,
                      fn_vec[2] * force_scale,
                      color='blue', linewidth=2, label='Local Normal' if i == 0 else "")
            ax.quiver(p_c4[0], p_c4[1], p_c4[2], ft_vec[0] * force_scale, ft_vec[1] * force_scale,
                      ft_vec[2] * force_scale,
                      color='deepskyblue', linewidth=2, label='Local Tangential' if i == 0 else "")

            # Global Red Vectors
            ax.quiver(p_c4[0], p_c4[1], p_c4[2], fz_global[0] * force_scale, fz_global[1] * force_scale,
                      fz_global[2] * force_scale,
                      color='red', linewidth=1.5, linestyle='--', label='Global Lift (Z)' if i == 0 else "")
            ax.quiver(p_c4[0], p_c4[1], p_c4[2], fx_global[0] * force_scale, fx_global[1] * force_scale,
                      fx_global[2] * force_scale,
                      color='crimson', linewidth=1.5, linestyle='--', label='Global Drag (X)' if i == 0 else "")

    draw_wing()

    # View Controls
    max_range = df['y'].max() / 2.0
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(df['y'].mean() - max_range, df['y'].mean() + max_range)
    ax.set_zlim(-max_range, max_range)
    ax.set_xlabel('Global X');
    ax.set_ylabel('Global Y');
    ax.set_zlabel('Global Z')
    ax.legend(loc='upper left', fontsize='small')

    class Index:
        def front(self, event): ax.view_init(elev=0, azim=0); plt.draw()

        def side(self, event): ax.view_init(elev=0, azim=-90); plt.draw()

        def top(self, event): ax.view_init(elev=90, azim=-90); plt.draw()

        def iso(self, event): ax.view_init(elev=30, azim=-45); plt.draw()

    callback = Index()
    ax_front = plt.axes([0.1, 0.05, 0.18, 0.06]);
    ax_side = plt.axes([0.3, 0.05, 0.18, 0.06])
    ax_top = plt.axes([0.5, 0.05, 0.18, 0.06]);
    ax_iso = plt.axes([0.7, 0.05, 0.18, 0.06])
    btn_front = Button(ax_front, 'Front (YZ)');
    btn_side = Button(ax_side, 'Side (XZ)')
    btn_top = Button(ax_top, 'Top (XY)');
    btn_iso = Button(ax_iso, 'Isometric')
    btn_front.on_clicked(callback.front);
    btn_side.on_clicked(callback.side)
    btn_top.on_clicked(callback.top);
    btn_iso.on_clicked(callback.iso)

    plt.show()


if __name__ == "__main__":
    visualize_comparison('xflr5_parsed.csv')