import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def generate_nvm_diagrams(csv_file):
    df = pd.read_csv(csv_file)

    # Sort Tip  to Root (Low Y)
    df = df.sort_values('y', ascending=False).reset_index(drop=True)

    n = len(df)
    # Result arrays for internal reactions
    vn, vt = np.zeros(n), np.zeros(n)
    mn, mt = np.zeros(n), np.zeros(n)
    ms = np.zeros(n)

    for i in range(n):
        y_station = df.loc[i, 'y']
        outboard = df.iloc[:i + 1]

        # INTERNAL SHEAR: Opposes the sum of external forces
        vn[i] = -outboard['Fn'].sum()
        vt[i] = -outboard['Ft'].sum()

        # INTERNAL BENDING MOMENTS: Oppose the moments from external forces
        lever_arms = outboard['y'] - y_station
        mt[i] = -(outboard['Fn'] * lever_arms).sum()
        mn[i] = -(outboard['Ft'] * lever_arms).sum()

        # INTERNAL TORSION: Opposes external pitching moments
        ms[i] = -outboard['Torsion'].sum()

    df['Vn'], df['Vt'] = vn, vt
    df['Mn'], df['Mt'] = mn, mt
    df['Ms_Torsion'] = ms

    # Plotting Logic
    fig, axes = plt.subplots(5, 1, figsize=(10, 16), sharex=True)
    plt.subplots_adjust(hspace=0.4)

    plot_data = [
        ('Vn', 'Normal (Vn)', 'blue'),
        ('Vt', 'Tangential (Vt)', 'green'),
        ('Mt', 'Bending (Mt)', 'red'),
        ('Mn', 'Bending (Mn)', 'orange'),
        ('Ms_Torsion', 'Torsion (Ms)', 'purple')
    ]

    for i, (col, label, color) in enumerate(plot_data):
        axes[i].plot(df['y'], df[col], color=color, lw=2, marker='.', markersize=4)
        axes[i].set_ylabel(f"{label}")
        axes[i].grid(True, alpha=0.3)
        axes[i].axhline(0, color='black', lw=0.8)

    axes[-1].set_xlabel("Spanwise Position y [m]")
    plt.show()

generate_nvm_diagrams('xflr5_parsed.csv')