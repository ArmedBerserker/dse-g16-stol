import matplotlib.pyplot as plt

# ============================================
# DATA BASED ON LATEX
# ============================================

MTOW = 1870.0  # kg

concepts = {
    "Concept 1 (Piston TD)": {
        "x": [3.30, 3.24, 3.22, 3.26, 3.30],
        "m": [1046.2, 1750.2, 1870.0, 1166.0, 1046.2]
    },

    "Concept 2 (TP TD)": {
        "x": [3.31, 3.25, 3.23, 3.27, 3.31],
        "m": [1046.2, 1750.2, 1870.0, 1166.0, 1046.2]
    },

    "Concept 3 (Piston Tri)": {
        "x": [3.24, 3.20, 3.19, 3.21, 3.24],
        "m": [1046.2, 1750.2, 1870.0, 1166.0, 1046.2]
    },

    "Concept 4 (TP Tri)": {
        "x": [3.25, 3.21, 3.19, 3.22, 3.25],
        "m": [1046.2, 1750.2, 1870.0, 1166.0, 1046.2]
    }
}

# ============================================
# PLOT
# ============================================

plt.figure(figsize=(10,6))

markers = ['o', 's', '^', 'D']

for i, (name, data) in enumerate(concepts.items()):

    xcg = data["x"]
    mass_fraction = [m / MTOW for m in data["m"]]

    plt.plot(
        xcg,
        mass_fraction,
        marker=markers[i],
        linewidth=2,
        markersize=7,
        label=name
    )

# ============================================
# FIGURE SETTINGS
# ============================================

plt.xlabel(r'$x_{cg}$ [m]', fontsize=12)
plt.ylabel(r'Mass fraction, $M/M_{TO}$ [-]', fontsize=12)

plt.xlim(3.1, 3.35)
plt.ylim(0.5, 1.05)

plt.grid(True)
plt.legend()
plt.tight_layout()

plt.title('Class I Loading Diagram')

plt.show()