import matplotlib.pyplot as plt

# ============================================
# DATA
# ============================================

MTOW = 1870.0  # kg

concepts = {

    "Taildragger Configuration": {
        "x": [2.90, 3.23, 3.19, 2.88, 2.90],
        "m": [1046.0, 1750.0, 1870.0, 1166.0, 1046.0],
        "labels": ["OEW", "OEW+WP", "OEW+WP+WF", "OEW+WF", ""]
    },

    "Tricycle Configuration": {
        "x": [2.85, 3.20, 3.16, 2.83, 2.85],
        "m": [1046.0, 1750.0, 1870.0, 1166.0, 1046.0],
        "labels": ["OEW", "OEW+WP", "OEW+WP+WF", "OEW+WF", ""]
    }
}

# ============================================
# PLOT
# ============================================

plt.figure(figsize=(11,7))

markers = ['o', 's']

for i, (name, data) in enumerate(concepts.items()):

    xcg = data["x"]
    mass_fraction = [m / MTOW for m in data["m"]]
    labels = data["labels"]

    plt.plot(
        xcg,
        mass_fraction,
        marker=markers[i],
        linewidth=2.5,
        markersize=8,
        label=name
    )

    # ============================================
    # LABEL POSITIONS
    # ============================================

    if name == "Taildragger Configuration":

        offsets = [
            (25, -10),   # OEW
            (38, -5),    # OEW+WP
            (15, -30),   # OEW+WP+WF
            (-70, -15),  # OEW+WF
            (0, 0)
        ]

    else:  # Tricycle

        offsets = [
            (-60, 15),    # OEW
            (55, 20),     # OEW+WP
            (-120, 20),   # OEW+WP+WF
            (-85, 5),     # OEW+WF
            (0, 0)
        ]

    # ============================================
    # ANNOTATIONS WITH ARROWS
    # ============================================

    for x, y, label, offset in zip(xcg, mass_fraction, labels, offsets):

        if label != "":

            dx, dy = offset

            plt.annotate(
                label,
                xy=(x, y),
                xytext=(dx, dy),
                textcoords='offset points',
                fontsize=10,
                arrowprops=dict(
                    arrowstyle='->',
                    lw=1
                )
            )

    # ============================================
    # MOST AFT CG
    # ============================================

    most_aft_cg = max(xcg)

    print(f"{name}:")
    print(f"  Most aft CG = {most_aft_cg:.2f} m\n")

# ============================================
# FIGURE SETTINGS
# ============================================

plt.xlabel(r'$x_{cg}$ [m]', fontsize=16)
plt.ylabel(r'Mass fraction, $M/M_{TO}$ [-]', fontsize=16)

plt.xlim(2.75, 3.30)
plt.ylim(0.5, 1.05)

plt.grid(True, alpha=0.4)

plt.legend(fontsize=12)

plt.title('Class I Loading Diagram', fontsize=20)

plt.tight_layout()

plt.show()