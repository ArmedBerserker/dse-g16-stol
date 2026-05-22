import matplotlib.pyplot as plt
import numpy as np

# ============================================================
# EMPENNAGE PARAMETERS
# ============================================================

configs = {
    "Taildragger": {
        "horizontal": {
            "S": 7.13,
            "AR": 4.5,
            "taper": 0.7,
            "quarter_chord_sweep_deg": 0
        },
        "vertical": {
            "S": 5.12,
            "AR": 1.4,
            "taper": 0.5,
            "LE_sweep_deg": 30
        }
    },

    "Tricycle": {
        "horizontal": {
            "S": 7.09,
            "AR": 4.5,
            "taper": 0.7,
            "quarter_chord_sweep_deg": 0
        },
        "vertical": {
            "S": 5.09,
            "AR": 1.4,
            "taper": 0.5,
            "LE_sweep_deg": 30
        }
    }
}


# ============================================================
# HORIZONTAL STABILIZER GEOMETRY
# (Quarter-chord sweep specified)
# ============================================================

def horizontal_geometry(S, AR, taper, quarter_chord_sweep_deg):

    # Span
    b = np.sqrt(S * AR)

    # Root chord
    c_root = 2 * S / (b * (1 + taper))

    # Tip chord
    c_tip = taper * c_root

    # Quarter-chord sweep
    sweep_rad = np.radians(quarter_chord_sweep_deg)

    # LE offset required for zero quarter-chord sweep
    dx = (
        (b / 2) * np.tan(sweep_rad)
        + (c_root - c_tip) / 4
    )

    # MAC
    MAC = (
        2 / 3
        * c_root
        * ((1 + taper + taper**2) / (1 + taper))
    )

    # MAC spanwise location
    y_mac = (
        b / 6
        * ((1 + 2 * taper) / (1 + taper))
    )

    # LE of MAC
    x_mac_le = (
        y_mac * np.tan(sweep_rad)
        + (c_root - MAC) / 4
    )

    return {
        "b": b,
        "c_root": c_root,
        "c_tip": c_tip,
        "dx": dx,
        "MAC": MAC,
        "y_mac": y_mac,
        "x_mac_le": x_mac_le
    }


# ============================================================
# VERTICAL STABILIZER GEOMETRY
# (Leading-edge sweep specified)
# ============================================================

def vertical_geometry(S, AR, taper, LE_sweep_deg):

    # Span
    b = np.sqrt(S * AR)

    # Root chord
    c_root = 2 * S / (b * (1 + taper))

    # Tip chord
    c_tip = taper * c_root

    # LE sweep
    sweep_rad = np.radians(LE_sweep_deg)

    # LE offset
    dx = b * np.tan(sweep_rad)

    # MAC
    MAC = (
        2 / 3
        * c_root
        * ((1 + taper + taper**2) / (1 + taper))
    )

    # MAC vertical location
    z_mac = (
        b / 3
        * ((1 + 2 * taper) / (1 + taper))
    )

    # LE of MAC
    x_mac_le = z_mac * np.tan(sweep_rad)

    return {
        "b": b,
        "c_root": c_root,
        "c_tip": c_tip,
        "dx": dx,
        "MAC": MAC,
        "z_mac": z_mac,
        "x_mac_le": x_mac_le
    }


# ============================================================
# CREATE FIGURE
# ============================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 12))


# ============================================================
# LOOP THROUGH CONFIGURATIONS
# ============================================================

for col, (config_name, config) in enumerate(configs.items()):

    # ========================================================
    # HORIZONTAL STABILIZER
    # ========================================================

    h = config["horizontal"]

    g = horizontal_geometry(
        h["S"],
        h["AR"],
        h["taper"],
        h["quarter_chord_sweep_deg"]
    )

    b = g["b"]
    c_root = g["c_root"]
    c_tip = g["c_tip"]
    dx = g["dx"]

    MAC = g["MAC"]
    y_mac = g["y_mac"]
    x_mac_le = g["x_mac_le"]

    ax = axes[0, col]

    # Outline
    x_outline = [
        0,
        c_root,
        dx + c_tip,
        dx,
        0
    ]

    y_outline = [
        0,
        0,
        b / 2,
        b / 2,
        0
    ]

    y_outline_mirror = [
        0,
        0,
        -b / 2,
        -b / 2,
        0
    ]

    # Fill
    ax.fill(x_outline, y_outline, alpha=0.3)
    ax.fill(x_outline, y_outline_mirror, alpha=0.3)

    # Outline
    ax.plot(x_outline, y_outline, 'k')
    ax.plot(x_outline, y_outline_mirror, 'k')

    # Root chord
    ax.plot(
        [0, c_root],
        [0, 0],
        linewidth=3,
        label="Root chord"
    )

    # Tip chord
    ax.plot(
        [dx, dx + c_tip],
        [b / 2, b / 2],
        linewidth=3,
        label="Tip chord"
    )

    # MAC
    ax.plot(
        [x_mac_le, x_mac_le + MAC],
        [y_mac, y_mac],
        'r',
        linewidth=4,
        label="MAC"
    )

    # LE marker
    ax.scatter(
        x_mac_le,
        y_mac,
        color='blue',
        zorder=5,
        label="LE of MAC"
    )

    ax.text(
        x_mac_le,
        y_mac + 0.08,
        "LE",
        color='blue',
        fontsize=10,
        fontweight='bold'
    )

    # Quarter chord point
    ax.scatter(
        x_mac_le + 0.25 * MAC,
        y_mac,
        color='red',
        zorder=5
    )

    # Centerline
    ax.axhline(
        0,
        linestyle='--',
        color='gray'
    )

    ax.set_title(
        f"{config_name} Horizontal Stabilizer"
    )

    ax.set_ylabel("y [m]")

    ax.axis('equal')
    ax.grid(True)
    ax.legend()

    # ========================================================
    # VERTICAL STABILIZER
    # ========================================================

    v = config["vertical"]

    g_v = vertical_geometry(
        v["S"],
        v["AR"],
        v["taper"],
        v["LE_sweep_deg"]
    )

    b_v = g_v["b"]
    c_root_v = g_v["c_root"]
    c_tip_v = g_v["c_tip"]
    dx_v = g_v["dx"]

    MAC_v = g_v["MAC"]
    z_mac = g_v["z_mac"]
    x_mac_le_v = g_v["x_mac_le"]

    ax2 = axes[1, col]

    # Outline
    xv = [
        0,
        c_root_v,
        dx_v + c_tip_v,
        dx_v,
        0
    ]

    zv = [
        0,
        0,
        b_v,
        b_v,
        0
    ]

    # Fill
    ax2.fill(
        xv,
        zv,
        alpha=0.3,
        color='orange'
    )

    # Outline
    ax2.plot(xv, zv, 'k')

    # Root chord
    ax2.plot(
        [0, c_root_v],
        [0, 0],
        linewidth=3,
        label="Root chord"
    )

    # Tip chord
    ax2.plot(
        [dx_v, dx_v + c_tip_v],
        [b_v, b_v],
        linewidth=3,
        label="Tip chord"
    )

    # MAC
    ax2.plot(
        [x_mac_le_v, x_mac_le_v + MAC_v],
        [z_mac, z_mac],
        'r',
        linewidth=4,
        label="MAC"
    )

    # LE marker
    ax2.scatter(
        x_mac_le_v,
        z_mac,
        color='blue',
        zorder=5,
        label="LE of MAC"
    )

    ax2.text(
        x_mac_le_v,
        z_mac + 0.08,
        "LE",
        color='blue',
        fontsize=10,
        fontweight='bold'
    )

    # Quarter chord point
    ax2.scatter(
        x_mac_le_v + 0.25 * MAC_v,
        z_mac,
        color='red',
        zorder=5
    )

    ax2.set_title(
        f"{config_name} Vertical Stabilizer"
    )

    ax2.set_xlabel("x [m]")
    ax2.set_ylabel("z [m]")

    ax2.axis('equal')
    ax2.grid(True)
    ax2.legend()


# ============================================================
# FINAL LAYOUT
# ============================================================

plt.tight_layout()
plt.show()