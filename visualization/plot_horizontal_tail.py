import matplotlib.pyplot as plt
import numpy as np

# ============================================================
# EMPENNAGE PARAMETERS (Single Configuration)
# ============================================================

config_name = "Taildragger"
config = {
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
}

# ============================================================
# HORIZONTAL STABILIZER GEOMETRY
# (Quarter-chord sweep specified)
# ============================================================

def horizontal_geometry(S, AR, taper, quarter_chord_sweep_deg):
    b = np.sqrt(S * AR)
    c_root = 2 * S / (b * (1 + taper))
    c_tip = taper * c_root
    sweep_rad = np.radians(quarter_chord_sweep_deg)
    
    # LE offset required for zero quarter-chord sweep
    dx = ((b / 2) * np.tan(sweep_rad) + (c_root - c_tip) / 4)
    
    MAC = (2 / 3 * c_root * ((1 + taper + taper**2) / (1 + taper)))
    y_mac = (b / 6 * ((1 + 2 * taper) / (1 + taper)))
    x_mac_le = (y_mac * np.tan(sweep_rad) + (c_root - MAC) / 4)
    
    return {
        "b": b, "c_root": c_root, "c_tip": c_tip, 
        "dx": dx, "MAC": MAC, "y_mac": y_mac, "x_mac_le": x_mac_le
    }


# ============================================================
# VERTICAL STABILIZER GEOMETRY
# (Leading-edge sweep specified)
# ============================================================

def vertical_geometry(S, AR, taper, LE_sweep_deg):
    b = np.sqrt(S * AR)
    c_root = 2 * S / (b * (1 + taper))
    c_tip = taper * c_root
    sweep_rad = np.radians(LE_sweep_deg)
    
    # LE offset
    dx = b * np.tan(sweep_rad)
    
    MAC = (2 / 3 * c_root * ((1 + taper + taper**2) / (1 + taper)))
    z_mac = (b / 3 * ((1 + 2 * taper) / (1 + taper)))
    x_mac_le = z_mac * np.tan(sweep_rad)
    
    return {
        "b": b, "c_root": c_root, "c_tip": c_tip, 
        "dx": dx, "MAC": MAC, "z_mac": z_mac, "x_mac_le": x_mac_le
    }


# ============================================================
# CREATE FIGURE (1x2 Grid)
# ============================================================

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# ========================================================
# HORIZONTAL STABILIZER
# ========================================================

h = config["horizontal"]
g = horizontal_geometry(h["S"], h["AR"], h["taper"], h["quarter_chord_sweep_deg"])

b = g["b"]
c_root = g["c_root"]
c_tip = g["c_tip"]
dx = g["dx"]
MAC = g["MAC"]
y_mac = g["y_mac"]
x_mac_le = g["x_mac_le"]

ax1 = axes[0]

# Outline
x_outline = [0, c_root, dx + c_tip, dx, 0]
y_outline = [0, 0, b / 2, b / 2, 0]
y_outline_mirror = [0, 0, -b / 2, -b / 2, 0]

# Fill
ax1.fill(x_outline, y_outline, alpha=0.3)
ax1.fill(x_outline, y_outline_mirror, alpha=0.3)

# Outline
ax1.plot(x_outline, y_outline, 'k')
ax1.plot(x_outline, y_outline_mirror, 'k')

# Chords
ax1.plot([0, c_root], [0, 0], linewidth=3, label="Root chord")
ax1.plot([dx, dx + c_tip], [b / 2, b / 2], linewidth=3, label="Tip chord")

# MAC
ax1.plot([x_mac_le, x_mac_le + MAC], [y_mac, y_mac], 'r', linewidth=4, label="MAC")

# LE marker
ax1.scatter(x_mac_le, y_mac, color='blue', zorder=5, label="LE of MAC")
ax1.text(x_mac_le, y_mac + 0.08, "LE", color='blue', fontsize=10, fontweight='bold')

# Quarter chord point
c4_x = x_mac_le + 0.25 * MAC
ax1.scatter(c4_x, y_mac, color='red', zorder=5, label="c/4 Point")
ax1.text(c4_x, y_mac - 0.15, "c/4", color='red', fontsize=10, fontweight='bold')

# Centerline
ax1.axhline(0, linestyle='--', color='gray')

ax1.set_title(f"{config_name} Horizontal Stabilizer")
ax1.set_xlabel("x [m]")
ax1.set_ylabel("y [m]")
ax1.axis('equal')
ax1.grid(True)
ax1.legend(loc="upper right")

# ========================================================
# VERTICAL STABILIZER
# ========================================================

v = config["vertical"]
g_v = vertical_geometry(v["S"], v["AR"], v["taper"], v["LE_sweep_deg"])

b_v = g_v["b"]
c_root_v = g_v["c_root"]
c_tip_v = g_v["c_tip"]
dx_v = g_v["dx"]
MAC_v = g_v["MAC"]
z_mac = g_v["z_mac"]
x_mac_le_v = g_v["x_mac_le"]

ax2 = axes[1]

# Outline
xv = [0, c_root_v, dx_v + c_tip_v, dx_v, 0]
zv = [0, 0, b_v, b_v, 0]

# Fill
ax2.fill(xv, zv, alpha=0.3, color='orange')

# Outline
ax2.plot(xv, zv, 'k')

# Chords
ax2.plot([0, c_root_v], [0, 0], linewidth=3, label="Root chord")
ax2.plot([dx_v, dx_v + c_tip_v], [b_v, b_v], linewidth=3, label="Tip chord")

# MAC
ax2.plot([x_mac_le_v, x_mac_le_v + MAC_v], [z_mac, z_mac], 'r', linewidth=4, label="MAC")

# LE marker
ax2.scatter(x_mac_le_v, z_mac, color='blue', zorder=5, label="LE of MAC")
ax2.text(x_mac_le_v, z_mac + 0.08, "LE", color='blue', fontsize=10, fontweight='bold')

# Quarter chord point
c4_x_v = x_mac_le_v + 0.25 * MAC_v
ax2.scatter(c4_x_v, z_mac, color='red', zorder=5, label="c/4 Point")
ax2.text(c4_x_v, z_mac - 0.15, "c/4", color='red', fontsize=10, fontweight='bold')

ax2.set_title(f"{config_name} Vertical Stabilizer")
ax2.set_xlabel("x [m]")
ax2.set_ylabel("z [m]")
ax2.axis('equal')
ax2.grid(True)
ax2.legend(loc="upper right")

# ============================================================
# FINAL LAYOUT & SAVE
# ============================================================

plt.tight_layout(pad=2.0)

# Save the figure as a PDF
plt.savefig(f"{config_name}_Empennage.pdf", format='pdf', bbox_inches='tight')

plt.show()