import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# INPUT DATA
# ============================================================

L = 11.0
g = 9.81
safety_factor = 1.5

W_fuse = (38 + 32.1 + 43.1 + 5.5 + 49.9 + 152) * g * 3.8
W_v = 33.5 * g * 3.8
W_h = 20.0 * g * 3.8
F_row = 77 * 2 * g * 3.8
F_cargo = 200 * g * 3.8
F_h = 300 # TBD

x_nose_cone = 1.0
x1 = 3.15
x2 = 4.05
x3 = 5.20
x_cargo = 7.10
x_tail_cone = 7.50
x_h = 10.70
x_v = 9.924
x_w = 5.30
x_hf = 10.80

R_y = -(
    W_fuse + W_v + W_h +
    3 * F_row + F_cargo + F_h
)

M_z = -(
    - 0.1 * W_fuse * 0.5 * x_nose_cone
    - F_row * (x1 - x_nose_cone)
    - F_row * (x2 - x_nose_cone)
    - F_row * (x3 - x_nose_cone)
    - F_cargo * (x_cargo - x_nose_cone)
    - 0.65 * W_fuse * (x_tail_cone - x_nose_cone) / 2
    - W_h * (x_h - x_tail_cone)
    - W_v * (x_v - x_tail_cone)
    - 0.25 * W_fuse * 0.5 * (L - x_tail_cone)
    - F_h * (x_hf - x_tail_cone)
    - R_y * (x_w - x_nose_cone)
    - (0.25 * W_fuse + W_v + W_h + F_h) * (x_tail_cone - x_nose_cone)
)

# ============================================================
# LOADS
# ============================================================

loads = {
    "point_forces": [

        {"z": 0.0,
         "Fx": 0.1 * W_fuse},

        {"z": x1 - x_nose_cone,
         "Fx": F_row},

        {"z": x2 - x_nose_cone,
         "Fx": F_row},

        {"z": x3 - x_nose_cone,
         "Fx": F_row},

        {"z": x_cargo - x_nose_cone,
         "Fx": F_cargo},

        {"z": x_tail_cone - x_nose_cone,
         "Fx": 0.25 * W_fuse + W_v + W_h + F_h},

        {"z": x_w - x_nose_cone,
         "Fx": R_y},
    ],

    "distributed_loads": [

        {
            "z1": 0.0,
            "z2": x_tail_cone - x_nose_cone,
            "qx": 0.65 * W_fuse /
                  (x_tail_cone - x_nose_cone)
        }

    ],

    "moments": [

        {
            "z": 0.0,
            "My": -0.1 * W_fuse * 0.5 * x_nose_cone
        },

        {
            "z": x_w - x_nose_cone,
            "My": M_z
        },

        {
            "z": x_tail_cone - x_nose_cone,
            "My": -F_h * (x_hf - x_tail_cone)
        },

        {
            "z": x_tail_cone - x_nose_cone,
            "My": -0.25 * W_fuse * 0.5 * (L - x_tail_cone)
        },

        {
            "z": x_tail_cone - x_nose_cone,
            "My": -W_v * (x_v - x_tail_cone)
        },

        {
            "z": x_tail_cone - x_nose_cone,
            "My": -W_h * (x_h - x_tail_cone)
        }

    ]
}

# ============================================================
# DISCRETIZATION
# ============================================================

beam_length = x_tail_cone - x_nose_cone

z = np.linspace(0, beam_length, 10000)

V = np.zeros_like(z)
M = np.zeros_like(z)

# ============================================================
# POINT LOADS
# ============================================================

for p in loads["point_forces"]:

    zp = p["z"]
    Fx = p["Fx"]

    active = z >= zp

    V += np.where(active, Fx, 0.0)

    M += np.where(
        active,
        Fx * (z - zp),
        0.0
    )

# ============================================================
# DISTRIBUTED LOADS
# ============================================================

for q in loads["distributed_loads"]:

    z1 = q["z1"]
    z2 = q["z2"]
    w = q["qx"]

    inside = (z >= z1) & (z <= z2)
    after = z > z2

    # -------------------------
    # SHEAR CONTRIBUTION
    # -------------------------

    V += np.where(
        inside,
        w * (z - z1),
        0.0
    )

    V += np.where(
        after,
        w * (z2 - z1),
        0.0
    )

    # -------------------------
    # MOMENT CONTRIBUTION
    # -------------------------

    M += np.where(
        inside,
        0.5 * w * (z - z1) ** 2,
        0.0
    )

    M += np.where(
        after,
        w * (z2 - z1) *
        (z - (z1 + z2) / 2),
        0.0
    )

# ============================================================
# APPLIED MOMENTS
# ============================================================

for m in loads["moments"]:

    zm = m["z"]
    My = m["My"]

    M += np.where(
        z >= zm,
        My,
        0.0
    )

# ============================================================
# EQUILIBRIUM CHECKS
# ============================================================

sum_forces = (
    sum(p["Fx"] for p in loads["point_forces"])
    +
    sum(
        q["qx"] * (q["z2"] - q["z1"])
        for q in loads["distributed_loads"]
    )
)

print("Force equilibrium:")
print(f"ΣF = {sum_forces:.6f} N\n")
print("Moment at beam end:")
print(f"M(L) = {M[-1]:.6f} Nm\n")

# ============================================================
# PLOTS
# ============================================================

V = V/1000
M = M/1000

fig, axes = plt.subplots(
    2,
    1,
    figsize=(12, 8),
    sharex=True
)

# ------------------------------------------------------------
# Shear
# ------------------------------------------------------------

axes[0].plot(z, V, lw=2)

axes[0].fill_between(
    z,
    0,
    V,
    alpha=0.3
)

axes[0].axhline(0, color='k', lw=0.8)

axes[0].set_ylabel("V [kN]")
axes[0].grid(True)

# ------------------------------------------------------------
# Moment
# ------------------------------------------------------------

axes[1].plot(z, M, lw=2)

axes[1].fill_between(
    z,
    0,
    M,
    alpha=0.3
)

axes[1].axhline(0, color='k', lw=0.8)

axes[1].set_ylabel("M [kNm]")
axes[1].set_xlabel("z [m]")

axes[1].grid(True)

x_marks = {
    "nose cone start": 0.0,
    "passenger row 1": x1 - x_nose_cone,
    "passenger row 2": x2 - x_nose_cone,
    "passenger row 3": x3 - x_nose_cone,
    "cargo hold": x_cargo - x_nose_cone,
    "tail cone start": x_tail_cone - x_nose_cone,
    "wing": x_w - x_nose_cone
    # "horizontal tail": x_h - x_nose_cone,
    # "vertical tail": x_v - x_nose_cone,
}

x_ticks = list(x_marks.values())
x_labels = list(x_marks.keys())

for ax in axes:
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels, rotation=90, ha="center")

plt.tight_layout()
plt.show()

print(f"Maximum shear: {max(V)} kN")
print(f"Maximum moment: {max(M)} kNm")
