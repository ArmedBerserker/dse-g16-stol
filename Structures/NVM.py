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
axes[1].set_xlabel("x (nosecone to tailcone)")

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

M *= 1.5 * 3.8
V *= 1.5 * 3.8
# Stringer location
y_stringers = np.array([-0.175,-0.36123634, -0.52543939, -0.64822498, -0.71509851, -0.72457458, -0.72469634, -0.72481809, -0.72493985, -0.71342573, -0.62681825, -0.47341679, -0.28708811, -0.09708821,  0.0929117,   0.2829116,   0.46952622,  0.62403519, 0.71234222,  0.72442211,  0.7245884,   0.72475468,  0.72492096,  0.71589258, 0.6503415,   0.52863656,  0.36514891,  0.17917666, -0.01082317]) 
z_stringers = np.array([ 0.85, 0.81744688,  0.72377284,  0.58004631,  0.40324118,  0.21386786, 0.0238679,  -0.16613206, -0.35613202, -0.54520295, -0.71228349, -0.82126519, -0.84938475, -0.84957573, -0.84976671, -0.84995769, -0.82281082, -0.71541079, -0.54924665, -0.36030855, -0.17030862,  0.01969131,  0.20969123,  0.39913536, 0.57643957,  0.72107721,  0.81597034,  0.84952987,  0.84978207])

Booms = np.arange(0, len(y_stringers), 1)
if len(Booms) != len(y_stringers):
    raise ValueError(f"Booms != y_stringers")

t_skin = 0.3*1e-3  # m
stringer_pitch = 0.2  # m
ratio1 = np.zeros_like(z_stringers)
ratio2 = np.zeros_like(z_stringers)
n = len(y_stringers)
for i, z in enumerate(z_stringers):
    ratio1[i] = z_stringers[(i + 1) % n] / z
    ratio2[i] = z_stringers[(i - 1) % n] / z
B = 70*1e-6 * np.ones_like(y_stringers) + t_skin * stringer_pitch / 6 * ((2 + ratio1) + (2 + ratio2))
# B = 70*1e-6 * np.ones_like(y_stringers) + t_skin * stringer_pitch / 6 * (2 + ratio1)
D_Iyy = B * z_stringers**2
M_max = np.max(M) * 1000
sigma_x = M_max * z_stringers / np.sum(D_Iyy)
# print(sigma_x)
sigma_yield = 345000000

print(f'max normal stress: {np.max(np.abs(sigma_x))}')

if np.max(np.abs(sigma_x)) > sigma_yield:
    raise ValueError(f'Structure yields')
else:
    min_margin_stress = (sigma_yield - np.max(np.abs(sigma_x))) / np.max(np.abs(sigma_x)) * 100
    print(f'Minimum stress margin {min_margin_stress}%')

dq_max = -np.max(V*1000) / np.sum(D_Iyy) * B * z_stringers
dq_min = -np.min(V*1000) / np.sum(D_Iyy) * B * z_stringers

q_max = np.zeros_like(dq_max)
q_min = np.zeros_like(dq_min)

T = 5750 * (2.04 / 2 + 0.8)  # VT force * arm
q_torque = T / 2 / 2.26649330114264

for i in range(14, 41):
    q_max[i % n] = q_max[(i - 1) % n] + dq_max[i % n]
    q_min[i % n] = q_min[(i - 1) % n] + dq_min[i % n]

q1 = q_max + q_torque
q2 = q_max - q_torque
q3 = q_min + q_torque
q4 = q_min - q_torque

q_max_overall = 1000 * np.max([np.max(np.abs(q1)), np.max(np.abs(q2)), np.max(np.abs(q3)), np.max(np.abs(q4)), np.max(np.abs(q_max)), np.max(np.abs(q_min)), q_torque])

tau_yield = 207000000
if np.abs(q_max_overall) > tau_yield:
    raise ValueError(f'Structure yields under shear stress')
else:
    min_margin_shear_stress = (tau_yield - q_max_overall) / q_max_overall * 100
    print(f'Minimum shear stress margin {min_margin_shear_stress}%')

print(f'\n \t OVERVIEW: \n \n max normal stress: {np.max(np.abs(sigma_x))} \n margin normal stress: {min_margin_stress} \n max shear stress: {q_max_overall} \n margin shear stress: {min_margin_shear_stress}')

