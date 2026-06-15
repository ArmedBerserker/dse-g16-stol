"""
Fuselage Barrel Mesh Generator
================================
Generates a 3D mesh of the aircraft fuselage barrel section (excluding nose and tail cone).
Cross-section: rectangle with different upper and lower corner radii.

Outputs:
  - fuselage_barrel.obj   : Wavefront OBJ mesh file
  - fuselage_barrel.png   : 3D visualization with frame positions highlighted
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import os

# =============================================================================
# PARAMETERS — edit these to match your aircraft
# =============================================================================

l_fus = 11.0
l_nc  = 3.4650000000000007
l_tc  = 4.429687500000001
l_section = l_fus - l_nc - l_tc   # barrel length (computed automatically)
x_main_spar = 4.4 + 0.19 * 1.87
x_rear_spar = 4.4 + 0.56 * 1.87
x_mlg = 6.01

# Cross-section geometry
h       = 1.7    # total height  [m]
w       = 1.45   # total width   [m]
r_upper = 0.55   # upper corner radius [m]
r_lower = 0.4    # lower corner radius [m]

# Mesh resolution
N_CIRC = 128     # circumferential points around the cross-section
# Stringers
STRINGER_SPACING = 0.19   # [m]

# Frame positions along the barrel (x-axis, measured from barrel start)
# Replace this array with your own computed positions.
# Example: uniform spacing of 0.45 m
FRAME_SPACING = 0.45   # [m]  used only when FRAME_POSITIONS is None
FRAME_POSITIONS = None  # set to e.g. np.array([0.0, 0.45, 0.90, ...]) to override
MINIMUM_FRAME_POSITIONS = np.array([l_nc, x_main_spar, x_rear_spar, x_mlg, l_tc])
FRAME_POSITIONS = np.array([3.465,  3.8951, 4.3252, 4.7553, 5.10125, 5.4472, 5.7286, 6.01, 6.29, 6.57])
# FRAMES_TO_ADD = []
# for i, frame in enumerate(MINIMUM_FRAME_POSITIONS):
#     if i < len(MINIMUM_FRAME_POSITIONS):
#         next_frame = MINIMUM_FRAME_POSITIONS[i+1]
#         if i == 0:
#             support_frame = next_frame - FRAME_SPACING / 2
#             spacing = np.arange(frame, support_frame+0.01, FRAME_SPACING)
#             FRAMES_TO_ADD.append(list(spacing))
#             FRAMES_TO_ADD.append(next_frame)
#         else:
#             fwd_support_frame = ...


# Output paths
OBJ_PATH = "fuselage_barrel.obj"
PNG_PATH = "fuselage_barrel.png"

# =============================================================================
# CROSS-SECTION PROFILE
# =============================================================================

def build_cross_section(w, h, r_upper, r_lower, n_points=128):
    """
    Build one closed cross-section polygon with:
      - upper-left  corner radius r_upper
      - upper-right corner radius r_upper
      - lower-right corner radius r_lower
      - lower-left  corner radius r_lower

    The shape is centred at (0, 0).
    Returns arrays (y, z) with shape (n_points,).
    """
    hw = w / 2.0
    hh = h / 2.0

    # Corner centres
    corners = [
        # (cy, cz,  angle_start, angle_end,  radius)  — CCW order starting bottom-right
        ( hw - r_lower, -hh + r_lower,  -np.pi/2, 0,          r_lower),   # bottom-right
        ( hw - r_upper,  hh - r_upper,   0,        np.pi/2,   r_upper),   # top-right
        (-hw + r_upper,  hh - r_upper,   np.pi/2,  np.pi,     r_upper),   # top-left
        (-hw + r_lower, -hh + r_lower,   np.pi,    3*np.pi/2, r_lower),   # bottom-left
    ]

    # Distribute points proportionally to arc length
    arc_lengths = [r * (end - start) for (_, _, start, end, r) in corners]
    total = sum(arc_lengths)
    pts_per_corner = [max(2, int(round(n_points * a / total))) for a in arc_lengths]
    # Adjust to hit exactly n_points
    diff = n_points - sum(pts_per_corner)
    pts_per_corner[0] += diff

    y_pts, z_pts = [], []
    for (cy, cz, a_start, a_end, r), n in zip(corners, pts_per_corner):
        angles = np.linspace(a_start, a_end, n, endpoint=False)
        y_pts.append(cy + r * np.cos(angles))
        z_pts.append(cz + r * np.sin(angles))

    y = np.concatenate(y_pts)
    z = np.concatenate(z_pts)
    return y, z

def build_stringers(y_cs, z_cs, spacing):

    # close the curve
    y_closed = np.append(y_cs, y_cs[0])
    z_closed = np.append(z_cs, z_cs[0])

    # locate top centre
    top_idx = np.argmax(z_cs)

    # rotate arrays so top centre is first
    y_rot = np.roll(y_cs, -top_idx)
    z_rot = np.roll(z_cs, -top_idx)

    y_rot = np.append(y_rot, y_rot[0])
    z_rot = np.append(z_rot, z_rot[0])

    # segment lengths
    ds = np.sqrt(
        np.diff(y_rot)**2 +
        np.diff(z_rot)**2
    )

    s = np.concatenate([[0], np.cumsum(ds)])

    perimeter = s[-1]

    s_targets = np.arange(
        0,
        perimeter,
        spacing
    )

    y_stringers = np.interp(
        s_targets,
        s,
        y_rot
    )

    z_stringers = np.interp(
        s_targets,
        s,
        z_rot
    )

    return y_stringers, z_stringers

# =============================================================================
# MESH CONSTRUCTION
# =============================================================================

def build_mesh(y_cs, z_cs, x_stations):
    """
    Extrude the cross-section along x_stations.
    Returns vertices (N,3) and quad faces (M,4) as index arrays.
    """
    n_circ = len(y_cs)
    n_stations = len(x_stations)

    # Vertices: shape (n_stations, n_circ, 3)
    verts = np.zeros((n_stations, n_circ, 3))
    for i, x in enumerate(x_stations):
        verts[i, :, 0] = x
        verts[i, :, 1] = y_cs
        verts[i, :, 2] = z_cs

    vertices = verts.reshape(-1, 3)

    # Quad faces connecting adjacent rings
    faces = []
    for i in range(n_stations - 1):
        for j in range(n_circ):
            j_next = (j + 1) % n_circ
            a = i     * n_circ + j
            b = i     * n_circ + j_next
            c = (i+1) * n_circ + j_next
            d = (i+1) * n_circ + j
            faces.append([a, b, c, d])

    return vertices, np.array(faces)


# =============================================================================
# OBJ EXPORT
# =============================================================================

def export_obj(vertices, faces, path):
    with open(path, "w") as f:
        f.write("# Fuselage barrel mesh\n")
        f.write(f"# Vertices: {len(vertices)}  Faces: {len(faces)}\n\n")
        for v in vertices:
            f.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
        f.write("\n")
        for face in faces:
            # OBJ is 1-indexed
            f.write("f " + " ".join(str(idx + 1) for idx in face) + "\n")
    print(f"[OBJ] Saved → {os.path.abspath(path)}")


# =============================================================================
# VISUALISATION
# =============================================================================

def visualize(vertices, 
              faces,
              y_cs,
              z_cs,
              x_stations,
              frame_positions,
              y_stringers,
              z_stringers,
              path):
    fig = plt.figure(figsize=(16, 7))

    # --- 3D view ---
    ax3d = fig.add_subplot(1, 2, 1, projection="3d")

    # Draw a subset of quads for speed (every 4th longitudinal strip)
    n_circ = len(y_cs)
    n_st   = len(x_stations)
    quad_verts = []
    # step = max(1, n_circ // 32)   # thin the circumferential quads for display
    step = max(1, n_circ // 128)   # was n_circ // 32 — shows twice as many strips
    # for i in range(0, n_st - 1, max(1, n_st // 40)):
    for i in range(n_st - 1):
        for j in range(0, n_circ, step):
            j_next = (j + 1) % n_circ
            a = i     * n_circ + j
            b = i     * n_circ + j_next
            c = (i+1) * n_circ + j_next
            d = (i+1) * n_circ + j
            quad_verts.append(vertices[[a, b, c, d]])

    # poly = Poly3DCollection(quad_verts, alpha=0.25, linewidth=0.0,
    #                         facecolor="steelblue", edgecolor="none")
    # poly = Poly3DCollection(quad_verts, alpha=0.20, linewidth=0.3,
    #                     facecolor="steelblue", edgecolor="steelblue")
    # poly = Poly3DCollection(
    #     quad_verts,
    #     alpha=1.0,
    #     linewidth=0,
    #     facecolor="steelblue",
    #     edgecolor="none"
    # )
    # poly = Poly3DCollection(
    #     quad_verts,
    #     facecolor="steelblue",
    #     edgecolor="none"
    # )
    poly = Poly3DCollection(
        quad_verts,
        facecolor="steelblue",
        edgecolor="none",
        alpha=1.0
    )

    poly.set_alpha(1.0)
    # ax3d.add_collection3d(poly)
    X = np.tile(x_stations[:,None], (1,n_circ))
    Y = np.tile(y_cs, (len(x_stations),1))
    Z = np.tile(z_cs, (len(x_stations),1))

    ax3d.plot_surface(
        X,
        Y,
        Z,
        color='steelblue',
        edgecolor='none',
        alpha=0.3
    )

    # Overlay frames
    # for xf in frame_positions:
    #     ax3d.plot(np.full(n_circ + 1, xf),
    #               np.append(y_cs, y_cs[0]),
    #               np.append(z_cs, z_cs[0]),
    #               color="tomato", linewidth=1.2)
        
    for xf in frame_positions:

        mask = y_cs < -0.05

        ax3d.plot(
            np.full(mask.sum(), xf),
            y_cs[mask],
            z_cs[mask],
            color='orangered',
            linewidth=2
        )

    # Overlay stringers

    # for ys, zs in zip(y_stringers, z_stringers):

    #     ax3d.plot(
    #         x_stations,
    #         np.full_like(x_stations, ys),
    #         np.full_like(x_stations, zs),
    #         color='limegreen',
    #         linewidth=1,
    #         alpha=0.4
    #     )

    n_half = len(y_stringers)//2

    for ys, zs in zip(
            y_stringers[:n_half],
            z_stringers[:n_half]
    ):

        ax3d.plot(
            x_stations,
            np.full_like(x_stations, ys),
            np.full_like(x_stations, zs),
            color='fuchsia',
            linewidth=1,
            alpha=0.4
        )

    # Nose and tail reference lines
    ax3d.plot(np.full(n_circ + 1, l_nc),
              np.append(y_cs, y_cs[0]),
              np.append(z_cs, z_cs[0]),
              color="gold", linewidth=1.5, linestyle="--", label="Barrel ends")
    ax3d.plot(np.full(n_circ + 1, x_stations[-1]),
              np.append(y_cs, y_cs[0]),
              np.append(z_cs, z_cs[0]),
              color="gold", linewidth=1.5, linestyle="--")

    ax3d.set_xlabel("x [m] (from nose)")
    ax3d.set_ylabel("y [m]")
    ax3d.set_zlabel("z [m]")
    ax3d.set_title("Fuselage barrel — 3D view\n(frames in red)", fontsize=11)
    # ax3d.set_box_aspect([x_stations[-1] - x_stations[0], w, h])
    ax3d.set_box_aspect((3, 1, 1))
    ax3d.view_init(elev=18, azim=-65)
    margin = 0.2

    ax3d.set_xlim(l_nc - 0.5, l_nc + l_section + 0.5)
    ax3d.set_ylim(-w/2 - margin, w/2 + margin)
    ax3d.set_yticks(np.linspace(-0.75,0.75,3))
    ax3d.tick_params(axis='y', labelsize=10)
    ax3d.set_zlim(-h/2 - margin, h/2 + margin)

    # --- Cross-section + frame spacing ---
    ax2 = fig.add_subplot(1, 2, 2)

    # Cross-section
    cs_closed_y = np.append(y_cs, y_cs[0])
    cs_closed_z = np.append(z_cs, z_cs[0])
    ax2.plot(cs_closed_y, cs_closed_z, "steelblue", linewidth=2)
    ax2.fill(y_cs, z_cs, alpha=0.15, color="steelblue")
    ax2.scatter(
        y_stringers,
        z_stringers,
        color='limegreen',
        s=35,
        zorder=10,
        label='Stringers'
    )

    # Annotate radii
    ax2.annotate(f"r_upper = {r_upper} m", xy=(w/2 - r_upper, h/2),
                 xytext=(w/2 - 0.7, h/2 -0.2),
                 arrowprops=dict(arrowstyle="->", color="gray"), fontsize=8, color="gray")
    ax2.annotate(f"r_lower = {r_lower} m", xy=(w/2 - r_lower, -h/2),
                 xytext=(w/2 - 0.7, -h/2 + 0.12),
                 arrowprops=dict(arrowstyle="->", color="gray"), fontsize=8, color="gray")
    ax2.annotate("", xy=(w/2, 0), xytext=(-w/2, 0),
                 arrowprops=dict(arrowstyle="<->", color="dimgray"))
    ax2.text(0, 0.04, f"w = {w} m", ha="center", fontsize=8, color="dimgray")
    ax2.annotate("", xy=(0, h/2), xytext=(0, -h/2),
                 arrowprops=dict(arrowstyle="<->", color="dimgray"))
    ax2.text(0.05, -h/6, f"h = {h} m", ha="left", va="center", fontsize=8, color="dimgray",
             rotation=90)

    ax2.set_aspect("equal")
    ax2.set_title("Cross-section  +  frame spacing diagram", fontsize=11)
    ax2.set_xlabel("y [m]")
    ax2.set_ylabel("z [m]")
    ax2.grid(True, linestyle=":", alpha=0.4)

    # # Inset: frame spacing bar along barrel length
    # ax_ins = ax2.inset_axes([0.0, -0.42, 1.0, 0.28])
    # ax_ins.set_xlim(0, x_stations[-1])
    # ax_ins.set_ylim(-0.5, 0.5)
    # ax_ins.axhline(0, color="steelblue", linewidth=3, alpha=0.4)
    # for i, xf in enumerate(frame_positions):
    #     ax_ins.axvline(xf, color="tomato", linewidth=1.5)
    #     if i % max(1, len(frame_positions) // 10) == 0:
    #         ax_ins.text(xf, 0.3, f"{xf:.2f}", fontsize=6, ha="center", color="tomato")
    # ax_ins.set_xlabel("x [m] along barrel", fontsize=8)
    # ax_ins.set_yticks([])
    # ax_ins.set_title(f"Frames ({len(frame_positions)} total)", fontsize=8)

    plt.tight_layout(rect=[0, 0.0, 1, 1])
    plt.savefig(path, dpi=150, bbox_inches="tight")
    print(f"[PNG] Saved → {os.path.abspath(path)}")
    plt.show()


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    print(f"Barrel length: {l_section:.4f} m")

    # Build cross-section
    y_cs, z_cs = build_cross_section(w, h, r_upper, r_lower, n_points=N_CIRC)
    y_stringers, z_stringers = build_stringers(
        y_cs,
        z_cs,
        STRINGER_SPACING
    )

    print(f"Stringer count: {len(y_stringers)}")

    # Frame positions
    if FRAME_POSITIONS is not None:
        frame_pos = np.asarray(FRAME_POSITIONS)
    else:
        frame_pos = np.arange(0, l_section + 1e-9, FRAME_SPACING)

    print(f"Frame count  : {len(frame_pos)}")
    print(f"Frame positions [m]: {np.round(frame_pos, 4)}")

    # Longitudinal stations = frame positions + fine fill for smooth mesh
    N_LONG = max(len(frame_pos) * 4, 80)
    x_fill = np.linspace(l_nc, l_nc + l_section, N_LONG)
    x_stations = np.unique(np.concatenate([x_fill, frame_pos]))

    # Build mesh
    vertices, faces = build_mesh(y_cs, z_cs, x_stations)
    print(f"Mesh: {len(vertices)} vertices, {len(faces)} quads")

    # Export
    export_obj(vertices, faces, OBJ_PATH)

    # Visualise
    # visualize(vertices, faces, y_cs, z_cs, x_stations, frame_pos, PNG_PATH)
    visualize(
        vertices,
        faces,
        y_cs,
        z_cs,
        x_stations,
        frame_pos,
        y_stringers,
        z_stringers,
        PNG_PATH
    )
