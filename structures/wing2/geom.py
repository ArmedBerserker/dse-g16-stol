import numpy as np
from airfoil import generate_megson_idealization
import matplotlib.pyplot as plt


def analyze_section_properties(nodes, elements, t_skin, t_web, A_spar, A_stringer, meta):
    # Identify which nodes are spar caps vs regular stringers using meta indices
    spar_nodes = {
        meta["idx_front_spar_bottom"], meta["idx_front_spar_top"],
        meta["idx_main_spar_bottom"], meta["idx_main_spar_top"]
    }

    # --- STEP 1: Calculate Centroid based on Direct Stress Carrying Booms ---
    total_boom_area = 0.0
    sum_Qx = 0.0  # First moment about X
    sum_Qy = 0.0  # First moment about Y

    boom_properties = []

    for idx, node in enumerate(nodes):
        x, y = node[0], node[1]
        # Assign area based on whether it is a heavy spar boom or a regular stringer
        area = A_spar if idx in spar_nodes else A_stringer

        total_boom_area += area
        sum_Qy += x * area
        sum_Qx += y * area

        boom_properties.append({
            "idx": idx,
            "x_global": x,
            "y_global": y,
            "area": area
        })

    # This finds the centroid relative to the bottom left spar cap (which is 0,0)
    X_bar = sum_Qy / total_boom_area
    Y_bar = sum_Qx / total_boom_area

    # --- STEP 2: Calculate Centroidal Moments of Inertia & Form Table ---
    Ixx = 0.0
    Iyy = 0.0
    Ixy = 0.0

    print("\n" + "=" * 85)
    print(
        f"{'Boom':<6} | {'x (mm)':<12} | {'y (mm)':<12} | {'B (mm2)':<12} | {'Ixx = By^2 (mm4)':<18} | {'Iyy = Bx^2 (mm4)':<18}")
    print("=" * 85)

    for bp in boom_properties:
        # Shift to true centroidal coordinate system (Origin at X_bar, Y_bar)
        x_c = bp["x_global"] - X_bar
        y_c = bp["y_global"] - Y_bar
        B = bp["area"]

        # Calculate contributions
        Ixx_boom = B * (y_c ** 2)
        Iyy_boom = B * (x_c ** 2)
        Ixy_boom = B * x_c * y_c

        Ixx += Ixx_boom
        Iyy += Iyy_boom
        Ixy += Ixy_boom

        # Convert to mm and mm2 for the tabular output
        x_mm = x_c * 1e3
        y_mm = y_c * 1e3
        B_mm2 = B * 1e6
        Ixx_mm4 = Ixx_boom * 1e12
        Iyy_mm4 = Iyy_boom * 1e12

        print(f"{bp['idx']:<6} | {x_mm:>12.2f} | {y_mm:>12.2f} | {B_mm2:>12.2f} | {Ixx_mm4:>18.2e} | {Iyy_mm4:>18.2e}")

    print("=" * 85)

    # Global property outputs in standard SI and mm units
    print(f"  Centroid (from Bottom-Left Spar): X_bar = {X_bar:.4f} m, Y_bar = {Y_bar:.4f} m")
    print(f"  Total Boom Area: {total_boom_area * 1e6:.2f} mm2")
    print(f"  Ixx = {Ixx * 1e12:.4e} mm4 ({Ixx:.4e} m4)")
    print(f"  Iyy = {Iyy * 1e12:.4e} mm4 ({Iyy:.4e} m4)")
    print(f"  Ixy = {Ixy * 1e12:.4e} mm4 ({Ixy:.4e} m4)")
    print("=" * 85 + "\n")

    # --- STEP 3: Multicell Torsion Properties Accumulation ---
    Ae1 = 0.0
    Ae2 = 0.0

    # Line integrals around individual loops:oint(ds/t)
    integral_ds_t_1 = 0.0
    integral_ds_t_2 = 0.0
    integral_ds_t_shared = 0.0

    for el in elements:
        # Determine geometric wall thickness
        t_element = t_web if "spar" in el["type"] else t_skin
        ds_over_t = el["length"] / t_element

        # --- Enclosed Area Accumulations ---
        if 1 in el["cells"]:
            # Skin flows naturally counter-clockwise, internal web front spar closes it running downwards
            Ae1 += el["swept_area"] if el["type"] == "skin" else -el["swept_area"]
            integral_ds_t_1 += ds_over_t

        if 2 in el["cells"]:
            if el["type"] == "skin":
                Ae2 += el["swept_area"]
            elif el["type"] == "spar_front":
                Ae2 += el["swept_area"]  # Front spar runs upwards relative to Cell 2 CCW flow
            elif el["type"] == "spar_main":
                Ae2 -= el["swept_area"]  # Main spar closes Cell 2 on the right heading downwards
            integral_ds_t_2 += ds_over_t

        # --- Track shared wall line integrals ---
        if el["type"] == "spar_front":
            integral_ds_t_shared = ds_over_t

    Ae1 = abs(Ae1)
    Ae2 = abs(Ae2)

    total_perimeter = meta["airfoil_perimeter"]
    h_front = meta["front_spar_height"]
    h_main = meta["main_spar_height"]
    A_frontspar = h_front * t_web
    A_rearspar = h_main * t_web

    # Secondary print block for torsion diagnostics
    print(f"  Cell 1 Enclosed Area (Ae1): {Ae1:.6f} m2  ({Ae1 * 1e6:.2f} mm2)")
    print(f"  Cell 2 Enclosed Area (Ae2): {Ae2:.6f} m2  ({Ae2 * 1e6:.2f} mm2)")
    print(f"  Loop 1 Line Integral  oint(ds/t): {integral_ds_t_1:.4f}")
    print(f"  Loop 2 Line Integral  oint(ds/t): {integral_ds_t_2:.4f}")
    print(f"  Shared Wall Integral  int(ds/t) : {integral_ds_t_shared:.4f}")
    print(total_perimeter)

    return X_bar, Y_bar, total_boom_area, Ixx, Iyy, Ixy, Ae1, Ae2, integral_ds_t_1, integral_ds_t_2, integral_ds_t_shared, A_frontspar, A_rearspar


if __name__ == "__main__":
    dat_file = "NACA23012.dat"
    chord_length = 1.87

    nodes, elements, meta, _ = generate_megson_idealization(
        filepath=dat_file,
        chord=chord_length,
        front_spar_pct=0.16,
        main_spar_pct=0.56,
        num_box_stringers=8,
        num_le_booms=4
    )

    thickness_skin = 0.001
    thickness_web = 0.003
    area_spar_cap = 0.0005
    area_stringer = 7e-5

    # Call structural analysis with updated multi-variable unpack
    outputs = analyze_section_properties(
        nodes=nodes,
        elements=elements,
        t_skin=thickness_skin,
        t_web=thickness_web,
        A_spar=area_spar_cap,
        A_stringer=area_stringer,
        meta=meta
    )

    (
        X_cg, Y_cg, total_mat_area, Ixx, Iyy, Ixy,
        Ae1, Ae2, ds_t1, ds_t2, ds_t_shared,
        Afs, Ars
    ) = outputs

    # --- Plotting Configuration ---
    plt.figure(figsize=(14, 6))

    skin_labeled = False
    web_labeled = False
    for el in elements:
        p1, p2 = nodes[el["nodes"][0]], nodes[el["nodes"][1]]
        if 'spar' in el["type"]:
            plt.plot([p1[0], p2[0]], [p1[1], p2[1]], color='orange', linewidth=2, linestyle='--',
                     label='Spar Webs' if not web_labeled else "")
            web_labeled = True
        else:
            plt.plot([p1[0], p2[0]], [p1[1], p2[1]], color='gray', linewidth=1, alpha=0.7,
                     label='Skin Panels' if not skin_labeled else "")
            skin_labeled = True

    spar_nodes = {
        meta["idx_front_spar_bottom"], meta["idx_front_spar_top"],
        meta["idx_main_spar_bottom"], meta["idx_main_spar_top"]
    }

    spar_x, spar_y, spar_sizes = [], [], []
    str_x, str_y, str_sizes = [], [], []

    for idx, (x, y) in enumerate(nodes):
        if idx in spar_nodes:
            spar_x.append(x)
            spar_y.append(y)
            area_mm2 = area_spar_cap * 1e6
            spar_sizes.append(area_mm2)
        else:
            str_x.append(x)
            str_y.append(y)
            area_mm2 = area_stringer * 1e6
            str_sizes.append(area_mm2)

    plt.scatter(str_x, str_y, s=str_sizes, color='royalblue', edgecolor='k', zorder=3)
    plt.scatter(spar_x, spar_y, s=spar_sizes, color='red', edgecolor='k', zorder=4)

    plt.plot(X_cg, Y_cg, 'gX', markersize=14, markeredgecolor='black', zorder=5,
             label=f'Centroid ({X_cg * 1e3:.1f}, {Y_cg * 1e3:.1f}) mm')

    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.axis('equal')
    plt.legend(loc='upper right', scatterpoints=1)
    plt.tight_layout()
    plt.show()