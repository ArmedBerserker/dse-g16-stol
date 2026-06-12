import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


def load_and_process_airfoil(filepath, chord):
    x_raw = []
    y_raw = []

    with open(filepath, 'r') as f:
        lines = f.readlines()

    for line in lines:
        parts = line.split()
        if len(parts) != 2:
            continue
        try:
            x_raw.append(float(parts[0]))
            y_raw.append(float(parts[1]))
        except ValueError:
            continue

    x_raw = np.array(x_raw) * chord
    y_raw = np.array(y_raw) * chord

    # Find the leading edge (minimum X coordinate) to split upper and lower surfaces
    le_idx = np.argmin(x_raw)

    # Split into upper and lower surfaces
    x_up, y_up = x_raw[:le_idx + 1], y_raw[:le_idx + 1]
    x_lo, y_lo = x_raw[le_idx:], y_raw[le_idx:]

    # Sort by X coordinate for interpolation purposes
    up_sort = np.argsort(x_up)
    lo_sort = np.argsort(x_lo)

    return x_up[up_sort], y_up[up_sort], x_lo[lo_sort], y_lo[lo_sort]


def generate_megson_idealization(filepath, chord, front_spar_pct, main_spar_pct, num_box_stringers, num_le_booms):
    """
    Generates a wing section idealization matching Megson's model.
    Outputs structured arrays of nodes (coordinates) and elements (connections)
    suitable for downstream structural calculations.
    """
    x_up, y_up, x_lo, y_lo = load_and_process_airfoil(filepath, chord)

    # Create interpolation functions
    interp_up = interp1d(x_up, y_up, kind='linear', fill_value="extrapolate")
    interp_lo = interp1d(x_lo, y_lo, kind='linear', fill_value="extrapolate")

    x_le_nose = float(np.min(x_up))
    x_front_spar = front_spar_pct * chord
    x_main_spar = main_spar_pct * chord

    # --- 1. Identify Spar Booms ---
    boom_tf = np.array([x_front_spar, float(interp_up(x_front_spar))])
    boom_bf = np.array([x_front_spar, float(interp_lo(x_front_spar))])
    boom_tm = np.array([x_main_spar, float(interp_up(x_main_spar))])
    boom_bm = np.array([x_main_spar, float(interp_lo(x_main_spar))])

    # --- 2. Generate Leading Edge (D-Cell) Booms ---
    # We step strictly from the nose back to the front spar
    x_le_space = np.linspace(x_le_nose, x_front_spar, num_le_booms + 1, endpoint=False)

    # Exclude index 0 for the lower nose array so we don't duplicate the exact nose point (0,0)
    le_upper_booms = np.array([[x, float(interp_up(x))] for x in x_le_space])
    le_lower_booms = np.array([[x, float(interp_lo(x))] for x in x_le_space])[1:]

    # --- 3. Generate Main Box Stringer Booms ---
    x_box_space = np.linspace(x_front_spar, x_main_spar, num_box_stringers + 2)[1:-1]
    box_upper_stringers = np.array([[x, float(interp_up(x))] for x in x_box_space])
    box_lower_stringers = np.array([[x, float(interp_lo(x))] for x in x_box_space])

    # --- 4. Apply Shift Vector (Origin at boom_bf) ---
    shift_vector = boom_bf.copy()

    boom_tf -= shift_vector
    boom_bf -= shift_vector  # (0,0)
    boom_tm -= shift_vector
    boom_bm -= shift_vector

    le_upper_booms -= shift_vector
    if len(le_lower_booms) > 0:
        le_lower_booms -= shift_vector

    if num_box_stringers > 0:
        box_upper_stringers -= shift_vector
        box_lower_stringers -= shift_vector

    # --- 5. Construct Global Ordered Nodes ---
    # We loop counter-clockwise: lower rear -> lower front -> nose -> upper front -> upper rear
    # This forms a single open profile outline, closed at the rear spar web.

    nodes_list = []

    # Track critical indices for internal elements (spar webs)
    idx_bm = len(nodes_list);
    nodes_list.append(boom_bm)

    if num_box_stringers > 0:
        for b in reversed(box_lower_stringers):
            nodes_list.append(b)

    idx_bf = len(nodes_list);
    nodes_list.append(boom_bf)

    if len(le_lower_booms) > 0:
        for b in reversed(le_lower_booms):
            nodes_list.append(b)

    # The nose point (index 0 of upper nose sequence)
    nodes_list.append(le_upper_booms[0])

    for b in le_upper_booms[1:]:
        nodes_list.append(b)

    idx_tf = len(nodes_list);
    nodes_list.append(boom_tf)

    if num_box_stringers > 0:
        for b in box_upper_stringers:
            nodes_list.append(b)

    idx_tm = len(nodes_list);
    nodes_list.append(boom_tm)

    nodes = np.array(nodes_list)

    # --- 6. Construct Elements (Connections) ---
    elements = []

    # Skin Panels along the outer skin profile chain
    for i in range(len(nodes) - 1):
        elements.append({"type": "skin", "nodes": [i, i + 1]})

    # Internal Webs
    elements.append({"type": "spar_front", "nodes": [idx_bf, idx_tf]})
    elements.append({"type": "spar_main", "nodes": [idx_bm, idx_tm]})

    # Optional closure skin panel for a multi-cell torque setup (Rear skin closure)
    # elements.append({"type": "skin", "nodes": [idx_tm, idx_bm]})

    # Keep track of structural identifiers for downstream script
    meta = {
        "idx_front_spar_bottom": idx_bf,
        "idx_front_spar_top": idx_tf,
        "idx_main_spar_bottom": idx_bm,
        "idx_main_spar_top": idx_tm
    }

    airfoil_profile = {
        "x_up": x_up - shift_vector[0],
        "y_up": y_up - shift_vector[1],
        "x_lo": x_lo - shift_vector[0],
        "y_lo": y_lo - shift_vector[1]
    }

    return nodes, elements, meta, airfoil_profile


if __name__ == "__main__":
    # --- CONFIGURATION ---
    dat_file = "NACA23012.dat"
    chord_length = 1.87
    front_spar_loc = 0.19
    main_spar_loc = 0.56
    num_le_booms = 3
    num_box_stringers = 4

    # Generate structured structural layout
    nodes, elements, meta, airfoil = generate_megson_idealization(
        filepath=dat_file,
        chord=chord_length,
        front_spar_pct=front_spar_loc,
        main_spar_pct=main_spar_loc,
        num_box_stringers=num_box_stringers,
        num_le_booms=num_le_booms
    )

    # --- PRINT SAMPLE FOR DOWNSTREAM SCRIPT ---
    print(f"Generated {len(nodes)} distinct boom coordinates.")
    print("Sample Nodes (First 3):")
    print(nodes[:3])
    print("\nSample Connectivity Elements (First 3):")
    print(elements[:3])

    # --- PLOTTING DATA ---
    plt.figure(figsize=(12, 5))
    plt.plot(airfoil['x_up'], airfoil['y_up'], 'gray', linestyle=':', alpha=0.5, label='Actual Airfoil Outline')
    plt.plot(airfoil['x_lo'], airfoil['y_lo'], 'gray', linestyle=':', alpha=0.5)

    # Plot segments dynamically using the connectivity table
    skin_labeled = False
    spar_labeled = False

    for el in elements:
        n_idx = el["nodes"]
        pt1, pt2 = nodes[n_idx[0]], nodes[n_idx[1]]

        if "skin" in el["type"]:
            plt.plot([pt1[0], pt2[0]], [pt1[1], pt2[1]], 'b-o', markersize=4,
                     label='Skin Panels / Booms' if not skin_labeled else "")
            skin_labeled = True
        else:
            plt.plot([pt1[0], pt2[0]], [pt1[1], pt2[1]], 'r-', linewidth=2,
                     label='Spar Webs' if not spar_labeled else "")
            spar_labeled = True

    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.axis('equal')
    plt.legend(loc='upper right')
    plt.show()