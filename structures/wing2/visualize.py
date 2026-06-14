import numpy as np
import pyvista as pv
from scipy.interpolate import interp1d

# Import your custom structural layout engine from airfoil.py
from airfoil import generate_megson_idealization, load_and_process_airfoil


def build_3d_wingbox_mesh(chord, front_spar_pct, main_spar_pct, filepath, wing_bays, y_root=0.725):
    """
    Processes the wing bay-by-bay to build 3D surfaces for spars,
    and line sets for longitudinal stringers. (Skin is now handled globally).
    """
    plotter = pv.Plotter()
    plotter.set_background("white")

    current_y = y_root
    front_spar_track = []
    main_spar_track = []

    # Loop through each bay to construct local meshes
    for idx, bay in enumerate(wing_bays):
        y_start = current_y
        y_end = bay['y_end']

        # Call the layout engine from airfoil.py
        nodes_2d, elements, meta, _ = generate_megson_idealization(
            filepath=filepath,
            chord=chord,
            front_spar_pct=front_spar_pct,
            main_spar_pct=main_spar_pct,
            num_box_stringers=bay['num_box_stringers'],
            num_le_booms=bay['num_le_booms']
        )

        num_nodes = len(nodes_2d)

        # X = Chordwise, Y = Spanwise, Z = Thickness
        pts_start = np.zeros((num_nodes, 3))
        pts_start[:, 0] = nodes_2d[:, 0]
        pts_start[:, 1] = y_start
        pts_start[:, 2] = nodes_2d[:, 1]

        pts_end = np.zeros((num_nodes, 3))
        pts_end[:, 0] = nodes_2d[:, 0]
        pts_end[:, 1] = y_end
        pts_end[:, 2] = nodes_2d[:, 1]

        # Track spar centerlines
        if idx == 0:
            front_spar_track.append(pts_start[meta["idx_front_spar_top"]])
            main_spar_track.append(pts_start[meta["idx_main_spar_top"]])
        front_spar_track.append(pts_end[meta["idx_front_spar_top"]])
        main_spar_track.append(pts_end[meta["idx_main_spar_top"]])

        bay_pts = np.vstack([pts_start, pts_end])

        front_spar_faces = []
        main_spar_faces = []

        # Build 3D Quad elements connecting bay panels
        for el in elements:
            n1, n2 = el["nodes"]
            quad = [4, n1, n2, n2 + num_nodes, n1 + num_nodes]

            if el["type"] == "spar_front":
                front_spar_faces.extend(quad)
            elif el["type"] == "spar_main":
                main_spar_faces.extend(quad)

        # Render Spar Webs
        if front_spar_faces:
            f_mesh = pv.PolyData(bay_pts, front_spar_faces)
            plotter.add_mesh(f_mesh, color="red", opacity=0.85, show_edges=True)
        if main_spar_faces:
            m_mesh = pv.PolyData(bay_pts, main_spar_faces)
            plotter.add_mesh(m_mesh, color="darkred", opacity=0.85, show_edges=True)

        # Render Stringers / Booms
        stringer_lines = []
        for i in range(num_nodes):
            stringer_lines.extend([2, i, i + num_nodes])
        stringers_mesh = pv.PolyData(bay_pts)
        stringers_mesh.lines = np.array(stringer_lines)
        plotter.add_mesh(stringers_mesh, color="black", line_width=2)

        current_y = y_end

    return plotter, np.array(front_spar_track), np.array(main_spar_track)


def add_full_wing_skin(plotter, chord, front_spar_pct, filepath, y_root, y_tip):

    # 1. Load the raw aerodynamic profile
    x_up, y_up, x_lo, y_lo = load_and_process_airfoil(filepath, chord)

    # 2. Calculate the exact same origin shift used for the structural wingbox
    interp_lo = interp1d(x_lo, y_lo, kind='linear', fill_value="extrapolate")
    shift_x = front_spar_pct * chord
    shift_z = float(interp_lo(shift_x))

    # 3. Stitch upper and lower surfaces into a continuous loop
    x_profile = np.concatenate([x_up[::-1], x_lo[1:]])
    z_profile = np.concatenate([y_up[::-1], y_lo[1:]])

    # Apply the structural shift
    x_shifted = x_profile - shift_x
    z_shifted = z_profile - shift_z
    num_pts = len(x_shifted)

    # 4. Create 3D points at the root and the absolute tip
    pts_root = np.column_stack((x_shifted, np.full(num_pts, y_root), z_shifted))
    pts_tip = np.column_stack((x_shifted, np.full(num_pts, y_tip), z_shifted))
    all_pts = np.vstack((pts_root, pts_tip))

    # 5. Build Quad faces connecting the root profile to the tip profile
    faces = []
    for i in range(num_pts - 1):
        n1 = i
        n2 = i + 1
        n3 = n2 + num_pts
        n4 = n1 + num_pts
        faces.extend([4, n1, n2, n3, n4])

    # 6. Render the full Outer Mold Line (OML) as a glassy shell
    full_skin_mesh = pv.PolyData(all_pts, faces)
    plotter.add_mesh(
        full_skin_mesh,
        color="gray",
        opacity=0.3,
        show_edges=False,
        label="Outer Wing Skin"
    )


def add_subsystems_and_masses(plotter, chord, front_spar_pct, main_spar_pct, filepath,
                              actuators, engine, fuel_tanks):
    x_up, y_up, x_lo, y_lo = load_and_process_airfoil(filepath, chord)
    interp_up = interp1d(x_up, y_up, kind='linear', fill_value="extrapolate")
    interp_lo = interp1d(x_lo, y_lo, kind='linear', fill_value="extrapolate")

    x_front_spar = front_spar_pct * chord
    shift_x = x_front_spar
    shift_z = float(interp_lo(x_front_spar))

    def transform_to_3d(x_global, y_global, z_offset=0.0):
        x_3d = x_global - shift_x

        # Clamping x_global to [0, chord] prevents Z-axis scale destruction
        x_clamped = np.clip(x_global, 0.0, chord)

        z_mid = 0.5 * (float(interp_up(x_clamped)) + float(interp_lo(x_clamped)))
        z_3d = z_mid - shift_z + z_offset
        return [x_3d, y_global, z_3d]

    added_labels = set()

    # Add Actuators (tucked safely underneath the control surfaces)
    for act in actuators:
        loc = transform_to_3d(act['x_loc'], act['y'], z_offset=-0.04)
        sphere = pv.Sphere(radius=0.06, center=loc)

        lbl = "Actuator" if "Actuator" not in added_labels else None
        plotter.add_mesh(sphere, color="chartreuse", label=lbl)
        added_labels.add("Actuator")

    # Add Engine (manually hang the engine pod below the wing using z_offset)
    eng_loc = transform_to_3d(engine['x_loc'], engine['y_loc'], z_offset=-0.25)
    engine_cylinder = pv.Cylinder(center=eng_loc, radius=0.15, height=0.5, direction=(0, 1, 0))
    plotter.add_mesh(engine_cylinder, color="darkorange", label="Engine")

    # Add Fuel Tanks (centered nicely inside the shifted spar limits)
    for tank in fuel_tanks:
        x_start_box = 0.0
        x_end_box = (main_spar_pct - front_spar_pct) * chord

        tank_box = pv.Box(bounds=(x_start_box, x_end_box, tank['y_start'], tank['y_end'], -0.05, 0.05))
        plotter.add_mesh(tank_box, color="teal", opacity=0.35, label="Fuel Tank Space")


# --- RUNTIME EXECUTION ---
if __name__ == "__main__":
    # Parameters
    chord = 1.87
    front_spar_pct = 0.15
    main_spar_pct = 0.50
    filepath = "NACA23012.dat"

    wing_bays = [
        {'y_end': 1.5, 'num_box_stringers': 14, 'num_le_booms': 2},
        {'y_end': 3.0, 'num_box_stringers': 13, 'num_le_booms': 2},
        {'y_end': 4.5, 'num_box_stringers': 9, 'num_le_booms': 1},
        {'y_end': 5.8, 'num_box_stringers': 6, 'num_le_booms': 0},
        {'y_end': 7.0, 'num_box_stringers': 2, 'num_le_booms': 0},
        {'y_end': 8.4, 'num_box_stringers': 1, 'num_le_booms': 1}
    ]

    flap_actuators = [{'y': 1.5, 'x_loc': 0.6 * chord}, {'y': 2.5, 'x_loc': 0.6* chord},
                      {'y': 3.5, 'x_loc': 0.6 * chord}]
    aileron_actuators = [{'y': 6.0, 'x_loc': 0.7 * chord}, {'y': 7.0, 'x_loc': 0.7 * chord}]
    spoiler_actuators = [{'y': 4.2, 'x_loc': 0.3 * chord}, {'y': 5.0, 'x_loc': 0.3 * chord}]
    all_actuators = flap_actuators + aileron_actuators + spoiler_actuators

    engine_config = {'y_loc': 2.2, 'x_loc': -0.1 * chord}
    fuel_tanks = [{'y_start': 0.725, 'y_end': 4.0}]
    y_root = 0.725

    # 1. Build the Internal Structure
    plotter, front_spar_line, main_spar_line = build_3d_wingbox_mesh(
        chord=chord,
        front_spar_pct=front_spar_pct,
        main_spar_pct=main_spar_pct,
        filepath=filepath,
        wing_bays=wing_bays,
        y_root=y_root
    )

    # 2. Add the Complete Aerodynamic Skin
    y_tip = wing_bays[-1]['y_end']
    add_full_wing_skin(
        plotter=plotter,
        chord=chord,
        front_spar_pct=front_spar_pct,
        filepath=filepath,
        y_root=y_root,
        y_tip=y_tip
    )

    # 3. Add Subsystems
    add_subsystems_and_masses(
        plotter=plotter, chord=chord, front_spar_pct=front_spar_pct, main_spar_pct=main_spar_pct,
        filepath=filepath, actuators=all_actuators, engine=engine_config, fuel_tanks=fuel_tanks
    )

    # 4. Render Setup
    plotter.add_lines(front_spar_line, color="red", width=3, label="Front Spar Axis", connected=True)
    plotter.add_lines(main_spar_line, color="maroon", width=3, label="Rear Spar Axis", connected=True)

    plotter.add_legend(bcolor='white', border=True)
    #plotter.add_axes()
    #plotter.show_grid()

    plotter.show()