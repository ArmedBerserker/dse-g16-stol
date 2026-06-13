import numpy as np
import pyvista as pv

class FuselageGenerator:
    def __init__(self, radius=0.8, length=2.85):
        self.R = radius
        self.L = length
        self.points = None
        self.mesh = pv.MultiBlock()

    def generate(self, n_stringers=16, rib_positions=None,
                 sub_theta=4, sub_z=10, rib_width=0.2,
                 skin_props=None,
                 stringer_props=None,
                 rib_props=None):

        if rib_positions is None:
            rib_positions = [0, 3, 6, 9, 12]
        rib_pos = np.sort(np.array(rib_positions))
        n_ribs = len(rib_pos)

        # 1. Grid Setup
        n_t = n_stringers * sub_theta
        theta = np.linspace(-np.pi / 2, 1.5 * np.pi, n_t, endpoint=False)

        z_coords = []
        for i in range(n_ribs - 1):
            z_segment = np.linspace(rib_pos[i], rib_pos[i + 1], sub_z, endpoint=False)
            z_coords.extend(z_segment)
        z_coords.append(rib_pos[-1])
        z_coords = np.array(z_coords)
        n_z = len(z_coords)

        # 2. Point Generation
        Theta, Z = np.meshgrid(theta, z_coords, indexing="ij")
        outer_pts = np.column_stack(((self.R * np.cos(Theta)).ravel(), (self.R * np.sin(Theta)).ravel(), Z.ravel()))

        inner_r = self.R - rib_width
        rib_z_indices = np.arange(0, n_z, sub_z)
        inner_pts = np.column_stack(((inner_r * np.cos(Theta[:, rib_z_indices])).ravel(),
                                     (inner_r * np.sin(Theta[:, rib_z_indices])).ravel(),
                                     Z[:, rib_z_indices].ravel()))
        self.points = np.vstack((outer_pts, inner_pts))

        outer_map = np.arange(len(outer_pts)).reshape((n_t, n_z))
        inner_map = np.arange(len(outer_pts), len(self.points)).reshape((n_t, n_ribs))

        # 3. SKIN
        skin_cells = []
        for i in range(n_t):
            for j in range(n_z - 1):
                skin_cells.extend([4, outer_map[i, j], outer_map[(i + 1) % n_t, j],
                                   outer_map[(i + 1) % n_t, j + 1], outer_map[i, j + 1]])
        self.mesh["Skin"] = pv.PolyData(self.points, faces=skin_cells)
        self.mesh["Skin"].cell_data["E"] = np.full(self.mesh["Skin"].n_cells, skin_props['E'])
        self.mesh["Skin"].cell_data["thickness"] = np.full(self.mesh["Skin"].n_cells, skin_props['t'])

        # 4. STRINGERS
        stringer_cells = []
        s_E, s_A, s_I11, s_I22, s_J = [], [], [], [], []
        stringer_theta_indices = np.arange(0, n_t, sub_theta)

        for idx, i in enumerate(stringer_theta_indices):
            # Property retrieval with fallbacks if properties are uniform vs varying
            curr_E = stringer_props['E'][idx] if isinstance(stringer_props['E'], (list, np.ndarray)) else \
                stringer_props['E']
            curr_A = stringer_props['A'][idx] if isinstance(stringer_props['A'], (list, np.ndarray)) else \
                stringer_props['A']
            curr_I11 = stringer_props['I11'][idx] if isinstance(stringer_props.get('I11', 0),
                                                                (list, np.ndarray)) else stringer_props.get('I11', 1e-8)
            curr_I22 = stringer_props['I22'][idx] if isinstance(stringer_props.get('I22', 0),
                                                                (list, np.ndarray)) else stringer_props.get('I22', 1e-8)
            curr_J = stringer_props['J'][idx] if isinstance(stringer_props.get('J', 0),
                                                            (list, np.ndarray)) else stringer_props.get('J', 1e-9)

            for j in range(n_z - 1):
                stringer_cells.extend([2, outer_map[i, j], outer_map[i, j + 1]])
                s_E.append(curr_E)
                s_A.append(curr_A)
                s_I11.append(curr_I11)
                s_I22.append(curr_I22)
                s_J.append(curr_J)

        self.mesh["Stringers"] = pv.PolyData(self.points, lines=stringer_cells)
        self.mesh["Stringers"].cell_data["E"] = np.array(s_E)
        self.mesh["Stringers"].cell_data["area"] = np.array(s_A)
        self.mesh["Stringers"].cell_data["I11"] = np.array(s_I11)
        self.mesh["Stringers"].cell_data["I22"] = np.array(s_I22)
        self.mesh["Stringers"].cell_data["J"] = np.array(s_J)

        # 5. RIBS
        rib_cells = []
        r_E, r_t = [], []
        for r_idx, z_idx in enumerate(rib_z_indices):
            curr_E = rib_props['E'][r_idx] if isinstance(rib_props['E'], (list, np.ndarray)) else rib_props['E']
            curr_t = rib_props['t'][r_idx] if isinstance(rib_props['t'], (list, np.ndarray)) else rib_props['t']

            for i in range(n_t):
                rib_cells.extend([4, outer_map[i, z_idx], outer_map[(i + 1) % n_t, z_idx],
                                  inner_map[(i + 1) % n_t, r_idx], inner_map[i, r_idx]])
                r_E.append(curr_E)
                r_t.append(curr_t)

        self.mesh["Ribs"] = pv.PolyData(self.points, faces=rib_cells)
        self.mesh["Ribs"].cell_data["E"] = np.array(r_E)
        self.mesh["Ribs"].cell_data["thickness"] = np.array(r_t)

    def show(self):
        pl = pv.Plotter()

        pl.add_mesh(
            self.mesh["Skin"],
            color="gray",
            opacity=0.4,
            show_edges=True,
            edge_color="black",
            line_width=1
        )

        pl.add_mesh(
            self.mesh["Ribs"],
            scalars="thickness",
            cmap="viridis",
            show_edges=True,
            edge_color="black"
        )

        pl.add_mesh(
            self.mesh["Stringers"],
            color="red",
            line_width=4,
            render_lines_as_tubes=True
        )

        pl.show()

gen = FuselageGenerator()

my_rib_positions = [0, 0.50, 1, 1.45, 2, 2.45, 2.85]
my_rib_thicknesses = [0.020, 0.005, 0.005, 0.005, 0.005, 0.010, 0.020]

gen.generate(
    n_stringers=16,
    rib_positions=my_rib_positions,
    sub_theta=4, sub_z=10, rib_width=0.05,
    rib_props={'E': 70e9, 't': my_rib_thicknesses},
    skin_props={'E': 70e9, 't': 0.0015},
    stringer_props={
        'E': 70e9,
        'A': 0.0004,
        'I11': 1.5e-8,  # Radial Direction (strong axis)
        'I22': 2.0e-8,  # Tangential Direction (weak axis)
        'J': 1.0e-9
    }
)


gen.show()