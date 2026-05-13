import numpy as np
import pandas as pd
import os

###############################################################################
# COORDINATE SYSTEM & SIGN CONVENTION DEFINITIONS
# -----------------------------------------------------------------------------
# GLOBAL FRAME (XFLR5):
#   X+ : Aft (Trailing Edge)
#   Y+ : Right (Outboard)
#   Z+ : Up (Global Lift)
#
# LOCAL STATION FRAME (Right-Handed System):
#   c_u (Chord)  : Points AFT       [Similar with Global X+]
#   s_u (Span)   : Points OUTBOARD  [Similar with Global Y+]
#   n_u (Normal) : Points DOWN      [Opposite to Global Z+] -> (s_u x c_u)
#
# OUTPUT SIGN CONVENTIONS:
#   Fn (Normal Force) : Positive (+) is DOWN, Negative (-) is UP (Lift).
#   Torsion           : Positive (+) is NOSE-UP.
#   Fx, Fy, Fz        : Preserved in Global orientation.
###############################################################################

class XFLR5Parser:
    def __init__(self, filepath, air_density=1.225):
        self.filepath = filepath
        self.rho = air_density
        self.v_inf = 0.0
        self.station_chords = {}
        self.strips_raw = []

    def parse(self):
        if not os.path.exists(self.filepath):
            return None

        with open(self.filepath, 'r', encoding='utf-8-sig', errors='ignore') as f:
            lines = f.readlines()

        mode = None
        current_strip = []
        for line in lines:
            line = line.replace('\xa0', ' ').strip()
            if not line: continue

            if "QInf" in line:
                try:
                    self.v_inf = float([p for p in line.split() if p.replace('.', '', 1).isdigit()][0])
                except:
                    pass
            elif "Main Wing" in line and "Cp" not in line and "y-span" not in line:
                mode = "GEOM"
            elif "Cp Coefficients" in line:
                mode = "CP"

            if mode == "GEOM":
                parts = line.split()
                if len(parts) >= 11:
                    try:
                        self.station_chords[float(parts[0])] = float(parts[1])
                    except:
                        pass
            elif mode == "CP":
                if line.startswith("Strip"):
                    if current_strip:
                        self.strips_raw.append(np.array(current_strip))
                    current_strip = []
                else:
                    parts = line.split()
                    if len(parts) == 9:
                        try:
                            int(parts[0])
                            current_strip.append([float(p) for p in parts[1:]])
                        except:
                            pass
        if current_strip:
            self.strips_raw.append(np.array(current_strip))

        return self._process_results()

    def _process_results(self):
        pre_data = []
        q_dyn = 0.5 * self.rho * (self.v_inf ** 2)

        for panels in self.strips_raw:
            y_avg = np.mean(panels[:, 1])
            if y_avg < 0.01: continue

            p_te = (panels[0, 0:3] + panels[-1, 0:3]) / 2.0
            mid = len(panels) // 2
            p_le = (panels[mid - 1, 0:3] + panels[mid, 0:3]) / 2.0

            forces = -panels[:, 7, np.newaxis] * q_dyn * panels[:, 6, np.newaxis] * panels[:, 3:6]

            pre_data.append({
                'y': y_avg, 'p_le': p_le, 'p_te': p_te,
                'f_total_vec': np.sum(forces, axis=0),
                'p_coords': panels[:, 0:3], 'p_forces': forces
            })

        pre_data.sort(key=lambda x: x['y'])
        final_data = []
        num_stations = len(pre_data)

        for i in range(num_stations):
            station = pre_data[i]
            c_vec = station['p_te'] - station['p_le']
            c_u = c_vec / np.linalg.norm(c_vec)

            if i < num_stations - 1:
                s_vec = pre_data[i + 1]['p_le'] - station['p_le']
            else:
                s_vec = station['p_le'] - pre_data[i - 1]['p_le']

            s_u_raw = s_vec / np.linalg.norm(s_vec)

            # Define Normal as Span cross Chord (Points UP)
            n_u = np.cross(s_u_raw, c_u)
            n_u /= np.linalg.norm(n_u)

            # Redefine Span as Chord cross Normal (Points OUTBOARD)
            s_u = np.cross(c_u, n_u)
            s_u /= np.linalg.norm(s_u)

            closest_y_geom = min(self.station_chords.keys(), key=lambda k: abs(k - station['y']))
            chord_len = self.station_chords[closest_y_geom]

            r_c4 = station['p_le'] + 0.25 * chord_len * c_u
            r_rel = station['p_coords'] - r_c4
            m_total_vec = np.sum(np.cross(r_rel, station['p_forces']), axis=0)
            torsion = np.dot(m_total_vec, s_u)
            f_normal = np.dot(station['f_total_vec'], n_u)
            f_tangential = np.dot(station['f_total_vec'], c_u)
            f_spanwise =np.dot(station['f_total_vec'], s_u)

            final_data.append({
                'y': station['y'],
                'chord': chord_len,
                'Fx': station['f_total_vec'][0],
                'Fy': station['f_total_vec'][1],
                'Fz': station['f_total_vec'][2],
                'Fn': f_normal,
                'Ft': f_tangential,
                'Fs': f_spanwise,
                'Torsion': torsion,
                'cx': c_u[0], 'cy': c_u[1], 'cz': c_u[2],
                'sx': s_u[0], 'sy': s_u[1], 'sz': s_u[2],
                'nx': n_u[0], 'ny': n_u[1], 'nz': n_u[2],
                'p_le_x': station['p_le'][0], 'p_le_y': station['p_le'][1], 'p_le_z': station['p_le'][2],
                'r_c4_x': r_c4[0],
                'r_c4_y': r_c4[1],
                'r_c4_z': r_c4[2]
            })

        return pd.DataFrame(final_data)


if __name__ == "__main__":
    parser = XFLR5Parser("MainWing_a=5.00_v=75.00ms.txt")
    df_results = parser.parse()

    if df_results is not None:
        df_results.to_csv("xflr5_parsed.csv", index=False)
        print("Success")