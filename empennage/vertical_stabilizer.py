import numpy as np


class VerticalTail:

    def __init__(self, config):

        self.config = config

        # Wing data
        self.S_w = config["wing"]["geometry"]["area"]
        self.b_w = config["wing"]["geometry"]["span"]

        # CG
        self.x_cg_aft = config["cg"]["x_cg_aft"]

        # VT inputs
        vt = config["empennage"]["vertical_tail"]

        self.V_v = vt["volume_coefficient"]

        self.AR_v = vt["geometry"]["aspect_ratio"]
        self.taper_v = vt["geometry"]["taper_ratio"]
        self.sweep_v = vt["geometry"]["sweep_deg"]

        self.x_ac_v = vt["position"]["x_ac"]

    # --------------------------------------------------
    # Tail arm
    # --------------------------------------------------

    def tail_arm(self):

        return self.x_ac_v - self.x_cg_aft

    # --------------------------------------------------
    # Tail area
    # --------------------------------------------------

    def area(self):

        l_v = self.tail_arm()

        S_v = (
            self.V_v
            * self.S_w
            * self.b_w
            / l_v
        )

        return S_v

    # --------------------------------------------------
    # Geometry
    # --------------------------------------------------

    def geometry(self):

        S_v = self.area()

        h_v = np.sqrt(self.AR_v * S_v)

        c_root = (
            2 * S_v
            / ((1 + self.taper_v) * h_v)
        )

        c_tip = self.taper_v * c_root

        MAC_v = (
            (2/3)
            * c_root
            * (
                (1 + self.taper_v + self.taper_v**2)
                / (1 + self.taper_v)
            )
        )

        return {
            "S_v": S_v,
            "h_v": h_v,
            "c_root_v": c_root,
            "c_tip_v": c_tip,
            "MAC_v": MAC_v
        }
 
    # --------------------------------------------------
    # Summary
    # --------------------------------------------------

    def summary(self):

        geom = self.geometry()

        print("\n--- Vertical Tail ---")

        for key, value in geom.items():
            print(f"{key}: {value:.3f}")