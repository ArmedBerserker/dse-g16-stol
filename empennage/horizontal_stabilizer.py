import numpy as np


class HorizontalTail:

    def __init__(self, config):

        self.config = config

        # Wing data
        self.S_w = config["wing"]["geometry"]["area"]
        self.MAC = config["wing"]["geometry"]["mac"]

        # CG
        self.x_cg_aft = config["cg"]["x_cg_aft"]

        # HT inputs
        ht = config["empennage"]["horizontal_tail"]

        self.V_h = ht["volume_coefficient"]

        self.AR_h = ht["geometry"]["aspect_ratio"]
        self.taper_h = ht["geometry"]["taper_ratio"]
        self.sweep_h = ht["geometry"]["sweep_deg"]

        self.x_ac_h = ht["position"]["x_ac"]

    # --------------------------------------------------
    # Tail arm
    # --------------------------------------------------

    def tail_arm(self):

        return self.x_ac_h - self.x_cg_aft

    # --------------------------------------------------
    # Tail area
    # --------------------------------------------------

    def area(self):

        l_h = self.tail_arm()

        S_h = (
            self.V_h
            * self.S_w
            * self.MAC
            / l_h
        )

        return S_h

    # --------------------------------------------------
    # Geometry
    # --------------------------------------------------

    def geometry(self):

        S_h = self.area()

        b_h = np.sqrt(self.AR_h * S_h)

        c_root = (
            2 * S_h
            / ((1 + self.taper_h) * b_h)
        )

        c_tip = self.taper_h * c_root

        MAC_h = (
            (2/3)
            * c_root
            * (
                (1 + self.taper_h + self.taper_h**2)
                / (1 + self.taper_h)
            )
        )

        return {
            "S_h": S_h,
            "b_h": b_h,
            "c_root_h": c_root,
            "c_tip_h": c_tip,
            "MAC_h": MAC_h
        }

    # --------------------------------------------------
    # Summary
    # --------------------------------------------------

    def summary(self):

        geom = self.geometry()

        print("\n--- Horizontal Tail ---")

        for key, value in geom.items():
            print(f"{key}: {value:.3f}")