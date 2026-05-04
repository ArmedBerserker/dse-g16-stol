import numpy as np

class CONSTS:
    g0 = 9.80665        # gravitational acceleration [m/s^2]
    R = 287.05287      # specific gas constant for air [J/(kg·K)]
    T0 = 288.15        # sea-level standard temperature [K]
    P0 = 101325.0      # sea-level standard pressure [Pa]
    rho0 = 1.225       # sea-level standard density [kg/m^3]
    L = -0.0065        # temperature lapse rate [K/m]
    r = 6356766.0      # nominal Earth radius [m], if using geopotential height


class Atmosphere:
    def __init__(self, height, delta_T = 0.):
        self.height = np.asarray(height)

        if self.height.ndim == 0:
            self.height = self.height[None]

        if np.any(self.height < 0) or np.any(self.height > 10000):
            raise ValueError('Height must be between 0m and 10,000m.')
        
        self.delta_T = delta_T

        self.H = CONSTS.r * self.height / (CONSTS.r + self.height)

        self.temp_isa = CONSTS.T0 + CONSTS.L * self.H

        self.temp = self.temp_isa + self.delta_T

        self.pressure = CONSTS.P0 * (self.temp_isa / CONSTS.T0) ** (
            -CONSTS.g0 / (CONSTS.L * CONSTS.R)
        )

        self.density = self.pressure / (CONSTS.R * self.temp)


        self.density_ratio = self.density / CONSTS.rho0
        self.pressure_ratio = self.pressure / CONSTS.P0

    def __str__(self):
        return(f"ISA results:\n"
        f"Altitude: {self.height} m\n"
        f"Temperature: {self.temp} K\n"
        f"Pressure: {self.pressure} Pa\n"
        f"Density: {self.density} kg/m^3\n"
        f"Density Ratio: {self.density_ratio}\n"
        f"Pressure Ratio: {self.pressure_ratio}")
    
