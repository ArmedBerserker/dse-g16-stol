import pyvista as pv
import numpy as np
import pandas as pd
from pyvista import examples

# -------------------------
# Load Earth
# -------------------------
mesh = examples.planets.load_earth()
texture = examples.load_globe_texture()

# -------------------------
# Locations
# -------------------------
locations = pd.DataFrame({
    "City": ["Amsterdam", "New York", "Tokyo", "Sydney"],
    "Lat": [52.3676, 40.7128, 35.6762, -33.8688],
    "Lon": [4.9041, -74.0060, 139.6503, 151.2093]
})

def latlon_to_xyz(lat, lon, r=1.00):
    lat = np.radians(lat)
    lon = np.radians(lon + 180)

    x = r * np.cos(lat) * np.cos(lon)
    y = r * np.cos(lat) * np.sin(lon)
    z = r * np.sin(lat)

    return x, y, z

plotter = pv.Plotter(off_screen=True)

# Earth
actor = plotter.add_mesh(mesh, texture=texture, smooth_shading=True)

# -------------------------
# Add city markers
# -------------------------
for _, row in locations.iterrows():
    x, y, z = latlon_to_xyz(row["Lat"], row["Lon"])

    # red sphere marker
    marker = pv.Sphere(radius=0.001, center=(x, y, z))
    plotter.add_mesh(marker, color="red")

    # label
    plotter.add_point_labels(
        [(x, y, z)],
        [row["City"]],
        font_size=12,
        text_color="white",
        point_size=0
    )

# -------------------------
# Show
# -------------------------
plotter.open_gif("rotator")
pts =