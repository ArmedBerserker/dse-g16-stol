from math import radians, sin, cos, sqrt, atan2
import numpy as np
import pandas as pd
from sklearn.neighbors import BallTree
from tqdm import tqdm
import matplotlib.pyplot as plt

def distance(lat1, lon1, lat2, lon2):
    R = 3440.065  # Radius of Earth in nm

    lat1, lon1 = radians(lat1), radians(lon1)
    lat2, lon2 = radians(lat2), radians(lon2)

    dlat = lat2 - lat1
    dlon = lon2 - lon1

    a = sin(dlat/2)**2 + cos(lat1)*cos(lat2)*sin(dlon/2)**2
    c = 2 * atan2(sqrt(a), sqrt(1-a))

    return R * c * 1.852

def filter_airports(csv1, csv2):
    airports1 = pd.read_csv(csv1)
    airports2 = pd.read_csv(csv2)
    airports = airports1.merge(
        airports2,
        on="ICAO",
        how="inner",
        suffixes=("_airport", "_runway")
    )
    # Check correct number of airports
    common = set(airports1["ICAO"]) & set(airports2["ICAO"])
    if len(common) != len(airports):
        raise ValueError(f'Merged dataframe has {len(airports)} airports, expecting {len(common)}')
    
    # Check if we need margins
    airports_left = airports[
        (airports["le_elevation_ft"] <= 2000) |
        (airports["le_elevation_ft"].isna())
    ]
    airports_left = airports_left[
        (airports_left["length_ft"].notna()) &
        (airports_left["length_ft"] <= 200)
    ]
    airports_left = airports_left[
        (airports_left["longitude_deg"].notna())
    ]
    airports_left = airports_left[
        (airports_left["latitude_deg"].notna())
    ]
    types_to_remove = [
        "heliport",
        "seaplane_base",
        "balloonport",
        "closed"
    ]
    airports_left = airports_left[
        ~airports_left["type"].isin(types_to_remove)
    ]
    surface_types_to_remove = ...
    airports_left = airports_left[
        ~airports_left["surface"].isin(surface_types_to_remove)
    ]
    return airports_left

def reachable_airports(airport_idx,
                       max_range_km,
                       tree):
    EARTH_RADIUS = 6371.0  # km
    radius = max_range_km / EARTH_RADIUS

    id = tree.query_radius(
        coords[airport_idx:airport_idx+1],
        r=radius
    )[0]

    return airports.iloc[id]

if __name__ == "__main__":

    file = ...
    airports = filter_airports(file)

    EARTH_RADIUS = 6371.0  # km

    max_range_km = 4000

    # Convert coordinates to radians
    coords = np.radians(
        airports[['latitude', 'longitude']].values
    )

    # Build BallTree
    tree = BallTree(coords, metric='haversine')

    # Angular radius
    radius = max_range_km / EARTH_RADIUS

    routes = []

    for i in tqdm(range(len(airports))):

        neighbors = tree.query_radius(
            coords[i:i+1],
            r=radius,
            return_distance=True
        )

        idxs = neighbors[0][0]
        dists = neighbors[1][0] * EARTH_RADIUS

        for j, dist in zip(idxs, dists):

            if i == j:
                continue

            routes.append([
                airports.iloc[i]['ICAO'],
                airports.iloc[j]['ICAO'],
                dist
            ])

    routes = pd.DataFrame(
        routes,
        columns=[
            'origin',
            'destination',
            'distance_km'
        ]
    )

    print(routes.head())

    # lookup dictionaries:
    icao_to_lat = dict(
        zip(
            airports['ICAO'],
            airports['latitude']
        )
    )

    icao_to_lon = dict(
        zip(
            airports['ICAO'],
            airports['longitude']
        )
    )

    # Plot routes from one airport:
    airport = 'EHAM'

    subset = routes[
        routes.origin == airport
    ]

    plt.figure(figsize=(12,8))

    for _, row in subset.iterrows():

        plt.plot(
            [
                icao_to_lon[row.origin],
                icao_to_lon[row.destination]
            ],
            [
                icao_to_lat[row.origin],
                icao_to_lat[row.destination]
            ],
            alpha=0.15
        )

    plt.scatter(
        airports.longitude,
        airports.latitude,
        s=1
    )

    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title(f'Reachable airports from {airport}')
    plt.show()

    # route search for one airport
    tree = BallTree(coords, metric='haversine')

    reachable_airports(airport_idx=232758,
                       max_range_km=500,
                       tree=tree)
    