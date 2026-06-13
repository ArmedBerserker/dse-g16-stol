from math import radians, sin, cos, sqrt, atan2
import numpy as np
import pandas as pd
from sklearn.neighbors import BallTree
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
import cartopy.feature as cfeature


# ---------------------------------------------------------------------------
# Regions: edit name, lat/lon bounds and colour as needed
# Each entry: (name, lat_min, lat_max, lon_min, lon_max, colour)
# ---------------------------------------------------------------------------
REGIONS = [
    ("Western Europe",   35,  72,  -10,  30,  "#4e79a7"),
    ("Eastern Europe",   35,  72,   30,  60,  "#f28e2b"),
    ("North America",    15,  75, -170, -50,  "#e15759"),
    ("South America",   -60,  15,  -82, -34,  "#76b7b2"),
    ("Africa",          -35,  38,  -18,  52,  "#59a14f"),
    ("Middle East",      15,  42,   25,  65,  "#edc948"),
    ("South Asia",        5,  38,   60,  90,  "#b07aa1"),
    ("Southeast Asia",  -10,  28,   90, 145,  "#ff9da7"),
    ("East Asia",        20,  55,  100, 145,  "#9c755f"),
    ("Oceania",         -50,  -5,  110, 180,  "#bab0ac"),
]


def assign_region(lat, lon):
    """Return the first matching region name, or 'Other'."""
    for name, lat_min, lat_max, lon_min, lon_max, _ in REGIONS:
        if lat_min <= lat <= lat_max and lon_min <= lon <= lon_max:
            return name
    return "Other"


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
        left_on="ident",
        right_on="airport_ident",
        how="inner",
        suffixes=("_airport", "_runway")
    )
    common = set(airports1["ident"]) & set(airports2["airport_ident"])
    matched1 = airports["ident"].nunique()
    if matched1 != len(common):
        raise ValueError(
            f"Matched {matched1} unique airports, expected {len(common)}"
        )

    print(f'Loaded data for {len(airports)} airports')

    airports_left = airports[
        (airports["le_elevation_ft"] <= 2000) |
        (airports["le_elevation_ft"].isna())
    ]
    print(f'{len(airports_left)} / {len(airports)} left after max elevation selection')

    airports_left = airports_left[
        (airports_left["length_ft"].notna()) &
        (airports_left["length_ft"] >= 200)
    ]
    print(f'{len(airports_left)} / {len(airports)} left after runway length selection')

    airports_left = airports_left[airports_left["longitude_deg"].notna()]
    print(f'{len(airports_left)} / {len(airports)} left after longitude selection')

    airports_left = airports_left[airports_left["latitude_deg"].notna()]
    print(f'{len(airports_left)} / {len(airports)} left after latitude selection')

    types_to_remove = ["heliport", "seaplane_base", "balloonport", "closed"]
    airports_left = airports_left[~airports_left["type"].isin(types_to_remove)]
    print(f'{len(airports_left)} / {len(airports)} left after heliport elimination')

    airports_left["surface"] = (
        airports_left["surface"]
        .fillna("").astype(str).str.strip().str.lower()
    )
    surface_types_to_remove = [
        "ALUM", "ALUM-DECK", "ALUMINUM", "Aluminum rooftop", "CLOSED",
        "Deck", "Delete duplicate", "Delete", "G", "GG", "GRASS", "graas",
        "Gr", "GRA", "Gra", "GOOD GRASS", "Ice", "Ice - frozen lake", "L",
        "lakebed", "Snow", "SNO", "Snow/Ice", "U", "UNK", "Unknown", "Winter snow"
    ]
    word_lst = ["OIL", "Oil", "GRASS", "Grass", 'grass', "GRAAS", "Graas", 'graas']
    pattern = "|".join(w.lower() for w in word_lst)

    airports_left = airports_left[
        ~airports_left["surface"].isin([s.lower() for s in surface_types_to_remove])
        & (airports_left["surface"] != "")
        & ~airports_left["surface"].str.contains(pattern, regex=True, na=False)
    ]
    print(f'{len(airports_left)} / {len(airports)} left after surface selection')
    return airports_left


if __name__ == "__main__":

    file1 = 'Airport analysis/airports.csv'
    file2 = 'Airport analysis/runways.csv'
    airports = filter_airports(file1, file2)
    n_airports = len(airports)
    print(f'\nNumber of airports after filtering = {n_airports}')

    EARTH_RADIUS = 6371.0  # km
    max_range_km = 600

    coords = np.radians(airports[['latitude_deg', 'longitude_deg']].values)
    tree = BallTree(coords, metric='haversine')
    radius = max_range_km / EARTH_RADIUS

    # routes = []
    # for i in tqdm(range(len(airports))):
    #     neighbors = tree.query_radius(
    #         coords[i:i+1], r=radius, return_distance=True
    #     )
    #     idxs = neighbors[0][0]
    #     dists = neighbors[1][0] * EARTH_RADIUS

    #     for j, dist in zip(idxs, dists):
    #         if i >= j:  # unique pairs only — no A→B + B→A duplicates
    #             continue
    #         routes.append([
    #             airports.iloc[i]['ident'],
    #             airports.iloc[j]['ident'],
    #             dist
    #         ])

    # routes = pd.DataFrame(routes, columns=['origin', 'destination', 'distance_km'])

    routes = pd.read_csv('Airport analysis/routes.csv')

    # ------------------------------------------------------------------
    # Remove airports with no routes (from both airports df and routes csv)
    # ------------------------------------------------------------------
    connected_idents = set(routes['origin']) | set(routes['destination'])
    airports = airports[airports['ident'].isin(connected_idents)].copy()
    print(f'\nAirports with at least one route: {len(airports)}')
    print(f'Airports removed (no routes): {n_airports - len(airports)}')

    routes.to_csv('Airport analysis/routes.csv', index=False)
    print(routes.head())

    # ------------------------------------------------------------------
    # Assign regions
    # ------------------------------------------------------------------
    airports['region'] = airports.apply(
        lambda r: assign_region(r['latitude_deg'], r['longitude_deg']), axis=1
    )

    icao_to_lat = dict(zip(airports['ident'], airports['latitude_deg']))
    icao_to_lon = dict(zip(airports['ident'], airports['longitude_deg']))
    icao_to_region = dict(zip(airports['ident'], airports['region']))

    # ------------------------------------------------------------------
    # Plot 1 — Routes from a single airport on a Cartopy map
    # ------------------------------------------------------------------
    focus_airport = 'EHAM'

    subset = routes[routes.origin == focus_airport]

    fig1, ax1 = plt.subplots(
        figsize=(14, 9),
        subplot_kw={'projection': ccrs.PlateCarree()}
    )
    ax1.set_global()
    ax1.add_feature(cfeature.OCEAN, facecolor='#d6eaf8', zorder=0)
    ax1.add_feature(cfeature.LAND, facecolor='#f0f0e8', zorder=0)
    ax1.add_feature(cfeature.COASTLINE, linewidth=0.4, zorder=1)
    ax1.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle=':', zorder=1)

    for _, row in subset.iterrows():
        ax1.plot(
            [icao_to_lon[row.origin], icao_to_lon[row.destination]],
            [icao_to_lat[row.origin], icao_to_lat[row.destination]],
            color='steelblue', alpha=0.15, linewidth=0.6,
            transform=ccrs.PlateCarree(), zorder=2
        )

    ax1.scatter(
        airports.longitude_deg, airports.latitude_deg,
        s=1, color='grey', alpha=0.4,
        transform=ccrs.PlateCarree(), zorder=3
    )
    ax1.scatter(
        [icao_to_lon[focus_airport]], [icao_to_lat[focus_airport]],
        s=40, color='red', zorder=4,
        transform=ccrs.PlateCarree()
    )

    ax1.set_title(f'Reachable airports from {focus_airport} (≤{max_range_km} km)',
                  fontsize=13)
    plt.tight_layout()
    plt.savefig('Airport analysis/routes_map.png', dpi=150)
    plt.show()

    # ------------------------------------------------------------------
    # Plot 2 — Bubble chart: % of airports and % of routes per region
    # ------------------------------------------------------------------

    # Count airports per region
    airport_counts = airports['region'].value_counts()
    total_airports = len(airports)

    # Count routes per region pair — a route belongs to a region if
    # *either* endpoint is in that region (count each route once per region)
    route_regions = pd.Series(
        list(routes['origin'].map(icao_to_region)) +
        list(routes['destination'].map(icao_to_region))
    )
    route_counts = route_regions.value_counts()
    total_route_endpoints = len(route_regions)

    all_region_names = [r[0] for r in REGIONS] + ['Other']
    region_colors = {r[0]: r[5] for r in REGIONS}
    region_colors['Other'] = '#cccccc'

    # Compute bubble centre as mean lat/lon of airports in each region
    region_centers = airports.groupby('region')[['latitude_deg', 'longitude_deg']].mean()

    fig2, ax2 = plt.subplots(
        figsize=(16, 9),
        subplot_kw={'projection': ccrs.PlateCarree()}
    )
    ax2.set_global()
    ax2.add_feature(cfeature.OCEAN, facecolor='#d6eaf8', zorder=0)
    ax2.add_feature(cfeature.LAND, facecolor='#f0f0e8', zorder=0)
    ax2.add_feature(cfeature.COASTLINE, linewidth=0.4, zorder=1)
    ax2.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle=':', zorder=1)

    # All airports as small background dots
    ax2.scatter(
        airports.longitude_deg, airports.latitude_deg,
        s=1, color='grey', alpha=0.3,
        transform=ccrs.PlateCarree(), zorder=2
    )

    SCALE = 800  # tweak to resize bubbles

    for region in all_region_names:
        if region not in region_centers.index:
            continue

        pct_airports = airport_counts.get(region, 0) / total_airports * 100
        pct_routes   = route_counts.get(region, 0) / total_route_endpoints * 100
        clat = region_centers.loc[region, 'latitude_deg']
        clon = region_centers.loc[region, 'longitude_deg']
        color = region_colors.get(region, '#cccccc')

        # Outer bubble = % of routes
        ax2.scatter(
            clon, clat,
            s=pct_routes * SCALE,
            color=color, alpha=0.35,
            transform=ccrs.PlateCarree(), zorder=3
        )
        # Inner bubble = % of airports
        ax2.scatter(
            clon, clat,
            s=pct_airports * SCALE,
            color=color, alpha=0.75,
            edgecolors='white', linewidths=0.8,
            transform=ccrs.PlateCarree(), zorder=4
        )
        # Label
        ax2.text(
            clon, clat,
            f'{region}\n{pct_airports:.1f}% apt\n{pct_routes:.1f}% rte',
            fontsize=6.5, ha='center', va='center',
            transform=ccrs.PlateCarree(), zorder=5,
            color='black'
        )

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor='grey',  alpha=0.75, label='Inner bubble = % airports'),
        mpatches.Patch(facecolor='grey',  alpha=0.35, label='Outer bubble = % routes'),
    ]
    ax2.legend(handles=legend_elements, loc='lower left', fontsize=8)
    ax2.set_title('Share of airports and routes by region', fontsize=13)

    plt.tight_layout()
    plt.savefig('Airport analysis/bubble_map.png', dpi=150)
    plt.show()

# from math import radians, sin, cos, sqrt, atan2
# import numpy as np
# import pandas as pd
# from sklearn.neighbors import BallTree
# from tqdm import tqdm
# import matplotlib.pyplot as plt

# def distance(lat1, lon1, lat2, lon2):
#     R = 3440.065  # Radius of Earth in nm

#     lat1, lon1 = radians(lat1), radians(lon1)
#     lat2, lon2 = radians(lat2), radians(lon2)

#     dlat = lat2 - lat1
#     dlon = lon2 - lon1

#     a = sin(dlat/2)**2 + cos(lat1)*cos(lat2)*sin(dlon/2)**2
#     c = 2 * atan2(sqrt(a), sqrt(1-a))

#     return R * c * 1.852

# def filter_airports(csv1, csv2):
#     airports1 = pd.read_csv(csv1)
#     airports2 = pd.read_csv(csv2)
#     airports = airports1.merge(
#         airports2,
#         left_on="ident",
#         right_on="airport_ident",
#         how="inner",
#         suffixes=("_airport", "_runway")
#     )
#     common = (
#         set(airports1["ident"])
#         & set(airports2["airport_ident"])
#     )

#     # if len(common) != len(airports):
#     #     raise ValueError(
#     #         f"Merged dataframe has {len(airports)} airports, "
#     #         f"expecting {len(common)}"
#     #     )
#     matched1 = airports["ident"].nunique()

#     if matched1 != len(common):
#         raise ValueError(
#             f"Matched {matched1} unique airports, expected {len(common)}"
#         )
#     # airports1 = pd.read_csv(csv1)
#     # airports2 = pd.read_csv(csv2)
#     # airports = airports1.merge(
#     #     airports2,
#     #     on="ICAO",
#     #     how="inner",
#     #     suffixes=("_airport", "_runway")
#     # )
#     # # Check correct number of airports
#     # common = set(airports1["ICAO"]) & set(airports2["ICAO"])
#     # if len(common) != len(airports):
#     #     raise ValueError(f'Merged dataframe has {len(airports)} airports, expecting {len(common)}')
    
#     # Check if we need margins
#     print(f'Loaded data for {len(airports)} airports')
#     airports_left = airports[
#         (airports["le_elevation_ft"] <= 2000) |
#         (airports["le_elevation_ft"].isna())
#     ]
#     print(f'{len(airports_left)} / {len(airports)} left after max elevation selection')
#     airports_left = airports_left[
#         (airports_left["length_ft"].notna()) &
#         (airports_left["length_ft"] >= 200)
#     ]
#     print(f'{len(airports_left)} / {len(airports)} left after runway length selection')
#     airports_left = airports_left[
#         (airports_left["longitude_deg"].notna())
#     ]
#     print(f'{len(airports_left)} / {len(airports)} left after longitude selection')
#     airports_left = airports_left[
#         (airports_left["latitude_deg"].notna())
#     ]
#     print(f'{len(airports_left)} / {len(airports)} left after latitude selection')
#     types_to_remove = [
#         "heliport",
#         "seaplane_base",
#         "balloonport",
#         "closed"
#     ]
#     airports_left = airports_left[
#         ~airports_left["type"].isin(types_to_remove)
#     ]
#     print(f'{len(airports_left)} / {len(airports)} left after heliport elimation')
#     airports_left["surface"] = (
#         airports_left["surface"]
#         .fillna("")
#         .astype(str)
#         .str.strip()
#         .str.lower()
#     )
#     surface_types_to_remove = [
#         "ALUM", "ALUM-DECK", "ALUMINUM", "Aluminum rooftop", "CLOSED", 
#         "Deck", "Delete duplicate", "Delete", "G", "GG", "GRASS", "graas", 
#         "Gr", "GRA", "Gra", "GOOD GRASS", "Ice", "Ice - frozen lake", "L", 
#         "lakebed", "Snow", "SNO", "Snow/Ice", "U", "UNK", "Unknown", "Winter snow"
#     ]
#     word_lst = ["OIL", "Oil", "GRASS", "Grass", 'grass', "GRAAS", "Graas", 'graas']
#     pattern = "|".join(w.lower() for w in word_lst)

#     airports_left = airports_left[
#         ~airports_left["surface"].isin(
#             [s.lower() for s in surface_types_to_remove]
#         )
#         & (airports_left["surface"] != "")
#         & ~airports_left["surface"].str.contains(pattern, regex=True, na=False)
#     ]
#     print(f'{len(airports_left)} / {len(airports)} left after surface selection')
#     return airports_left

# def reachable_airports(airport_idx,
#                        max_range_km,
#                        tree):
#     EARTH_RADIUS = 6371.0  # km
#     radius = max_range_km / EARTH_RADIUS

#     id = tree.query_radius(
#         coords[airport_idx:airport_idx+1],
#         r=radius
#     )[0]

#     return airports.iloc[id]

# if __name__ == "__main__":

#     file1 = 'Airport analysis/airports.csv'
#     file2 = 'Airport analysis/runways.csv'
#     airports = filter_airports(file1, file2)
#     # airports = airports[:20]
#     n_airports = len(airports)
#     print(f'\n Number of airports that we can land at = {n_airports}')

#     EARTH_RADIUS = 6371.0  # km

#     max_range_km = 600

#     # Convert coordinates to radians
#     coords = np.radians(
#         airports[['latitude_deg', 'longitude_deg']].values
#     )

#     # Build BallTree
#     tree = BallTree(coords, metric='haversine')

#     # Angular radius
#     radius = max_range_km / EARTH_RADIUS

#     routes = []

#     for i in tqdm(range(len(airports))):

#         neighbors = tree.query_radius(
#             coords[i:i+1],
#             r=radius,
#             return_distance=True
#         )

#         idxs = neighbors[0][0]
#         dists = neighbors[1][0] * EARTH_RADIUS

#         for j, dist in zip(idxs, dists):
#             if i >= j:  # skip self and already-seen reverse pairs
#                 continue

#             if i == j:
#                 continue

#             routes.append([
#                 airports.iloc[i]['ident'],
#                 airports.iloc[j]['ident'],
#                 dist
#             ])

#     routes = pd.DataFrame(
#         routes,
#         columns=[
#             'origin',
#             'destination',
#             'distance_km'
#         ]
#     )

#     routes.to_csv('Airport analysis/routes.csv')

#     print(routes.head())

#     routes = pd.read_csv('Airport analysis/routes.csv')

#     # lookup dictionaries:
#     icao_to_lat = dict(
#         zip(
#             airports['ident'],
#             airports['latitude_deg']
#         )
#     )

#     icao_to_lon = dict(
#         zip(
#             airports['ident'],
#             airports['longitude_deg']
#         )
#     )

#     # Plot routes from one airport:
#     airport = 'EHAM'

#     subset = routes[
#         routes.origin == airport
#     ]

#     plt.figure(figsize=(12,8))

#     for _, row in subset.iterrows():

#         plt.plot(
#             [
#                 icao_to_lon[row.origin],
#                 icao_to_lon[row.destination]
#             ],
#             [
#                 icao_to_lat[row.origin],
#                 icao_to_lat[row.destination]
#             ],
#             alpha=0.15
#         )

#     plt.scatter(
#         airports.longitude_deg,
#         airports.latitude_deg,
#         s=1
#     )

#     plt.xlabel('Longitude [deg]')
#     plt.ylabel('Latitude [deg]')
#     plt.title(f'Reachable airports from {airport}')
#     plt.show()

#     # # route search for one airport
#     # tree = BallTree(coords, metric='haversine')

#     # reachable_airports(airport_idx=232758,
#     #                    max_range_km=600,
#     #                    tree=tree)
    