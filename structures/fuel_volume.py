import numpy as np
from scipy.integrate import simpson, cumulative_trapezoid
from scipy.interpolate import interp1d

def wingbox_volume(
    c_root,
    c_tip,
    y_tip,
    xfs_pct,
    xrs_pct,
    hfs_norm,
    hrs_norm,
    n=500
):

    y = np.linspace(0.0, y_tip, n)

    chord = c_root + (c_tip - c_root) * (y / y_tip)

    width = (xrs_pct - xfs_pct) * chord

    h_front = hfs_norm * chord
    h_rear  = hrs_norm * chord

    area = 0.5 * (h_front + h_rear) * width

    volume = simpson(area, y)

    return volume

def wingbox_volume_location(
    target_volume,
    c_root,
    c_tip,
    y_tip,
    xfs_pct,
    xrs_pct,
    hfs_norm,
    hrs_norm,
    y_fuel_start=0.0,
    n=2000
):

    # full geometric span
    y = np.linspace(0.0, y_tip, n)

    # chord from aerodynamic roo
    chord = c_root + (c_tip - c_root) * (y / y_tip)

    width = (xrs_pct - xfs_pct) * chord
    h_front = hfs_norm * chord
    h_rear  = hrs_norm * chord

    area = 0.5 * (h_front + h_rear) * width

    # cumulative volume from aerodynamic root
    cumulative_volume = cumulative_trapezoid(area, y, initial=0)

    # volume at fuel start station
    v0 = np.interp(y_fuel_start, y, cumulative_volume)

    # slice arrays from fuel start onward
    mask = y >= y_fuel_start
    y_f = y[mask]
    v_f = cumulative_volume[mask] - v0

    total_fuel_volume = v_f[-1]

    if target_volume > total_fuel_volume:
        raise ValueError("Target exceeds available wingbox volume in fuel region")

    y_at_target = interp1d(v_f, y_f)(target_volume)

    return float(y_at_target)
