from scipy.optimize import differential_evolution, NonlinearConstraint
from size_wingbox import *

# --- Setup & Constants ---
S, AR, taper = 25.675, 9.0, 0.4
xflr5_file, airfoil_file = "MainWing_a=3.50_v=72.00ms.txt", "onze_airfoil.dat"
load_factor, y_eng, m_eng, thrust_eng = (3.8, 2.35, 112, 825)
Fn_eng, Ft_eng, dz_eng, dx_eng = m_eng * 9.81, -thrust_eng, -0.375, 0.8

span, y_tip, c_root, c_tip, _ = calculate_wing_geometry(S, AR, taper)
y_stations, T_c4, Fn_aero, Ft_aero, _ = load_aerodynamic_data(xflr5_file)
chords = c_root - (c_root - c_tip) * (y_stations / y_tip)
b_limit_val = span / 2 * 1000
N_STR_FIXED = 16

MAT_CATALOG = [
    {'name': 'AA6061-T6', 'E': 68.9e9, 'G': 26.0e9, 'rho': 2700},
    {'name': 'AA2024-T3', 'E': 73.1e9, 'G': 28.0e9, 'rho': 2780}
    #{'name': 'AA2099-T83 (Al-Li)',    'E': 75.0e9, 'G': 28.5e9, 'rho': 2630},
    #{'name': 'GFRP (E-Glass/Epoxy)',  'E': 22.0e9, 'G':  8.5e9, 'rho': 1950},
    #{'name': 'Sitka Spruce',          'E': 11.2e9, 'G':  0.75e9, 'rho': 400}
]


def compute_all(x):
    # 1. Unpack the expanded design vector
    xfs, xrs, tskin, tspar_web, wspar_cap, tspar_cap, Astr, idx_skin, idx_spar, idx_str = x

    # Fast fail for invalid geometry (front spar behind rear spar)
    # Return a high penalty mass and an MOS of 0 to immediately fail the constraint
    if xfs >= xrs: return 1e5, 1e5, 1e5, 1e5, 0.0

    m_skin = MAT_CATALOG[int(idx_skin)]
    m_spar = MAT_CATALOG[int(idx_spar)]
    m_str = MAT_CATALOG[int(idx_str)]

    materials = {'skin': m_skin, 'spar': m_spar, 'str': m_str}
    E_ref, G_ref = m_spar['E'], m_spar['G']

    try:
        # 2. Structural Properties & Mass
        Ixx, Izz, J, x_sc = compute_wingbox_properties(
            chords, airfoil_file, xfs, xrs, tskin,
            tspar_web, wspar_cap, tspar_cap,
            Astr, N_STR_FIXED, materials, E_ref, G_ref
        )
        half_mass, _, _, _, m_prime = compute_wingbox_mass(
            y_stations, chords, airfoil_file, xfs, xrs, tskin,
            tspar_web, wspar_cap, tspar_cap,
            Astr, N_STR_FIXED, materials
        )

        # 3. Internal Loads
        V_stations, Mx, Mz, T = compute_internal_loads(
            y_stations, chords, T_c4, Fn_aero, Ft_aero, m_prime, x_sc,
            load_factor, y_eng, Fn_eng, Ft_eng, dz_eng, dx_eng
        )

        # 4. Buckling Check
        hfs_norm, hrs_norm = get_airfoil_heights(airfoil_file, xfs, xrs)
        fspars = chords * hfs_norm
        rspars = chords * hrs_norm

        buckle_mos_front = check_buckle(fspars, tspar_web, V_stations, m_spar['E'])
        buckle_mos_rear = check_buckle(rspars, tspar_web, V_stations, m_spar['E'])
        min_mos = min(buckle_mos_front, buckle_mos_rear)

        # 5. Deflections
        twist_tip = abs(calculate_torsional_deflection(y_stations, T, G_ref, J)[-1])
        by_tip = abs(calculate_bending_deflection(y_stations, Mx, E_ref, Ixx)[-1])
        bz_tip = abs(calculate_bending_deflection(y_stations, Mz, E_ref, Izz)[-1])

        return (half_mass * 2, twist_tip, by_tip, bz_tip, min_mos)

    except Exception:
        # Return failing penalty values on math crash
        return 1e5, 1e5, 1e5, 1e5, 0.0


def objective(x):
    return compute_all(x)[0]


def constraint_functions(x):
    _, twist, b_y, b_z, min_mos = compute_all(x)
    return [b_y, b_z, twist, x[1] - x[0], min_mos]


if __name__ == '__main__':
    num_mats = len(MAT_CATALOG)

    # Order matches the unpacked variables in compute_all()
    bounds = [
        (0.10, 0.20),  # xfs
        (0.55, 0.65),  # xrs
        (0.0005, 0.0015),  # tskin
        (0.003, 0.005),  # tspar_web
        (0.010, 0.050),  # wspar_cap
        (0.002, 0.010),  # tspar_cap
        (1e-5, 5e-5),  # Astr
        (0, num_mats - 1),  # Skin Material Index
        (0, num_mats - 1),  # Spar Material Index
        (0, num_mats - 1)  # Stringer Material Index
    ]

    # Explicitly specify which indices are continuous (False) and discrete (True)
    integrality = [False, False, False, False, False, False, False, True, True, True]

    # Add the MOS bounds
    nl_constraints = NonlinearConstraint(
        constraint_functions,
        lb=[0.0, 0.0, 0.0, 0.10, 1.1],
        ub=[0.10 * b_limit_val, 0.05 * b_limit_val, 5.0, 1.0, np.inf]
    )

    res = differential_evolution(
        objective, bounds=bounds, constraints=nl_constraints,
        integrality=integrality, popsize=50, maxiter=100,
        disp=True, workers=10, polish=True
    )

    f_xfs, f_xrs, f_tskin, f_tspar_web, f_wspar_cap, f_tspar_cap, f_Astr, f_skin_idx, f_spar_idx, f_str_idx = res.x
    m, tw, by, bz, min_mos = compute_all(res.x)

    print("\n" + "=" * 50)
    print("--- OPTIMIZATION RESULTS ---")
    print(f"Minimum Mass: {m:.2f} kg")
    print("-" * 50)
    print(f"Materials:")
    print(f"  Skin Component        : {MAT_CATALOG[int(f_skin_idx)]['name']}")
    print(f"  Spar Component        : {MAT_CATALOG[int(f_spar_idx)]['name']}")
    print(f"  Stringer Component    : {MAT_CATALOG[int(f_str_idx)]['name']}")
    print("-" * 50)
    print(f"Design Layout Parameters:")
    print(f"  Front Spar (xfs)      : {f_xfs:.4f}")
    print(f"  Rear Spar (xrs)       : {f_xrs:.4f}")
    print(f"  Skin Thick. (tskin)   : {f_tskin * 1000:.3f} mm")
    print(f"  Spar Web (tweb)       : {f_tspar_web * 1000:.3f} mm")
    print(f"  Spar Cap Width (wcap) : {f_wspar_cap * 1000:.3f} mm")
    print(f"  Spar Cap Thick. (tcap): {f_tspar_cap * 1000:.3f} mm")
    print(f"  Stringer Area (Astr)  : {f_Astr:.3e} m^2")
    print("-" * 50)
    print("Deflections & Constraints:")
    print(f"  Bending y Deflection  : {by:.2f} / Limit: {0.10 * b_limit_val:.2f} mm")
    print(f"  Bending z Deflection  : {bz:.2f} / Limit: {0.05 * b_limit_val:.2f} mm")
    print(f"  Torsional Twist       : {tw:.4f}° / Limit: 5.000°")
    print(f"  Buckling Margin (MOS) : {min_mos:.3f} / Limit: > 1.1")
    print("=" * 50)