from nvm import *



if __name__ == "__main__":
    dat_file = "NACA23012.dat"
    chord_length = 1.87
    thickness_skin = 0.001
    thickness_web = 0.004
    area_spar_cap = 0.0005
    area_stringer = 0.0001

    # Define External Weights & Point Placements
    engine_spanwise_pos = 2.2 # Engine location along span (m)
    engine_mass_kg = 92.0  #  mass in kg
    wtd_mass_kg = 20.0  # Wing tip device mass in kg

    y_steps, shear, moment, torsion = calculate_and_plot_vmt(
        # Aerodynamic Arguments
        half_wing_lift_target=10000,
        b=16.8,
        y_root=0.725,
        c=chord_length,
        q=571.0,
        flap_start=1.1,
        flap_end=3.8,
        delta_cl_flap=1.25,
        aileron_start=5.6,
        aileron_end=7.5,
        delta_cl_aileron=0.9,

        # Structural Layout Arguments
        filepath=dat_file,
        front_spar_pct=0.16,
        main_spar_pct=0.56,
        num_box_stringers=8,
        num_le_booms=4,
        t_skin=thickness_skin,
        t_web=thickness_web,
        A_spar=area_spar_cap,
        A_stringer=area_stringer,

        # Auxiliary Point Load Configurations
        engine_mass=engine_mass_kg,
        engine_y_loc=engine_spanwise_pos,
        engine_x_loc=0.30 * chord_length,  # e.g., Engine CG at 30% local chord
        wtd_mass=wtd_mass_kg,
        wtd_x_loc=0.45 * chord_length,  # Wingtip device twist point
        n_load_factor=3.8*1.5
    )

    print(f"Root Shear Force:   {shear[0] / 1000:.2f} kN")
    print(f"Root Bending Moment:{moment[0] / 1000:.2f} kN*m")
    print(f"Root Torsion {torsion[0]:.2f} N*m")