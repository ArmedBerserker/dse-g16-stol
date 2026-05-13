import numpy as np
import matplotlib.pyplot as plt

def plot_empennage(emp_data: dict):
    """
    Renders precise 2D planform top-views of the horizontal stabilizer
    and side-views of the vertical stabilizer using dictionary inputs.
    """
    h = emp_data["horizontal_tail"]
    v = emp_data["vertical_tail"]
    
    # Create a unified figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # ==========================================================
    # 1. HORIZONTAL STABILIZER (Top-Down View)
    # ==========================================================
    b_h = h["b_h"]
    c_root_h = h["c_root"]
    c_tip_h = h["c_tip"]
    sweep_qc_h = h["sweep_qc"]
    
    half_span_h = b_h / 2
    
    # Derive X-coordinates based on the quarter-chord line geometry
    x_qc_root_h = 0.25 * c_root_h
    x_qc_tip_h = x_qc_root_h + half_span_h * np.tan(np.radians(sweep_qc_h))
    
    x_le_tip_h = x_qc_tip_h - 0.25 * c_tip_h
    x_te_tip_h = x_le_tip_h + c_tip_h
    
    # Starboard (Right) half-wing polygon array
    x_h_starboard = [0.0, x_le_tip_h, x_te_tip_h, c_root_h, 0.0]
    y_h_starboard = [0.0, half_span_h, half_span_h, 0.0, 0.0]
    
    # Plot Starboard side
    ax1.plot(x_h_starboard, y_h_starboard, color='#0066CC', lw=1.5, label='Horizontal Tail')
    ax1.fill(x_h_starboard, y_h_starboard, color='#0066CC', alpha=0.2)
    
    # Mirror Port (Left) side for full planform representation
    y_h_port = [-y for y in y_h_starboard]
    ax1.plot(x_h_starboard, y_h_port, color='#0066CC', lw=1.5)
    ax1.fill(x_h_starboard, y_h_port, color='#0066CC', alpha=0.2)
    
    # Centerline chord reference
    ax1.plot([0, c_root_h], [0, 0], color='black', linestyle='-.', lw=0.8, alpha=0.5)
    
    ax1.set_title("Horizontal Stabilizer Planform (Top View)", fontsize=12, fontweight='bold')
    ax1.set_xlabel("Chordwise Axis: X [m]")
    ax1.set_ylabel("Spanwise Axis: Y [m]")
    ax1.axis("equal")
    ax1.grid(True, linestyle=":", alpha=0.6)

    # ==========================================================
    # 2. VERTICAL STABILIZER (Profile Side-View)
    # ==========================================================
    b_v = v["b_v"]
    c_root_v = v["c_root"]
    c_tip_v = v["c_tip"]
    sweep_qc_v = v["sweep_qc"]
    
    # Derive X-coordinates based on the quarter-chord line geometry
    x_qc_root_v = 0.25 * c_root_v
    x_qc_tip_v = x_qc_root_v + b_v * np.tan(np.radians(sweep_qc_v))
    
    x_le_tip_v = x_qc_tip_v - 0.25 * c_tip_v
    x_te_tip_v = x_le_tip_v + c_tip_v
    
    # Vertical fin polygon array (Z acts as the spanwise axis standing on the fuselage)
    x_v_fin = [0.0, x_le_tip_v, x_te_tip_v, c_root_v, 0.0]
    z_v_fin = [0.0, b_v, b_v, 0.0, 0.0]
    
    # Plot Vertical Fin
    ax2.plot(x_v_fin, z_v_fin, color='#CC0000', lw=1.5, label='Vertical Fin')
    ax2.fill(x_v_fin, z_v_fin, color='#CC0000', alpha=0.2)
    
    # Fuselage mounting line reference
    ax2.axhline(0, color='black', linestyle='-', lw=1.2, label='Fuselage Crown')
    
    ax2.set_title("Vertical Stabilizer Planform (Side View)", fontsize=12, fontweight='bold')
    ax2.set_xlabel("Chordwise Axis: X [m]")
    ax2.set_ylabel("Vertical Axis: Z [m]")
    ax2.axis("equal")
    ax2.grid(True, linestyle=":", alpha=0.6)
    ax2.legend(loc='upper right')

    # Final formatting display
    plt.tight_layout()
    plt.show()

# Test execution stub if run directly
if __name__ == "__main__":
    # Mock output matching the pipeline data schema
    dummy_emp_data = {
        "horizontal_tail": {"b_h": 3.43, "c_root": 1.70, "c_tip": 0.68, "sweep_qc": 21.6},
        "vertical_tail": {"b_v": 2.53, "c_root": 2.40, "c_tip": 0.96, "sweep_qc": 20.3}
    }
    plot_empennage(dummy_emp_data)