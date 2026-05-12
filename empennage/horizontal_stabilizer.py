import numpy as np

# Input parameters
## Main (get these from YAML file)
X_aftcg = 1
l_fuselage = 8.8
MAC = 1.5
S_w = 23
c = 1.5
AR = 11

## Horizontal stabilizer parameters
X_h = 3
V_h = 0.786 # Average Roskam (see Excel)
AR_h = 4 
Sweep_h_LE = 25 #deg
Taper_h = 0.4
t_c_h = 0.12 # Thickness to chord ratio # AIRFOIL = 

# Calculations

## Planform geometry of the horizontal stabilizer

def horizontal_stabilizer(X_aftcg, l_fuselage, MAC, S_w, c, AR, X_h, V_h, AR_h, Sweep_h_LE, Taper_h, t_c_h):
    # Calculate the area of the horizontal stabilizer
    S_h = (V_h * S_w * c) / (X_h - X_aftcg)

    # Calculate the span of the horizontal stabilizer
    b_h = np.sqrt(AR_h * S_h)

    # Calculate the chord lengths at the root and tip of the horizontal stabilizer
    c_root = (2 * S_h) / (b_h * (1 + Taper_h))
    c_tip = Taper_h * c_root

    # Calculate the sweep angle at the quarter chord point
    Sweep_h_qc = Sweep_h_LE - np.arctan((c_root - c_tip) / b_h) * (180 / np.pi)

    return S_h, b_h, c_root, c_tip, Sweep_h_qc

def plot_horizontal_stabilizer(S_h, b_h, c_root, c_tip, Sweep_h_qc):
    # Create a simple plot of the horizontal stabilizer planform
    import matplotlib.pyplot as plt

    # Define the coordinates of the corners of the horizontal stabilizer
    x = [0, c_root, c_root + b_h * np.tan(np.radians(Sweep_h_qc)), b_h * np.tan(np.radians(Sweep_h_qc))]
    y = [0, 0, b_h, 0]

    plt.figure(figsize=(10, 5))
    plt.plot(x, y, 'b-')
    plt.fill(x, y, 'lightblue', alpha=0.5)
    plt.title('Horizontal Stabilizer Planform')
    plt.xlabel('Chordwise Distance (m)')
    plt.ylabel('Spanwise Distance (m)')
    plt.axis('equal')
    plt.grid()
    plt.show()

