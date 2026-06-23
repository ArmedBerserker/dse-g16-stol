import numpy as np
import matplotlib.pyplot as plt

# Cumulative production quantity
Q = np.arange(1, 251)

# Inputs
manufacturing_hours = 284470.9839
Q_temp = 180
learning_curve_exponent = 0.95

# Learning curve calculation
manufacturing_hours_per_unit = (
    manufacturing_hours / Q_temp
) * Q ** (np.log(learning_curve_exponent) / np.log(2))

# Plot
plt.figure(figsize=(10, 6))

plt.plot(
    Q,
    manufacturing_hours_per_unit,
    linewidth=2)

# plt.title('Manufacturing Learning Curve (95% Learning Rate)', fontsize=14)
plt.xlabel('Production Unit Number Q', fontsize=12)
plt.ylabel('Manufacturing Labor Hours per Unit', fontsize=12)
plt.ylim(bottom=0)
plt.xlim(left=0)
plt.grid(True, linestyle='--', alpha=0.7)
# plt.legend()
plt.tight_layout()

plt.savefig('learning_curve.png', dpi=800)
