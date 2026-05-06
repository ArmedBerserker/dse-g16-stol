import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
df = pd.read_csv('lookups/ref.csv')

OEW_list = df['OEW [kg]'].to_list()
MTOW_list = df['MTOW [kg]'].to_list()

x = np.array(MTOW_list)
y = np.array(OEW_list)

# Remove invalid values for log-log fit
valid = (x > 0) & (y > 0)
x_log = x[valid]
y_log = y[valid]

# -------------------------
# Linear fit: OEW = a MTOW + b
# -------------------------
a_lin, b_lin = np.polyfit(x, y, 1)

y_pred_lin = a_lin * x + b_lin

ss_res_lin = np.sum((y - y_pred_lin) ** 2)
ss_tot_lin = np.sum((y - np.mean(y)) ** 2)
r2_lin = 1 - ss_res_lin / ss_tot_lin

x_line = np.linspace(0, max(x), 100)
y_line_lin = a_lin * x_line + b_lin

# -------------------------
# Log-log fit: OEW = c MTOW^d
# -------------------------
log_x = np.log10(x_log)
log_y = np.log10(y_log)

d_log, log_c = np.polyfit(log_x, log_y, 1)
c_log = 10 ** log_c

log_y_pred = d_log * log_x + log_c

ss_res_log = np.sum((log_y - log_y_pred) ** 2)
ss_tot_log = np.sum((log_y - np.mean(log_y)) ** 2)
r2_log = 1 - ss_res_log / ss_tot_log

x_line_log = np.logspace(np.log10(min(x_log)), np.log10(max(x_log)), 100)
y_line_log = c_log * x_line_log ** d_log

# -------------------------
# Print results
# -------------------------
print("Linear fit:")
print(f"OEW = {a_lin:.6f} MTOW + {b_lin:.6f}")
print(f"R² = {r2_lin:.6f}")

print("\nLog-log fit:")
print(f"OEW = {c_log:.6f} MTOW^{d_log:.6f}")
print(f"R² = {r2_log:.6f}")

# -------------------------
# Subplots
# -------------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Linear subplot
axes[0].scatter(x, y, label='Aircraft data', )
axes[0].plot(
    x_line,
    y_line_lin,
    label=f'OEW = {a_lin:.3f} MTOW + {b_lin:.1f}\n$R^2$ = {r2_lin:.3f}'
)
axes[0].set_xlabel('MTOW [kg]')
axes[0].set_ylabel('OEW [kg]')
axes[0].set_title('Linear Fit')
axes[0].set_xlim(0, max(x) * 1.05)
axes[0].set_ylim(0, max(y) * 1.05)
axes[0].grid(True)
axes[0].legend()

# Log-log subplot
axes[1].scatter(x_log, y_log, label='Aircraft data')
axes[1].plot(
    x_line_log,
    y_line_log,
    label=f'OEW = {c_log:.3f} MTOW^{d_log:.3f}\n$R^2$ = {r2_log:.3f}'
)
axes[1].set_xscale('log')
axes[1].set_yscale('log')
axes[1].set_xlabel('MTOW [kg]')
axes[1].set_ylabel('OEW [kg]')
axes[1].set_title('Log-Log Fit')
axes[1].grid(True, which='both')
axes[1].legend()

plt.tight_layout()
plt.show()
