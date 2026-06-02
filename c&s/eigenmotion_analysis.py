import numpy as np
import matplotlib.pyplot as plt
from matrices import *
from pathlib import Path



Path("plots").mkdir(exist_ok = True)

eigenvals_sym, eigenvecs_sym = np.linalg.eig(A_sym)
eigenvals_asym, eigenvecs_asym = np.linalg.eig(A_Asym)

def plot_response(sys, x0, labels, t_final, title):
    t = np.linspace(0,t_final, 5000)

    t,y = ctr.initial_response(sys, t, x0)
    
    fig, axs = plt.subplots(len(labels), 1, figsize=(8,8),sharex=True)

    for i in range(len(labels)):
        axs[i].plot(t,y[i])
        axs[i].set_ylabel(labels[i])
        axs[i].grid()

    axs[-1].set_xlabel("Time [s]")
    fig.suptitle(title)

    save_dir = os.path.join(
    os.path.dirname(__file__),
    "eigenmotion_figures")

    os.makedirs(save_dir, exist_ok=True)

    save_path = os.path.join(save_dir, f"{title.lower()}.png")

    plt.savefig(save_path, dpi=300, bbox_inches="tight")
   
    # filename = f"plots/{title.lower().replace(' ', '_')}.png"
    


#short period
idx_sp = np.argmax(np.abs(np.imag(eigenvals_sym)))
x0_sp = np.real(eigenvecs_sym[:,idx_sp])
plot_response(sys_sym, x0_sp, ["u", 'alpha','theta','q'],t_final=10, title = "short period")

#phugoid
complex_modes = np.where(np.abs(np.imag(eigenvals_sym)) > 1e-6)[0]
idx_ph = complex_modes[np.argmin(np.abs(np.imag(eigenvals_sym[complex_modes])))]
x0_ph = np.real(eigenvecs_sym[:, idx_ph])
plot_response(sys_sym, x0_ph, ["u","alpha","theta","q"], t_final=200, title = 'phugoid')

#Dutch roll
idx_dr = np.argmax(np.abs(np.imag(eigenvals_asym)))
x0_dr = np.real(eigenvecs_asym[:, idx_dr])
plot_response(sys_asym, x0_dr, ["beta","phi","p","r"], t_final = 40, title = 'dutch roll')

#Aperiodic roll
real_modes = np.where(np.abs(np.imag(eigenvals_asym)) < 1e-6 )[0]
idx_ar = real_modes[np.argmin(np.real(eigenvals_asym[real_modes]))]
x0_ar = np.real(eigenvecs_asym[:, idx_ar])
plot_response(sys_asym, x0_ar,["beta","phi","p","r"], t_final=5, title = 'aperiodic roll')

#Spiral
spiral_options = [i for i in real_modes if i != idx_ar]
idx_spiral = spiral_options[0]
x0_spiral = np.real(eigenvecs_asym[:, idx_spiral])
plot_response(sys_asym, x0_spiral, ["beta","phi","p","r"], t_final = 500, title = 'spiral')

#disturbance inputs

#velocity disturbance
x0_u = np.array([1,0,0,0]) #phugoid
x0_alpha = np.array([0, np.deg2rad(2),0,0]) #short period
x0_beta = np.array([np.deg2rad(5),0,0,0]) #dutch roll
x0_p = np.array([0,0,np.deg2rad(10),0]) #aperiodic roll
x0_phi = np.array([0,np.deg2rad(5), 0,0]) #spiral

plot_response(sys_sym, x0_u, ["u","alpha","theta","q"], t_final=200, title = 'v disturbance, phugoid')
plot_response(sys_sym, x0_alpha, ["u", 'alpha','theta','q'],t_final=10, title = 'alpha disturbance, short period')
plot_response(sys_asym, x0_beta,  ["beta","phi","p","r"], t_final = 40, title = 'beta disturbance, dutch roll')
plot_response(sys_asym, x0_p,["beta","phi","p","r"], t_final=5, title = 'roll disturbance, aperiodic roll')
plot_response(sys_asym, x0_phi, ["beta","phi","p","r"], t_final = 500, title = 'phi disturbance, spiral')

def modal_char(eigenvals,idx):
    eigval = eigenvals[idx]
    sigma = np.real(eigval)
    omega = np.imag(eigval)
    return sigma, omega

def oscillatory(sigma, omega):
    wn = np.sqrt(sigma**2 + omega**2)
    zeta = -sigma /wn
    period = 2*np.pi /abs(omega)
    t_half = np.log(2)/(-sigma)
    return wn, zeta, period, t_half

def real_m(sigma):
    tau = -1/sigma
    t_half = np.log(2)/(-sigma)
    return tau, t_half


sigma_sp, omega_sp = modal_char(eigenvals_sym, idx_sp)
wn_sp, zeta_sp, T_sp, t_half_sp = oscillatory(sigma_sp, omega_sp)

sigma_ph, omega_ph = modal_char(eigenvals_sym, idx_ph)
wn_ph, zeta_ph, T_ph, t_half_ph = oscillatory(sigma_ph, omega_ph)

sigma_dr, omega_dr = modal_char(eigenvals_asym,idx_dr)
wn_dr, zeta_dr, T_dr, t_half_dr = oscillatory(sigma_dr, omega_dr)

sigma_ar, omega_ar = modal_char(eigenvals_asym, idx_ar)
tau_ar, t_half_ar = real_m(sigma_ar)

sigma_spiral, omega_spiral = modal_char(eigenvals_asym, idx_spiral)
tau_spiral, t_half_spiral = real_m(sigma_spiral)

print("sigma", sigma_spiral)
print("omega",omega_spiral)
print("wn", wn_dr)
print("zeta", zeta_dr)
print("period", T_dr)
print("half time period", t_half_spiral)
print("tau", tau_spiral)
