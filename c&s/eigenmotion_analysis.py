import numpy as np
import matplotlib.pyplot as plt
from matrices import *
from pathlib import Path
import os
import pandas as pd



Path("plots").mkdir(exist_ok = True)

eigenvals_sym, eigenvecs_sym = np.linalg.eig(A_sym)
eigenvals_asym, eigenvecs_asym = np.linalg.eig(A_Asym)
print(eigenvals_sym)
print(eigenvals_asym)
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
    plt.close(fig)
   

def plot_disturbance_response_sym(sys,t_final, labels, title, e_angle):
    t = np.linspace(0,t_final, 5000)
    u_e = np.zeros_like(t)
    u_e[(t >= 1) & (t < 2)] = np.deg2rad(e_angle) #elevator disturbance for 1 second

    t, y = ctr.forced_response(sys, t, u_e)
    fig, axs = plt.subplots(4, 1, figsize=(8, 8), sharex=True)

    for i in range(4):
        axs[i].plot(t, y[i])
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
    plt.close(fig)

def plot_disturbance_response_asym(sys, t_final, labels, title, a_angle, r_angle):
    t = np.linspace(0,t_final, 5000)
    u_r = np.zeros_like(t)
    u_r[(t >= 1) & (t < 2)] = np.deg2rad(r_angle) #rudder disturbance for 1 second
    u_a = np.zeros_like(t)
    u_a[(t >= 1) & (t < 2)] = np.deg2rad(a_angle) #aileron disturbance for 1 second
    u_ra= np.vstack((u_a, u_r))
    t, y = ctr.forced_response(sys,t, u_ra)

    fig, axs = plt.subplots(4, 1, figsize=(8, 8), sharex=True)

    for i in range(4):
        axs[i].plot(t, y[i])
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
    plt.close(fig)



#short period
idx_sp = np.argmax(np.abs(np.imag(eigenvals_sym))) #largest imaginary eigenvalue is sp
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
idx_spiral = real_modes[np.argmax(np.real(eigenvals_asym[real_modes]))]
x0_spiral = np.real(eigenvecs_asym[:, idx_spiral])
plot_response(sys_asym, x0_spiral, ["beta","phi","p","r"], t_final = 500, title = 'spiral')

#disturbance due to control input
plot_disturbance_response_sym(sys_sym, 10, ["u", 'alpha','theta','q'], title = "short period disturbance", e_angle=2)
plot_disturbance_response_sym(sys_sym, 200, ["u","alpha","theta","q"], title = 'phugoid disturbance',e_angle=2)
plot_disturbance_response_asym(sys_asym, 40,  ["beta","phi","p","r"], title = 'dutch roll disturbance', a_angle = 0, r_angle = 2)
plot_disturbance_response_asym(sys_asym, 50,["beta","phi","p","r"], title = 'aperiodic roll disturbance', a_angle = 2, r_angle = 0)  
plot_disturbance_response_asym(sys_asym, 50, ["beta","phi","p","r"], title = 'spiral disturbance', a_angle = 2, r_angle=0)

#modal characteristics
def modal_char(eigenvals,idx):
    eigval = eigenvals[idx]
    sigma = np.real(eigval)
    omega = np.imag(eigval)
    return sigma, omega

def oscillatory(sigma, omega):
    wn = np.sqrt(sigma**2 + omega**2)
    zeta = -sigma /wn
    period = 2*np.pi /abs(omega)
    t_half = np.log(2)/(abs(sigma))
    return wn, zeta, period, t_half

def real_m(sigma):
    tau = -1/sigma
    t_half = np.log(2)/(-sigma)
    return tau, t_half

modes = {
    "short period": (eigenvals_sym, idx_sp, "oscillatory"),
    "Phugoid": (eigenvals_sym, idx_ph, "oscillatory"),
    "Dutch roll":(eigenvals_asym, idx_dr, "oscillatory"),
    "Aperiodic Roll": (eigenvals_asym, idx_ar, "real"),
    "Spiral": (eigenvals_asym, idx_spiral, "real")
}
results =[]
for eig_mode, (eigvals, idx, mode_type) in modes.items():
    sigma, omega = modal_char(eigvals, idx)
    if mode_type =="oscillatory":
        wn, zeta, period, t_half = oscillatory(sigma, omega)
        results. append({"Eigenmode": eig_mode,
            "sigma": sigma,
            "omega": omega,
            "wn": wn,
            "zeta": zeta,
            "period": period,
            "t_half": t_half,
            "tau": np.nan
        })
    else: 
        tau, t_half = real_m(sigma)
        results.append({
            "Eigenmode": eig_mode,
            "sigma": sigma,
            "omega": omega,
            "wn": np.nan,
            "zeta": np.nan,
            "period": np.nan,
            "t_half": t_half,
            "tau": tau
        })
mode_characteristics = pd.DataFrame(results)
print(mode_characteristics.round(4).to_string(index=False))

#disturbances 
state_disturbances_sym = {
    "u_dist": np.array([1,0,0,0]),
    "alpha_dist": np.array([0, np.deg2rad(1),0,0]),
    "theta_dist": np.array([0,0,np.deg2rad(1),0]),
    "q_dist": np.array([0,0,0,np.deg2rad(5)])
}
state_disturbances_asym = {
    "beta_dist": np.array([np.deg2rad(5),0,0,0]),
    "phi_dist": np.array([0,np.deg2rad(5),0,0]),
    "p_dist": np.array([0,0,np.deg2rad(10),0]),
    "r_dist": np.array([0,0,0,np.deg2rad(10)])
}
for title, x0 in state_disturbances_sym.items():
    plot_response(sys_sym,x0, ["u", "alpha", "theta", "q"], t_final=200, title=title)
for title,x0 in state_disturbances_asym.items():
    plot_response(sys_asym, x0,["beta","phi","p","r"], t_final=200, title=title )

E = CL*(Clb*Cnr-Cnb*Clr)
print(E)
polynomial = np.poly(A_Asym)
A,B,C,D,E = polynomial
R = B*C*D - A*D**2 - B**2*E

print(R)