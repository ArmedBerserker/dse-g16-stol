import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import brentq

kgtoslug  = 0.00194032
ktastofps = 1.68781
mstofps   = 3.28084

rho_cruise=0.7812# kg
rho_TO= 1.0789# kg
D_ft= 5.66667# ft
P_bhp_cr=230/2# hp cruise per engine
P_bhp_to=160# hp takeoff per engine

# ── step 1: load JavaProp CSVs ─────────────────────────────────────────────
# each CSV has columns: J, Ct, Cp
# exported from JavaProp analysis mode at fixed beta

beta_angles = [10, 15, 20, 25, 30, 35, 40]   # degrees — match what you ran in JavaProp

def load_javaprop(beta_list):
    """
    Load one CSV per blade angle.
    Returns a common J grid, and Ct/Cp matrices shaped (n_beta, n_J).
    """
    all_J  = []
    all_Ct = {}
    all_Cp = {}

    for beta in beta_list:
        df = pd.read_csv(
            rf"C:\JavaProp\beta{beta}1",
            sep="\t",
            skiprows=2,
            header=0,
        )
        df.columns = df.columns.str.strip().str.lower()
        # JavaProp may use 'j', 'ct', 'cp' or similar — check your export
        J_col  = [c for c in df.columns if 'j'  in c][0]
        Ct_col = [c for c in df.columns if 'ct' in c][0]
        Cp_col = [c for c in df.columns if 'cp' in c][0]

        all_J.append(df[J_col].values)
        all_Ct[beta] = (df[J_col].values, df[Ct_col].values)
        all_Cp[beta] = (df[J_col].values, df[Cp_col].values)

    # common J grid spanning the intersection of all beta ranges
    J_min = max(j.min() for j in all_J)
    J_max = min(j.max() for j in all_J)
    J_common = np.linspace(J_min, J_max, 200)

    return J_common, all_Ct, all_Cp

# ── step 2: build 2D interpolator ─────────────────────────────────────────
def build_interpolators(beta_list, J_common, all_Ct, all_Cp):
    """
    Interpolate each beta curve onto the common J grid,
    then build a 2D spline over (beta, J).
    """
    from scipy.interpolate import CubicSpline

    Ct_matrix = np.zeros((len(beta_list), len(J_common)))
    Cp_matrix = np.zeros((len(beta_list), len(J_common)))

    for i, beta in enumerate(beta_list):
        J_b,  Ct_b = all_Ct[beta]
        J_b2, Cp_b = all_Cp[beta]

        Ct_matrix[i, :] = CubicSpline(J_b,  Ct_b)(J_common)
        Cp_matrix[i, :] = CubicSpline(J_b2, Cp_b)(J_common)

    beta_arr = np.array(beta_list, dtype=float)

    Ct_interp = RectBivariateSpline(beta_arr, J_common, Ct_matrix)
    Cp_interp = RectBivariateSpline(beta_arr, J_common, Cp_matrix)

    return Ct_interp, Cp_interp

# ── step 3: governor logic — find beta that absorbs P_required ────────────
def solve_beta(J, Cp_required, Cp_interp, beta_min, beta_max):
    """
    Find β such that Cp(β, J) = Cp_required.
    Returns nan if no solution exists in [beta_min, beta_max].
    """
    def residual(beta):
        return float(Cp_interp(beta, J)) - Cp_required

    try:
        f_min = residual(beta_min)
        f_max = residual(beta_max)
        if f_min * f_max > 0:
            return np.nan
        return brentq(residual, beta_min, beta_max, xtol=0.01)
    except Exception:
        return np.nan

# ── step 4: efficiency curve for CSP propeller ────────────────────────────
def csp_efficiency(D_ft, rho_si, P_bhp, RPM,
                   Ct_interp, Cp_interp,
                   beta_min=10, beta_max=40):
    """
    At fixed RPM (constant speed prop), sweep airspeed.
    Governor adjusts beta to always absorb P_bhp.
    Returns airspeed, efficiency, selected beta, Ct, Cp.
    """
    rho_slug  = rho_si * kgtoslug
    n         = RPM / 60.0                          # rev/s
    Cp_req    = P_bhp * 550 / (rho_slug * n**3 * D_ft**5)

    V_knots   = np.linspace(20, 150, 300)
    V_fps     = V_knots * ktastofps

    eta_out   = np.full_like(V_fps, np.nan)
    beta_out  = np.full_like(V_fps, np.nan)
    Ct_out    = np.full_like(V_fps, np.nan)

    for i, V in enumerate(V_fps):
        J   = V / (n * D_ft)
        lam = J / np.pi

        beta = solve_beta(J, Cp_req, Cp_interp, beta_min, beta_max)
        if np.isnan(beta):
            continue

        Ct  = float(Ct_interp(beta, J))
        if Ct < 0:
            continue

        # actuator disk efficiency with swirl correction
        A   = np.pi * (D_ft / 2)**2
        q   = 0.5 * rho_slug * V**2
        T   = Ct * rho_slug * n**2 * D_ft**4

        num = 2 * (1 - lam**2 * np.log(1 + 1/lam**2))
        den = 1 + np.sqrt(1 + T/(q*A)) - 2*lam**2 * np.log(1 + 1/lam**2)
        eta_out[i]  = num / den

        beta_out[i] = beta
        Ct_out[i]   = Ct

    return V_knots, eta_out, beta_out, Ct_out, Cp_req

# ── step 5: plot efficiency families at fixed beta (sanity check) ──────────
def plot_fixed_beta_families(J_common, Ct_interp, Cp_interp, beta_list):
    """
    Plot eta = Ct/Cp * J for each fixed beta.
    Use this to sanity-check your JavaProp data before running CSP analysis.
    """
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    for beta in beta_list:
        Ct = Ct_interp(beta, J_common).flatten()
        Cp = Cp_interp(beta, J_common).flatten()
        eta = np.where(Cp > 0.001, Ct / Cp * J_common, np.nan)

        axes[0].plot(J_common, Ct,  label=f"β={beta}°")
        axes[1].plot(J_common, Cp,  label=f"β={beta}°")
        axes[2].plot(J_common, eta, label=f"β={beta}°")

    for ax, ylabel, title in zip(axes,
            ["Ct", "Cp", "η"],
            ["Thrust coefficient", "Power coefficient", "Efficiency"]):
        ax.set_xlabel("J  [–]")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(fontsize=8)
        ax.grid(True, lw=0.4)

    plt.tight_layout()
    plt.show()

# ── step 6: plot CSP results ───────────────────────────────────────────────
def plot_csp_results(V_knots, eta, beta_sel, title="CSP propeller — cruise"):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 6), sharex=True)

    ax1.plot(V_knots, eta, color="#185FA5", lw=2)
    ax1.axvline(132, color="gray", lw=1, ls="--", label="V cruise = 132 kt")
    ax1.set_ylabel("Efficiency η  [–]")
    ax1.set_ylim(0, 1.0)
    ax1.legend()
    ax1.grid(True, lw=0.4)
    ax1.set_title(title)

    ax2.plot(V_knots, beta_sel, color="#0F6E56", lw=2)
    ax2.axvline(132, color="gray", lw=1, ls="--")
    ax2.set_xlabel("Airspeed  [knots]")
    ax2.set_ylabel("Blade angle β  [°]")
    ax2.set_title("Governor blade angle schedule")
    ax2.grid(True, lw=0.4)

    plt.tight_layout()
    plt.show()

# ── main ───────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    # 1. load data
    J_common, all_Ct, all_Cp = load_javaprop(beta_angles)

    # 2. build interpolators
    Ct_interp, Cp_interp = build_interpolators(
        beta_angles, J_common, all_Ct, all_Cp
    )

    # 3. sanity check — fixed beta families should look like textbook prop charts
    plot_fixed_beta_families(J_common, Ct_interp, Cp_interp, beta_angles)

    # 4. your propeller RPM — fix this to your actual gearbox output RPM
    RPM = 2300   #Max propeller RPM

    # 5. run CSP analysis at cruise
    V, eta, beta_sel, Ct, Cp_req = csp_efficiency(
        D_ft, rho_cruise, P_bhp_cr, RPM,
        Ct_interp, Cp_interp
    )
    plot_csp_results(V, eta, beta_sel, title="CSP — cruise (8000 ft, 137 hp)")

    # 6. same at takeoff
    V, eta, beta_sel, Ct, Cp_req = csp_efficiency(
        D_ft, rho_TO, P_bhp_to, RPM,
        Ct_interp, Cp_interp
    )
    plot_csp_results(V, eta, beta_sel, title="CSP — takeoff (2000 ft, 160 hp)")