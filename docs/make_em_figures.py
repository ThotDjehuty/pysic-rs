#!/usr/bin/env python3
"""Render the electromagnetism figures for the documentation.

Every curve is evaluated by pysic-rs itself — nothing here is plotted from
scipy or from a hand-written copy of the formula, so the figures show what the
library actually returns. Where a closed form is drawn for comparison it is
labelled as such and plotted as a dashed reference line.

Output: docs/source/_static/figures/em_*.svg, with a transparent background and
a mid-tone palette so they stay legible in both the light and dark Furo themes.

    python docs/make_em_figures.py
"""

import math
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import pysicrs as ps

OUT = pathlib.Path(__file__).parent / "source" / "_static" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# Same mid-tone palette as make_special_function_figures.py.
TEAL, GOLD, RED, BLUE, PURPLE, GREEN = (
    "#2a9d8f", "#c9a227", "#c0392b", "#4a7fb5", "#8e6bb5", "#4f9d5d",
)
GRID = "#9aa0a6"

plt.rcParams.update({
    "figure.dpi": 110,
    "savefig.transparent": True,
    "axes.edgecolor": GRID,
    "axes.labelcolor": GRID,
    "xtick.color": GRID,
    "ytick.color": GRID,
    "text.color": GRID,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.6,
    "legend.frameon": False,
    "font.size": 9,
})

# SI constants, passed explicitly into pysic-rs (the API takes them as arguments
# rather than baking them in, so the figures document the same units).
C = 2.997_924_58e8
EPS0 = 8.854_187_817e-12
MU0 = 4e-7 * math.pi
H = 6.626_070_15e-34
M_E = 9.109_383_7015e-31
Q_E = 1.602_176_634e-19


def finish(fig, name):
    fig.tight_layout()
    path = OUT / f"{name}.svg"
    fig.savefig(path, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"  {path.relative_to(OUT.parents[3])}")


# ------------------------------------------------------- Green's function
def fig_green():
    """G(r) = 1/(4πr) in the convention ∇²G = -δ³(r)."""
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.2))

    r = np.linspace(0.05, 4.0, 600)
    g = np.array([ps.green_fn_static(float(v)) for v in r])

    ax = axes[0]
    ax.plot(r, g, color=TEAL, lw=1.8, label=r"$G(r)=1/(4\pi r)$")
    ax.set_xlabel("r"); ax.set_ylabel("G(r)")
    ax.set_title(r"Static Green's function, $\nabla^2 G=-\delta^3$")
    ax.set_ylim(0, 1.7)
    ax.legend(loc="upper right")

    # Log-log: the 1/r slope is exactly -1.
    ax = axes[1]
    r2 = np.logspace(-1.3, 1.3, 400)
    g2 = np.array([ps.green_fn_static(float(v)) for v in r2])
    ax.loglog(r2, g2, color=TEAL, lw=1.8, label="pysic-rs")
    ax.loglog(r2, 1.0 / (4 * np.pi * r2), "--", color=GOLD, lw=1.2,
              label=r"$1/(4\pi r)$ (closed form)")
    slope = np.polyfit(np.log(r2), np.log(g2), 1)[0]
    ax.set_xlabel("r"); ax.set_ylabel("G(r)")
    ax.set_title(f"slope = {slope:.4f} (exact: $-1$)")
    ax.legend(loc="lower left")

    finish(fig, "em_green_static")


# ------------------------------------------------------------ Point charge
def fig_point_charge():
    """Coulomb field E(r) = Q/(4πε₀r²)."""
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.2))

    ax = axes[0]
    r = np.linspace(0.05, 2.0, 600)
    e = np.array([ps.field_strength_point_charge(Q_E, float(v), EPS0) for v in r])
    ax.plot(r, e, color=BLUE, lw=1.8)
    ax.set_xlabel("r  [m]"); ax.set_ylabel(r"$E_r$  [V/m]")
    ax.set_title(r"Coulomb field of one electron charge")
    ax.set_ylim(0, e[np.argmin(np.abs(r - 0.15))])

    ax = axes[1]
    r2 = np.logspace(-1.5, 1.5, 400)
    e2 = np.array([ps.field_strength_point_charge(Q_E, float(v), EPS0) for v in r2])
    ax.loglog(r2, e2, color=BLUE, lw=1.8, label="pysic-rs")
    ax.loglog(r2, Q_E / (4 * np.pi * EPS0 * r2**2), "--", color=RED, lw=1.2,
              label=r"$Q/(4\pi\varepsilon_0 r^2)$")
    slope = np.polyfit(np.log(r2), np.log(e2), 1)[0]
    ax.set_xlabel("r  [m]"); ax.set_ylabel(r"$E_r$  [V/m]")
    ax.set_title(f"slope = {slope:.4f} (exact: $-2$)")
    ax.legend(loc="lower left")

    finish(fig, "em_point_charge")


# ----------------------------------------------------------------- Larmor
def fig_larmor():
    """P = q²a²/(6πε₀c³): quadratic in both charge and acceleration."""
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.2))

    ax = axes[0]
    a = np.logspace(15, 22, 300)
    p = np.array([ps.larmor_formula(Q_E, float(v), EPS0, C) for v in a])
    ax.loglog(a, p, color=GOLD, lw=1.8, label="pysic-rs")
    slope = np.polyfit(np.log(a), np.log(p), 1)[0]
    ax.set_xlabel(r"$a$  [m/s$^2$]"); ax.set_ylabel("P  [W]")
    ax.set_title(f"Larmor: $P\\propto a^2$ (slope {slope:.3f})")
    ax.legend(loc="upper left")

    ax = axes[1]
    q = np.logspace(-19.5, -17.0, 300)
    p = np.array([ps.larmor_formula(float(v), 1e20, EPS0, C) for v in q])
    ax.loglog(q, p, color=PURPLE, lw=1.8, label="pysic-rs")
    slope = np.polyfit(np.log(q), np.log(p), 1)[0]
    ax.axvline(Q_E, color=GRID, ls=":", lw=0.9)
    ax.annotate("e", xy=(Q_E, p[0]), xytext=(Q_E * 1.15, p[0] * 3),
                color=GRID, fontsize=9)
    ax.set_xlabel("q  [C]"); ax.set_ylabel("P  [W]")
    ax.set_title(f"Larmor: $P\\propto q^2$ (slope {slope:.3f})")
    ax.legend(loc="upper left")

    finish(fig, "em_larmor")


# -------------------------------------------------------- Dipole radiation
def fig_dipole():
    """P ∝ ω⁴ and the sin²θ angular pattern."""
    fig = plt.figure(figsize=(9.5, 3.4))
    p0 = 1e-30

    ax = fig.add_subplot(1, 2, 1)
    omega = np.logspace(12, 16, 300)
    power = np.array([ps.dipole_radiation(p0, float(w), MU0, C) for w in omega])
    ax.loglog(omega, power, color=TEAL, lw=1.8, label="pysic-rs")
    slope = np.polyfit(np.log(omega), np.log(power), 1)[0]
    ax.set_xlabel(r"$\omega$  [rad/s]"); ax.set_ylabel("P  [W]")
    ax.set_title(f"$P\\propto\\omega^4$ (slope {slope:.3f})")
    ax.legend(loc="upper left")

    # Polar: dP/dΩ ∝ sin²θ, the classic doughnut with a null on the axis.
    ax = fig.add_subplot(1, 2, 2, projection="polar")
    theta = np.linspace(0, 2 * np.pi, 721)
    dpdo = np.array([
        ps.dipole_angular_distribution(float(t), p0, 1e15, MU0, C) for t in theta
    ])
    ax.plot(theta, dpdo, color=RED, lw=1.8)
    ax.fill(theta, dpdo, color=RED, alpha=0.12)
    ax.set_theta_zero_location("N")
    ax.set_title(r"$dP/d\Omega\propto\sin^2\theta$", pad=14)
    ax.set_yticklabels([])
    ax.grid(alpha=0.25)

    # Cross-check that the pattern integrates back to the total power.
    th = np.linspace(0, np.pi, 4001)
    integrand = np.array([
        ps.dipole_angular_distribution(float(t), p0, 1e15, MU0, C) for t in th
    ]) * 2 * np.pi * np.sin(th)
    integral = np.trapezoid(integrand, th) if hasattr(np, "trapezoid") else np.trapz(integrand, th)
    total = ps.dipole_radiation(p0, 1e15, MU0, C)
    print(f"    dipole check: ∫dP/dΩ = {integral:.6e} W vs P = {total:.6e} W "
          f"(rel. err {abs(integral-total)/total:.2e})")

    finish(fig, "em_dipole")


# ---------------------------------------------------------------- Compton
def fig_compton():
    """Δλ = λ_C(1 - cos θ), zero forward and 2λ_C backward."""
    fig, ax = plt.subplots(figsize=(5.4, 3.2))

    theta = np.linspace(0, np.pi, 400)
    dl = np.array([ps.compton_wavelength_shift(float(t), H, M_E, C) for t in theta])
    lam_c = H / (M_E * C)

    ax.plot(np.degrees(theta), dl / lam_c, color=GREEN, lw=1.8, label="pysic-rs")
    ax.axhline(2.0, color=GRID, ls=":", lw=0.9)
    ax.annotate(r"$2\lambda_C$ (backscatter)", xy=(120, 2.0), xytext=(58, 2.08),
                color=GRID, fontsize=8.5)
    ax.set_xlabel(r"scattering angle $\theta$  [deg]")
    ax.set_ylabel(r"$\Delta\lambda / \lambda_C$")
    ax.set_title(r"Compton shift, $\lambda_C=h/m_ec=%.4f$ pm" % (lam_c * 1e12))
    ax.set_xlim(0, 180); ax.set_ylim(0, 2.2)
    ax.set_xticks([0, 45, 90, 135, 180])
    ax.legend(loc="upper left")

    finish(fig, "em_compton")


if __name__ == "__main__":
    print("Rendering electromagnetism figures:")
    fig_green()
    fig_point_charge()
    fig_larmor()
    fig_dipole()
    fig_compton()
    print("Done.")
