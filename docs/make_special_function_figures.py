#!/usr/bin/env python3
"""Render the special-function figures for the documentation.

Every curve is evaluated by pysic-rs itself — nothing here is plotted from
scipy or numpy's own special functions, so the figures show what the library
actually returns.

Output: docs/source/_static/figures/sf_*.svg, with a transparent background and
a mid-tone palette so they stay legible in both the light and dark Furo themes.

    python docs/make_special_function_figures.py
"""

import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import pysicrs as ps

OUT = pathlib.Path(__file__).parent / "source" / "_static" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# Mid-tone palette: readable on white and on #131416 alike.
TEAL, GOLD, RED, BLUE, PURPLE, GREEN = (
    "#2a9d8f", "#c9a227", "#c0392b", "#4a7fb5", "#8e6bb5", "#4f9d5d",
)
CYCLE = [TEAL, GOLD, RED, BLUE, PURPLE, GREEN]
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


def finish(fig, name):
    fig.tight_layout()
    path = OUT / f"{name}.svg"
    fig.savefig(path, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"  {path.relative_to(OUT.parents[3])}")


def axline(ax):
    ax.axhline(0, color=GRID, lw=0.8, alpha=0.6)


# --------------------------------------------------------------------- Gamma
def fig_gamma():
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.2))

    # Poles at 0, -1, -2, ...: plot each continuous branch separately.
    ax = axes[0]
    for lo, hi in [(-4, -3), (-3, -2), (-2, -1), (-1, 0), (0, 4.2)]:
        x = np.linspace(lo + 1e-3, hi - 1e-3, 600)
        y = np.array(ps.gamma(x))
        y[np.abs(y) > 20] = np.nan
        ax.plot(x, y, color=TEAL, lw=1.8)
    for k in range(0, -4, -1):
        ax.axvline(k, color=RED, ls=":", lw=0.9, alpha=0.7)
    n = np.arange(1, 5)
    ax.plot(n + 1, [ps.gamma(float(v + 1)) for v in n], "o", color=GOLD, ms=6,
            label=r"$\Gamma(n{+}1)=n!$", zorder=5)
    axline(ax)
    ax.set_ylim(-20, 20); ax.set_xlim(-4.2, 4.2)
    ax.set_xlabel("x"); ax.set_ylabel(r"$\Gamma(x)$")
    ax.set_title(r"$\Gamma(x)$ — poles at $0,-1,-2,\dots$")
    ax.legend(loc="lower right")

    # ln|Gamma| stays finite where Gamma overflows.
    ax = axes[1]
    x = np.linspace(0.05, 40, 800)
    ax.plot(x, ps.ln_gamma(x), color=GOLD, lw=1.8, label=r"$\ln\Gamma(x)$")
    stirling = (x - 0.5) * np.log(x) - x + 0.5 * np.log(2 * np.pi)
    ax.plot(x, stirling, "--", color=BLUE, lw=1.2, label="Stirling")
    ax.set_xlabel("x"); ax.set_ylabel(r"$\ln\Gamma(x)$")
    ax.set_title("Log-gamma vs. Stirling's leading term")
    ax.legend(loc="upper left")
    finish(fig, "sf_gamma")


# ---------------------------------------------------------------------- Beta
def fig_beta():
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.2))

    ax = axes[0]
    a = np.linspace(0.3, 6, 300)
    for b, c in zip([0.5, 1.0, 2.0, 4.0], CYCLE):
        ax.plot(a, ps.beta(a, b), color=c, lw=1.7, label=f"b = {b}")
    ax.set_yscale("log")
    ax.set_xlabel("a"); ax.set_ylabel("B(a, b)")
    ax.set_title(r"$B(a,b)=\Gamma(a)\Gamma(b)/\Gamma(a{+}b)$")
    ax.legend()

    # Symmetry B(a,b) = B(b,a) as a heat map.
    ax = axes[1]
    g = np.linspace(0.4, 4.0, 90)
    A, B = np.meshgrid(g, g)
    Z = np.array([ps.beta(g, float(bi)) for bi in g])
    im = ax.contourf(A, B, np.log10(Z), levels=18, cmap="viridis")
    ax.plot([0.4, 4.0], [0.4, 4.0], "--", color="white", lw=1.0, alpha=0.8)
    ax.set_xlabel("a"); ax.set_ylabel("b")
    ax.set_title(r"$\log_{10}B(a,b)$ — symmetric about $a=b$")
    ax.grid(False)
    fig.colorbar(im, ax=ax, shrink=0.9)
    finish(fig, "sf_beta")


# -------------------------------------------------------------------- Bessel
def fig_bessel():
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.2))

    x = np.linspace(0.0, 25, 1400)
    ax = axes[0]
    ax.plot(x, ps.bessel_j0(x), color=TEAL, lw=1.7, label=r"$J_0(x)$")
    ax.plot(x, ps.bessel_j1(x), color=GOLD, lw=1.7, label=r"$J_1(x)$")
    env = np.sqrt(2 / (np.pi * np.maximum(x, 1e-9)))
    ax.plot(x, env, ":", color=GRID, lw=1.0)
    ax.plot(x, -env, ":", color=GRID, lw=1.0, label=r"$\pm\sqrt{2/\pi x}$")
    for z in [2.404825558, 5.520078110, 8.653727913, 11.79153444]:
        ax.plot(z, 0, "v", color=RED, ms=5)
    axline(ax); ax.set_ylim(-0.8, 1.05)
    ax.set_xlabel("x"); ax.set_title("First kind — zeros marked")
    ax.legend(ncol=2, loc="upper right")

    xs = np.linspace(0.05, 25, 1400)
    ax = axes[1]
    ax.plot(xs, ps.bessel_y0(xs), color=BLUE, lw=1.7, label=r"$Y_0(x)$")
    ax.plot(xs, ps.bessel_y1(xs), color=PURPLE, lw=1.7, label=r"$Y_1(x)$")
    axline(ax); ax.set_ylim(-2.2, 0.8)
    ax.set_xlabel("x"); ax.set_title(r"Second kind — singular as $x\to0^+$")
    ax.legend(loc="lower right")
    finish(fig, "sf_bessel")


# ------------------------------------------------------------------ Legendre
def fig_legendre():
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.2))

    x = np.linspace(-1, 1, 600)
    ax = axes[0]
    for l, c in zip(range(6), CYCLE):
        ax.plot(x, ps.legendre_p(l, x), color=c, lw=1.6, label=f"$P_{l}$")
    axline(ax); ax.set_ylim(-1.15, 1.15)
    ax.set_xlabel("x"); ax.set_title(r"$P_\ell(x)$ — $\ell$ interior zeros, $P_\ell(1)=1$")
    ax.legend(ncol=3, fontsize=8)

    ax = axes[1]
    for (l, m), c in zip([(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3)], CYCLE):
        ax.plot(x, ps.legendre_plm(l, m, x), color=c, lw=1.6, label=f"$P_{l}^{m}$")
    axline(ax)
    ax.set_xlabel("x"); ax.set_title(r"Associated $P_\ell^m(x)$")
    ax.legend(ncol=3, fontsize=8)
    finish(fig, "sf_legendre")


# --------------------------------------------------------- spherical harmonics
def fig_spherical():
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.2), subplot_kw={"projection": "polar"})
    theta = np.linspace(0, np.pi, 400)
    for ax, (l, m) in zip(axes, [(1, 0), (2, 0), (3, 2)]):
        re, im = ps.spherical_harmonic(l, m, theta, 0.0)
        r = np.hypot(re, im)
        full_t = np.concatenate([theta, 2 * np.pi - theta[::-1]])
        full_r = np.concatenate([r, r[::-1]])
        ax.plot(full_t, full_r, color=TEAL, lw=1.8)
        ax.fill(full_t, full_r, color=TEAL, alpha=0.18)
        ax.set_title(rf"$|Y_{{{l}}}^{{{m}}}(\theta,0)|$", pad=12)
        ax.set_yticklabels([]); ax.grid(alpha=0.2)
    finish(fig, "sf_spherical")


# ----------------------------------------------------------------- Chebyshev
def fig_chebyshev():
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.2))
    x = np.linspace(-1, 1, 800)
    ax = axes[0]
    for n, c in zip(range(1, 6), CYCLE):
        ax.plot(x, ps.chebyshev_t(n, x), color=c, lw=1.6, label=f"$T_{n}$")
    ax.axhline(1, color=GRID, ls=":", lw=0.8); ax.axhline(-1, color=GRID, ls=":", lw=0.8)
    axline(ax); ax.set_ylim(-1.25, 1.25)
    ax.set_xlabel("x"); ax.set_title(r"$T_n(x)=\cos(n\arccos x)$ — equioscillation")
    ax.legend(ncol=3, fontsize=8)

    ax = axes[1]
    for n, c in zip(range(1, 6), CYCLE):
        ax.plot(x, ps.chebyshev_u(n, x), color=c, lw=1.6, label=f"$U_{n}$")
    axline(ax)
    ax.set_xlabel("x"); ax.set_title(r"Second kind $U_n(x)$")
    ax.legend(ncol=3, fontsize=8)
    finish(fig, "sf_chebyshev")


# ----------------------------------------------------------------------- erf
def fig_erf():
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.2))
    x = np.linspace(-3.5, 3.5, 700)
    ax = axes[0]
    ax.plot(x, ps.erf(x), color=TEAL, lw=1.8, label=r"$\mathrm{erf}(x)$")
    ax.plot(x, ps.erfc(x), color=GOLD, lw=1.8, label=r"$\mathrm{erfc}(x)$")
    ax.axhline(1, color=GRID, ls=":", lw=0.8); ax.axhline(2, color=GRID, ls=":", lw=0.8)
    axline(ax)
    ax.set_xlabel("x"); ax.set_title(r"$\mathrm{erf}+\mathrm{erfc}=1$")
    ax.legend()

    ax = axes[1]
    xt = np.linspace(0, 6, 600)
    ax.semilogy(xt, np.maximum(ps.erfc(xt), 1e-18), color=GOLD, lw=1.8,
                label=r"$\mathrm{erfc}(x)$")
    bound = np.exp(-xt**2) / (xt * np.sqrt(np.pi) + 1e-12)
    ax.semilogy(xt[1:], bound[1:], "--", color=BLUE, lw=1.2,
                label=r"$e^{-x^2}/(x\sqrt{\pi})$")
    ax.set_ylim(1e-12, 3)
    ax.set_xlabel("x"); ax.set_title("Gaussian tail decay")
    ax.legend()
    finish(fig, "sf_erf")


# ---------------------------------------------------------------------- zeta
def fig_zeta():
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.2))
    ax = axes[0]
    s = np.linspace(1.08, 10, 700)
    ax.plot(s, ps.zeta(s), color=TEAL, lw=1.8)
    ax.axhline(1, color=GRID, ls=":", lw=0.9)
    ax.axvline(1, color=RED, ls=":", lw=0.9)
    ax.plot([2, 4, 6], [np.pi**2 / 6, np.pi**4 / 90, np.pi**6 / 945], "o",
            color=GOLD, ms=6, label=r"$\zeta(2),\zeta(4),\zeta(6)$")
    ax.set_ylim(0.9, 4)
    ax.set_xlabel("s"); ax.set_ylabel(r"$\zeta(s)$")
    ax.set_title(r"$\zeta(s)\to1$ as $s\to\infty$; pole at $s=1$")
    ax.legend()

    ax = axes[1]
    ns = np.array([2, 4, 6, 8, 10])
    exact = [np.pi**2 / 6, np.pi**4 / 90, np.pi**6 / 945, np.pi**8 / 9450,
             np.pi**10 / 93555]
    got = ps.zeta(ns.astype(float))
    err = [abs(g - e) / e for g, e in zip(got, exact)]
    ax.semilogy(ns, err, "o-", color=RED, lw=1.5)
    ax.set_xlabel("even n"); ax.set_ylabel("relative error")
    ax.set_title(r"Relative error vs the $\pi^{2n}$ closed forms")
    finish(fig, "sf_zeta")


# ---------------------------------------------------------------------- Airy
def fig_airy():
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.2))
    ax = axes[0]
    x = np.linspace(-14, 4, 1600)
    ax.plot(x, ps.airy_ai(x), color=TEAL, lw=1.7, label=r"$\mathrm{Ai}(x)$")
    ax.axvspan(0, 4, color=GOLD, alpha=0.12)
    ax.text(2.0, 0.42, "classically\nforbidden", ha="center", fontsize=8, color=GRID)
    ax.plot(0, ps.airy_ai(0.0), "o", color=GOLD, ms=6)
    ax.annotate(r"$\mathrm{Ai}(0)=3^{-2/3}/\Gamma(2/3)$", xy=(0, ps.airy_ai(0.0)),
                xytext=(-9.5, 0.55), fontsize=8, color=GRID,
                arrowprops=dict(arrowstyle="->", color=GRID, lw=0.8))
    axline(ax); ax.set_ylim(-0.6, 0.75)
    ax.set_xlabel("x"); ax.set_title("Airy Ai — oscillatory left, decaying right")
    ax.legend(loc="lower right")

    ax = axes[1]
    xr = np.linspace(0.5, 12, 600)
    ax.semilogy(xr, np.maximum(ps.airy_ai(xr), 1e-30), color=TEAL, lw=1.7,
                label=r"$\mathrm{Ai}(x)$")
    xi = 2 * xr**1.5 / 3
    asym = np.exp(-xi) / (2 * np.sqrt(np.pi) * xr**0.25)
    ax.semilogy(xr, asym, "--", color=BLUE, lw=1.2,
                label=r"$e^{-\xi}/(2\sqrt{\pi}x^{1/4})$")
    ax.set_xlabel("x"); ax.set_title(r"Exponential decay, $\xi=\frac{2}{3}x^{3/2}$")
    ax.legend()
    finish(fig, "sf_airy")


# ----------------------------------------------------------------- exp integral
def fig_expint():
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.2))
    ax = axes[0]
    x = np.linspace(0.02, 5, 700)
    for n, c in zip([1, 2, 3, 4, 5], CYCLE):
        ax.plot(x, ps.expint_en(n, x), color=c, lw=1.6, label=f"$E_{n}$")
    ax.set_ylim(0, 2.2)
    ax.set_xlabel("x"); ax.set_title(r"$E_n(x)=\int_1^\infty e^{-xt}t^{-n}dt$")
    ax.legend(ncol=3, fontsize=8)

    ax = axes[1]
    xl = np.linspace(0.02, 12, 700)
    ax.semilogy(xl, np.maximum(ps.expint_e1(xl), 1e-18), color=TEAL, lw=1.7,
                label=r"$E_1(x)$")
    ax.semilogy(xl, np.exp(-xl) / xl, "--", color=BLUE, lw=1.2, label=r"$e^{-x}/x$")
    ax.set_ylim(1e-8, 1e2)
    ax.set_xlabel("x"); ax.set_title(r"$E_1(x)\sim e^{-x}/x$ for large $x$")
    ax.legend()
    finish(fig, "sf_expint")


if __name__ == "__main__":
    print(f"pysic-rs {ps.__version__} -> {OUT}")
    for f in (fig_gamma, fig_beta, fig_bessel, fig_legendre, fig_spherical,
              fig_chebyshev, fig_erf, fig_zeta, fig_airy, fig_expint):
        f()
    print("done")
