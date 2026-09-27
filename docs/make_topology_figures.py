#!/usr/bin/env python3
"""Render the topology figures for the documentation.

Every invariant shown here is computed by pysic-rs itself. Closed forms, where
drawn, are dashed reference lines and labelled as such.

Output: docs/source/_static/figures/topo_*.svg, transparent background, mid-tone
palette so the figures stay legible in both the light and dark Furo themes.

    python docs/make_topology_figures.py
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


# ------------------------------------------------------------- Berry phase
def fig_berry():
    """γ = -½Ω = -π(1-cos θ), from both the closed form and the discrete
    Wilson-loop estimator applied to explicit spin-1/2 states."""
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.3))

    theta = np.linspace(0, np.pi, 240)

    # Closed form, straight from the library.
    gamma_closed = np.array([ps.berry_phase_bloch(float(t)) for t in theta])

    # Discrete Berry phase of |n(θ,φ)⟩ = (cos θ/2, sin θ/2 e^{iφ}) around a
    # φ: 0 → 2π loop, evaluated by pysic-rs from the state vectors themselves.
    def discrete_berry(t, n_phi=400):
        phis = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
        re = [[math.cos(t / 2), math.sin(t / 2) * math.cos(p)] for p in phis]
        im = [[0.0, math.sin(t / 2) * math.sin(p)] for p in phis]
        return ps.berry_phase(re, im)

    theta_s = np.linspace(0.05, np.pi - 0.05, 22)
    gamma_disc = np.array([discrete_berry(float(t)) for t in theta_s])
    # The Wilson loop returns the phase mod 2π; unwrap onto the same branch
    # as the closed form for a like-for-like comparison.
    gamma_disc = np.where(gamma_disc > 0, gamma_disc - 2 * np.pi, gamma_disc)

    ax = axes[0]
    ax.plot(theta, gamma_closed / np.pi, color=TEAL, lw=1.8,
            label=r"$\gamma=-\pi(1-\cos\theta)$")
    ax.plot(theta_s, gamma_disc / np.pi, "o", color=GOLD, ms=5,
            label=r"discrete $-\arg\prod\langle n_i|n_{i+1}\rangle$")
    ax.set_xlabel(r"cone half-angle $\theta$")
    ax.set_ylabel(r"$\gamma/\pi$")
    ax.set_title("Berry phase on the Bloch sphere")
    ax.set_xticks([0, np.pi / 2, np.pi])
    ax.set_xticklabels(["0", r"$\pi/2$", r"$\pi$"])
    ax.legend(loc="lower left", fontsize=8)

    err = np.abs(gamma_disc - np.array([ps.berry_phase_bloch(float(t)) for t in theta_s]))
    print(f"    berry: max |discrete - closed form| = {err.max():.3e} rad")

    # Geometric reading: the phase is half the enclosed solid angle.
    ax = axes[1]
    solid = 2 * np.pi * (1 - np.cos(theta))
    ax.plot(solid / np.pi, -gamma_closed / np.pi, color=PURPLE, lw=1.8,
            label=r"pysic-rs")
    ax.plot(solid / np.pi, solid / (2 * np.pi), "--", color=GRID, lw=1.2,
            label=r"$|\gamma| = \Omega/2$")
    ax.set_xlabel(r"enclosed solid angle $\Omega/\pi$")
    ax.set_ylabel(r"$|\gamma|/\pi$")
    ax.set_title(r"Phase is exactly half the solid angle")
    ax.legend(loc="upper left")

    finish(fig, "topo_berry_phase")


# ---------------------------------------------------------- Winding number
def fig_winding():
    """n = (1/2π)∮dθ, counted by pysic-rs for curves that wind n times."""
    ns = [1, 2, 3, -2]
    fig, axes = plt.subplots(1, 4, figsize=(11.0, 2.9))

    for ax, n, color in zip(axes, ns, CYCLE):
        t = np.linspace(0.0, 2.0 * np.pi, 3000, endpoint=False)
        # Radius modulation keeps |z| > 0 while separating the n loops visually.
        radius = 1.0 + 0.3 * np.cos(t)
        f1 = radius * np.cos(n * t)
        f2 = radius * np.sin(n * t)

        computed = ps.winding_number(list(f1), list(f2))

        ax.plot(f1, f2, color=color, lw=1.4)
        ax.plot([0], [0], "+", color=RED, ms=9, mew=1.6)
        ax.set_aspect("equal")
        ax.set_title(f"pysic-rs: n = {computed:+.0f}", fontsize=9.5)
        ax.set_xticks([]); ax.set_yticks([])
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_alpha(0.35)
        assert abs(computed - n) < 1e-9, f"winding mismatch: {computed} vs {n}"

    fig.suptitle("Winding number about the origin (+ marks the origin)",
                 y=1.04, fontsize=10)
    finish(fig, "topo_winding")


# ----------------------------------------------------- Chern / QWZ plateaus
def _qwz_unit_vector_grid(m, n_k):
    """d̂(k) for the Qi–Wu–Zhang model on a periodically padded BZ grid.

    h(k) = (sin kx, sin ky, m + cos kx + cos ky). The returned grid has one
    extra ring of points so that skyrmion_number's interior stencil (which
    skips the outermost ring) covers exactly one full Brillouin zone.
    """
    dk = 2.0 * np.pi / n_k
    idx = np.arange(-1, n_k + 1)          # -1 .. n_k  → padded by one each side
    ks = -np.pi + idx * dk
    grid = []
    for kx in ks:
        row = []
        for ky in ks:
            dx, dy = math.sin(kx), math.sin(ky)
            dz = m + math.cos(kx) + math.cos(ky)
            norm = math.sqrt(dx * dx + dy * dy + dz * dz)
            if norm < 1e-12:
                # d(k) = 0 means the gap closes: the band is degenerate at this
                # k-point and the Chern number is undefined. This happens only
                # at m = 0, ±2, the transition points themselves.
                raise ValueError(
                    f"QWZ model is gapless at m={m} (k=({kx:.4f},{ky:.4f})); "
                    "the Chern number is undefined at a transition point"
                )
            row.append([dx / norm, dy / norm, dz / norm])
        grid.append(row)
    return grid, dk


def fig_chern():
    """The Chern number of a two-band model is the skyrmion number of d̂(k)
    over the Brillouin zone — so the same pysic-rs routine computes both."""
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.3))

    # Plateaus: C = 0 for |m| > 2, and ∓1 inside.
    ax = axes[0]
    n_k = 96
    # Offset the sweep so it never lands exactly on m = 0, ±2, where the gap
    # closes and the invariant is undefined (see _qwz_unit_vector_grid).
    ms = np.linspace(-3.35, 3.35, 68)
    cs = []
    for m in ms:
        grid, dk = _qwz_unit_vector_grid(float(m), n_k)
        cs.append(ps.skyrmion_number(grid, dk, dk))
    cs = np.array(cs)

    ax.plot(ms, cs, color=TEAL, lw=1.8, label=f"pysic-rs ({n_k}×{n_k} BZ grid)")
    for level in (-1, 0, 1):
        ax.axhline(level, color=GRID, ls=":", lw=0.8, alpha=0.7)
    for edge in (-2, 0, 2):
        ax.axvline(edge, color=RED, ls="--", lw=0.9, alpha=0.6)
    ax.set_xlabel("mass parameter m")
    ax.set_ylabel("Chern number C")
    ax.set_title("Qi–Wu–Zhang: quantized plateaus")
    ax.set_ylim(-1.45, 1.45)
    ax.legend(loc="upper right", fontsize=8)

    for label, m in (("C=0", -3.0), ("C=+1", -1.0), ("C=-1", 1.0), ("C=0", 3.0)):
        ax.annotate(label, xy=(m, 0.0), xytext=(m - 0.42, 1.13),
                    color=GRID, fontsize=8.5)

    # Report how close the plateaus sit to integers.
    for m in (-3.0, -1.0, 1.0, 3.0):
        grid, dk = _qwz_unit_vector_grid(m, n_k)
        c = ps.skyrmion_number(grid, dk, dk)
        print(f"    QWZ m={m:+.1f}: C = {c:+.5f}")

    # Quantization sharpens as the k-grid is refined.
    ax = axes[1]
    for m, color, want in ((-1.0, GOLD, +1.0), (1.0, BLUE, -1.0)):
        grids = [24, 32, 48, 64, 96, 128]
        errs = []
        for nk in grids:
            grid, dk = _qwz_unit_vector_grid(m, nk)
            errs.append(abs(ps.skyrmion_number(grid, dk, dk) - want))
        ax.loglog(grids, errs, "o-", color=color, lw=1.5, ms=5,
                  label=f"m = {m:+.0f}  (C = {want:+.0f})")
    ax.set_xlabel(r"BZ grid points per axis $N_k$")
    ax.set_ylabel(r"$|C - C_{\rm exact}|$")
    ax.set_title("Quantization error vs resolution")
    ax.legend(loc="lower left", fontsize=8)

    finish(fig, "topo_chern_qwz")


# ------------------------------------------------ Chern from a curvature grid
def fig_chern_integrator():
    """chern_number() integrates a supplied Berry curvature and divides by 2π.
    Feeding a constant curvature Ω = 2πC/L² must return exactly C."""
    fig, ax = plt.subplots(figsize=(5.4, 3.2))

    n, l = 64, 1.0
    dk = l / n
    targets = np.arange(-3, 4)
    got = []
    for c_target in targets:
        curvature = [[2 * np.pi * c_target / (l * l)] * n for _ in range(n)]
        got.append(ps.chern_number(curvature, dk))

    ax.plot(targets, got, "o", color=TEAL, ms=7, label="pysic-rs")
    ax.plot(targets, targets, "--", color=GRID, lw=1.2, label="exact")
    ax.set_xlabel(r"prescribed $C$ via $\Omega = 2\pi C/L^2$")
    ax.set_ylabel(r"$C$ returned by `chern_number`")
    ax.set_title(r"$C=\frac{1}{2\pi}\int\Omega\,d^2k$ on a $64\times64$ grid")
    ax.legend(loc="upper left")
    resid = np.max(np.abs(np.array(got) - targets))
    print(f"    chern_number integrator: max |error| = {resid:.2e}")

    finish(fig, "topo_chern_integrator")


# ------------------------------------------------------------ Skyrmion
def _bp_grid(span, n):
    """Belavin–Polyakov profile n̂ = (2x, 2y, r²-1)/(r²+1), |Q| = 1."""
    coords = np.linspace(-span, span, n)
    h = float(coords[1] - coords[0])
    grid = []
    for x in coords:
        row = []
        for y in coords:
            r2 = float(x * x + y * y)
            den = r2 + 1.0
            row.append([2 * float(x) / den, 2 * float(y) / den, (r2 - 1.0) / den])
        grid.append(row)
    return grid, h, coords


def fig_skyrmion():
    """The hedgehog texture and the convergence of its topological charge."""
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.6))

    # Texture: in-plane components as arrows, n_z as background colour.
    span, n = 3.0, 33
    grid, h, coords = _bp_grid(span, n)
    arr = np.array(grid)                      # [ix][iy][3]
    nz = arr[:, :, 2].T
    ax = axes[0]
    im = ax.imshow(nz, origin="lower", cmap="coolwarm", vmin=-1, vmax=1,
                   extent=[-span, span, -span, span], alpha=0.85)
    s = 2
    ax.quiver(coords[::s], coords[::s],
              arr[::s, ::s, 0].T, arr[::s, ::s, 1].T,
              color="#20222a", scale=22, width=0.005)
    ax.set_title(r"Belavin–Polyakov texture ($n_z$ shaded)")
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax.grid(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label(r"$n_z$")
    cb.outline.set_edgecolor(GRID)

    # Convergence: the residual is the analytic truncation tail 1/(1+R²).
    ax = axes[1]
    spans = [4.0, 6.0, 8.0, 12.0, 16.0, 24.0]
    errs, tails = [], []
    for sp in spans:
        npts = int(2 * sp / 0.05) | 1          # hold h ≈ 0.05, odd point count
        g, hh, _ = _bp_grid(sp, npts)
        q = abs(ps.skyrmion_number(g, hh, hh))
        errs.append(abs(q - 1.0))
        tails.append(1.0 / (1.0 + sp * sp))
        print(f"    skyrmion R={sp:4.1f} (h={hh:.4f}): |Q| = {q:.6f}")

    ax.loglog(spans, errs, "o-", color=TEAL, lw=1.6, ms=6, label="pysic-rs error")
    ax.loglog(spans, tails, "--", color=GOLD, lw=1.4,
              label=r"analytic tail $1/(1+R^2)$")
    ax.set_xlabel("domain half-width R")
    ax.set_ylabel(r"$|\,|Q| - 1\,|$")
    # At small R the error IS the truncated tail; by R ≈ 20 the fixed h = 0.05
    # discretization error takes over and the two curves part company.
    ax.set_title(r"$|Q|\to1$: tail-limited, then $h$-limited")
    ax.legend(loc="lower left", fontsize=8)

    finish(fig, "topo_skyrmion")


if __name__ == "__main__":
    print("Rendering topology figures:")
    fig_berry()
    fig_winding()
    fig_chern()
    fig_chern_integrator()
    fig_skyrmion()
    print("Done.")
