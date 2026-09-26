# Changelog

All notable changes to **Pysic-rs** are documented here, following keep-a-changelog style
(`Added` / `Changed` / `Fixed`). The library is at **stable v0.1.0**.

---

## [0.1.0] — 2026-09-06

**Added.**

### Core modules

- **Special functions** — Gamma (Lanczos, $g=7$), ln Γ, Beta, Bessel $J_0,J_1,Y_0,Y_1$,
  Legendre $P_\ell, P_\ell^m$, spherical harmonics, Chebyshev $T,U$, Airy, erf/erfc,
  exponential integrals, Riemann ζ (Euler–Maclaurin).
- **Linear algebra** — Cholesky, LU with partial pivoting, inverse (Gauss–Jordan/LU),
  determinant with sign handling, trace/norm/cross/dot/outer, tensor raise/lower, metric
  signature, 3×3 closed-form helpers, Lie bracket.
- **Numerical calculus** — gradient, Hessian, Jacobian (central differences), trapezoid,
  Simpson, Gauss–Legendre and Romberg quadrature, linear/cubic-spline interpolation.

### ODE & PDE

- **ODE** — RK4, RK45 (Dormand–Prince, adaptive), backward Euler, Crank–Nicolson,
  leapfrog, velocity Verlet, Yoshida (4th-order symplectic), `integrate_hamiltonian`.
- **PDE** — Schrödinger (split-step FFT, time-dependent + eigenstates), Dirac (split-step),
  heat (Crank–Nicolson 1D/2D + ADI), wave (FDTD 1D/2D), Poisson (FFT + SOR), 3D Maxwell
  Yee-grid FDTD.

### Physics modules

- **Fourier** — complex FFT/IFFT, real-input helpers, power spectral density, Welch PSD,
  spectral derivatives.
- **Quantum mechanics** — Pauli matrices, spin operators, density matrices, partial trace,
  purity, inner products, projectors, coherent/number states, Rayleigh–Schrödinger
  perturbation (1st/2nd order), discretized path integral, effective potential, Wick
  rotation, Feynman propagators.
- **General relativity** — Schwarzschild, Kerr, FLRW, Minkowski metrics; Christoffel,
  Riemann, Ricci, Einstein tensors; geodesic RHS + integration; ADM 3+1 (lapse, shift,
  evolve, constraints); perfect-fluid and EM stress–energy.
- **Gauge theory** — SU(2)/SU(3) structure constants, Gell-Mann matrices, gauge connections,
  covariant derivatives, field strength, Yang–Mills action/EOM, instanton action.
- **Classical mechanics** — Hamiltonian, canonical equations, Poisson bracket, canonical
  transformations, Lagrange/Euler–Lagrange, rigid-body dynamics (Euler equations,
  torque-free tops, inertia tensors, angular momentum).
- **Topology** — Berry phase/curvature, Bloch-Lopez-type discrete evaluation, Chern/TKNN
  numbers (integral + discrete), winding numbers, skyrmion number.
- **Electromagnetism** — static/retarded/advanced Green's functions, Larmor and dipole
  radiation, angular distributions, magnetic-dipole radiation.
- **Casimir effect** — parallel-plates energy/force, sphere, zeta-regularized, cylinders,
  finite-T correction, Lifshitz formula, Polder potential.

### Infrastructure

- Python bindings via PyO3 (abi3-style, no NumPy dependency at the boundary).
- `python/pysicrs/__init__.py` — dynamic export of `_core` with clean unsuffixed aliases
  (e.g. `rk4_solve` local → `rk4_solve_py` binding).

**Added (documentation, this pass).**

- Full rewrite of the 13 algorithm pages to optimiz-rs depth with proven mathematical
  foundations (definitions, theorems, tables, implementation notes, references).
- New `theory/mathematical_foundations.md` — unifying variational/stability/index results.
- New detailed `quickstart.md` (every snippet executed against the wheel), `getting-started.md`,
  `installation.md` (wheels, features, troubleshooting), `contributing.md`, `changelog.md`.
- API pages reconciled with the real Python surface (exact `inspect.signature` output,
  Python-bound vs Rust-only status per function).
- `index.rst` rewritten — module map, features grid, quick example, toctrees for
  Algorithms/API/Theory/Advanced.
- Doc build is warning-free (`sphinx-build -b html`, Furo theme).
- Sphinx documentation on Furo theme (mermaid, copybutton, announcement bar, footer icons,
  intersphinx to Python/NumPy/SciPy).
- `.readthedocs.yaml` — ReadTheDocs rendering without a Rust build step on CI.
- `Makefile`, README, LICENSE (MIT).
- Test suite: **57 passing tests** (45 unit + 12 integration) + verified Python smoke tests.

**Changed.**

- `src/lib.rs` — cleaned `_core` module registration (O(N²) reflection removed; each module
  registers its own functions).
- `src/quantum/python_bindings.rs`, `src/gauge/python_bindings.rs` — PyO3 0.22 API fixes.

**Fixed.**

- ODE integration import (`integrate_hamiltonian` now re-exported).
- Special-function/linalg/calculus Python binding references removed in favor of
  module-scoped registration.
- Docs reconciled with the real Python surface (all snippets verified against the wheel).

---

## Future (unreleased)

**Planned.**

- Prebuilt wheels for `cp38-abi3` / `cp311` on macOS/Linux/Windows.
- General-order Bessel + complex special functions.
- Sparse linear-algebra backend.
- PML absorbing boundaries for the FDTD Maxwell solver.