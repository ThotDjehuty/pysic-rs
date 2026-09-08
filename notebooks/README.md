# pysic-rs — Companion Notebooks

Generic suite of 14 notebooks (kernel `rhftlab`), executed end-to-end
(0 errors, ≥1 figure per code cell), covering the entire `pysicrs` API (~48 functions).

Each code cell follows the pedagogical sandwich structure:
`Theorem / Pivot Equation / Demonstration` (PRE) followed by
`Expected Result / Graph Reading / Conclusion` (POST). Discrepancies
between binding and literature are marked **CONSTAT** (observed discrepancy, no failing assert).

| # | Notebook | Content | Status |
|---|----------|---------|--------|
| 00 | `00_quickstart_and_constants.ipynb` | API overview, CODATA constants | ✓ 7 cells · 0 err · 1 fig |
| 01 | `01_ode_solvers.ipynb` | RK4, backward Euler, leapfrog, RK45 | ✓ 19 · 0 · 5 |
| 02 | `02_pde_heat.ipynb` | Heat 1D CN / 2D ADI | ✓ 10 · 0 · 2 |
| 03 | `03_pde_wave.ipynb` | Wave FDTD 1D/2D | ✓ 10 · 0 · 2 |
| 04 | `04_pde_poisson.ipynb` | Poisson FFT 2D | ✓ 7 · 0 · 1 |
| 05 | `05_pde_quantum.ipynb` | Schrödinger split-step, eigen, Dirac | ✓ 16 · 0 · 4 |
| 06 | `06_fourier_analysis.ipynb` | FFT, Parseval, spectral derivative, Welch | ✓ 13 · 0 · 3 |
| 07 | `07_quantum_states.ipynb` | Pauli/su(2), Fock/coherent states, density | ✓ 13 · 0 · 3 |
| 08 | `08_propagators.ipynb` | Feynman propagator, free propagator | ✓ 10 · 0 · 2 |
| 09 | `09_gauge_theory.ipynb` | su(3) constants, instanton action | ✓ 10 · 0 · 2 |
| 10 | `10_topology.ipynb` | Berry, Chern, winding, skyrmion | ✓ 16 · 0 · 4 |
| 11 | `11_general_relativity.ipynb` | Schwarzschild, Christoffel, Ricci/Einstein, ADM | ✓ 16 · 0 · 4 |
| 12 | `12_classical_mechanics.ipynb` | Euler rotations, Euler equations, inertia | ✓ 13 · 0 · 3 |
| 13 | `13_electromagnetism_casimir.ipynb` | Green, Coulomb, Larmor, dipole, Compton, Casimir, Polder | ✓ 19 · 0 · 5 |

Total: **179 cells, 0 errors, 41 figures**.

## Documented CONSTATs (binding ≠ literature)

- **01 ODE**: `backward_euler` — fixed-point iteration, converges only if
  `h·|λ| < 1` (not unconditionally A-stable); `leapfrog_integrate` —
  energy drifts (2nd half-kick defective, cf. `src/ode/symplectic.rs`);
  `rk45_solve` — stub returning `None`.
- **04 Poisson**: `sol/u` ratio non-constant (233–641) in interior.
- **05 Quantum**: `schrodinger_eigen_1d` energies shifted (e.g., ∞ well:
  1301 vs 0.049); `dirac_split_step_1d` norm ×4 per step (probability
  not unitary).
- **06 Fourier**: spectral derivative order-1 ≈ 0 (only even orders
  work); Welch frequency in cycles/sample.
- **07 States**: `pauli_matrices` returns real parts only → σ_y zero matrix;
  `su2_structure_constants` f_123 = 1 (generators σᵢ/2).
- **09 Gauge**: su(3) table not fully antisymmetric (`f_147` = −1/2
  instead of ±1/2) — Rust fill bug.
- **10 Topology**: `skyrmion_number` = stub (0.0).
- **11 GR**: `ricci_scalar`/`einstein_tensor` Schwarzschild hardcoded M=1
  → R ≈ −2/r² in "vacuum" (instead of 0); ADM constraints in
  16πGρ / 8πGj convention.
- **13 EM**: `dipole_radiation` divides by `12πc³` instead of `12πc`
  (factor 1/c²).

## Execution

```bash
conda activate rhftlab
jupyter nbconvert --to notebook --inplace --execute notebooks/*.ipynb
```

Documentation rendering: see `docs/source/examples.rst` (nbsphinx, pre-executed notebooks).