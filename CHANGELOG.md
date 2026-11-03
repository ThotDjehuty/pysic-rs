# Changelog

All notable changes to `pysic-rs` are recorded here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning follows
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] — 2026-11-03

First release published to crates.io and PyPI.

This release repairs **eight incorrect numerical routines**, all of which
compiled cleanly, ran without panicking and returned finite, plausible-looking
numbers. Every repair is held in place by a test that fails on the previous
code. Accuracy figures below are maximum deviations measured against SciPy on
dense grids by the companion reproduction notebook.

### Fixed

- **`bessel_j0`, `bessel_j1`, `bessel_y0`, `bessel_y1`** — a four-term Taylor
  polynomial was used out to `|x| = 8`, and the large-argument branch applied
  `.sin()` to the coefficient polynomial instead of to the phase. `J0(7.98)`
  returned `+50.6` for a function bounded by 1. Replaced with the ascending
  series plus the Hankel asymptotic expansion, with the crossover chosen by
  measuring both branches rather than by assumption. Now `9.9e-12` / `3.4e-11`
  absolute.
- **`airy_ai`** — the large-`|x|` asymptotic form was applied for *all*
  `x >= -1`, including `x = 0`, where it divides by `zeta = 0` and returned
  `NaN`. Replaced with the Maclaurin series from `y'' = x y` on the central
  interval and the decaying/oscillatory asymptotics outside it. A `sin`/`cos`
  transposition against DLMF 9.7.9 in the oscillatory branch was also fixed.
  Now `7.7e-10` absolute on `(-20, 10]`.
- **`expint_en`** — the upward recurrence evaluated `exp(-(k-1))` where
  `exp(-x)` was intended, so every order was wrong, `n = 1` included (one order
  worse than previously documented). Replaced with the modified Lentz continued
  fraction for `x > 1` and the ascending series otherwise. Now `1.1e-14`
  relative for `n <= 6`.
- **`zeta`** — the Euler–Maclaurin correction built its rising factorial from
  `floor(s)` integers instead of the real product `s(s+1)...(s+2k-2)`, and the
  reflection formula truncated the exponent via `powi(s as i32)`. `zeta(2)` was
  correct to four digits only and `zeta(-1)` was wrong. Rewritten with the
  correct Bernoulli coefficients and functional equation; `zeta(-1) = -1/12`
  and the trivial zeros now hold. Now `1.1e-15` relative.
- **`ricci_scalar`** — `schwarzschild_metric` froze `g_phiphi = r^2` at
  `theta = pi/2`, so every angular derivative vanished and the cancellation
  that makes a vacuum solution Ricci-flat failed. The routine returned the
  analytic artefact `-2/r^2`. Added `schwarzschild_metric_at(r, theta, mass)`
  with `g_phiphi = r^2 sin^2(theta)` and used it in the binding's metric
  closure. Vacuum `|R|` is now at the finite-difference residual (`~5e-9`), on
  and off the equator.
- **`dirac_split_step_1d`** — the mass rotation was coded as
  `exp(i sin(theta/2))`, whose modulus is 1, instead of `i sin(theta/2)`. The
  step matrix degenerated to `[[-1,-1],[-1,1]]` with singular values 2 and 0,
  so the norm grew by a fixed factor per step *independently of `dt`*. Now
  unitary: the norm ratio is `1` to `4e-14` across four decades of `dt`.
- **`schrodinger_eigen_1d`** — the "Jacobi sweep" computed a malformed rotation
  angle and updated off-diagonals from themselves
  (`h_off[p-1] = c*h_off[p-1] - s*h_off[p-1]`), so it never deflated and
  returned values on the order of the matrix diagonal, `O(1/dx^2)` rather than
  the physical spectrum (`32771` where `4.93` was expected). A plane-rotation
  sweep is in any case the wrong shape for a tridiagonal. Replaced with QL
  with implicit shifts. `E_1` now matches the closed form
  `2/dx^2 sin^2(k pi / 2(N+1))` to `6.5e-15`.
- **`rk45_solve`** — the Python binding ran the solve and returned
  `Ok(py.None())` behind a `TODO`, so every caller received `None`. The
  integrator itself used Fehlberg's `41/840` rather than Dormand–Prince's
  `1/40`, evaluated `f` inside the per-component fold and took component `[0]`
  of that evaluation for every component, then compared an unnormalised L2 norm
  against a per-component mean of the tolerance — so the accept test was not
  dimensionally meaningful and the integrator stalled. Rewritten with the
  correct FSAL seventh stage, embedded fourth-order weights and a scaled RMS
  error. `dy/dt = -y` now reaches `y(1) = e^-1` in 5–43 steps with the error
  tracking `rtol`.
- **`gamma`** — formed `t.powf(x + 0.5)` as a separate factor, which overflows
  near `x = 143` although `Gamma(x)` is finite out to `x ~ 171.6`, so it
  reported `inf` for representable values. Regrouped into a single exponential.
- **`ln_gamma`** — computed `gamma(x).ln()`, inheriting that overflow and
  losing mantissa to a large exponential. Now evaluates Lanczos in logarithmic
  form directly. Improved from `3.0e-12` to `1.8e-15` relative.

### Added

- **The special-function suite is now reachable from Python.** All eighteen
  functions existed in Rust but `special_functions` was never registered with
  the `PyO3` module initialiser, so none of them were callable from Python —
  including the `bessel_j0` and `beta` calls the documented binding surface
  advertises. Added `special_functions_bindings` with the documented semantics:
  a float returns a float, any sequence returns a list, two-argument functions
  broadcast a scalar against a sequence, and a length mismatch between two
  sequences **raises** rather than truncating to the shorter.
- **Chern–Simons module** (`topology::chern_simons`), graded only by what each
  routine can honestly be graded against:
  - `linking_number` — the Gauss double integral, an integer-valued knot
    invariant. Hopf link `±1`, unlink `0`, an `n`-fold wrapping `±n`.
  - `sl2_exp`, `sl2_exp_trace` — the `sl(2,R)` matrix exponential in
    Cayley–Hamilton closed form; `det exp X = 1` exactly, since `tr X = 0`.
  - `btz_holonomy_traces`, `btz_horizons_from_traces`, `btz_mass_spin` — the
    BTZ black hole as a flat `SL(2,R) x SL(2,R)` Chern–Simons connection, with
    `tr H± = 2 cosh(pi (r+ ∓ r-) / l)` and a round trip back to `(M, J)`.
  - `chern_simons_triple` — the discrete cubic functional on a simplicial
    2-skeleton. Documented explicitly as a heuristic summary statistic, **not**
    a quantised invariant, and graded only by the exact algebraic identities
    `CS(A^T) = -CS(A)` and `CS(A) = 0` for symmetric `A`.
- Tests that can fail, replacing several that could not:
  - the Bessel Wronskian identity `J1 Y0 - J0 Y1 = 2/(pi x)`, which relates all
    four routines independently of any table;
  - the bound `|J_nu| <= 1` sampled across `|x| <= 40`;
  - agreement between the series and asymptotic branches evaluated at the same
    argument (comparing `f(j-eps)` to `f(j+eps)` only measures the function's
    own slope);
  - `Ai` finite across `[-20, 10]` and matching tabulated values;
  - the exponential-integral recurrence `n E_{n+1}(x) = e^-x - x E_n(x)`;
  - `rk45_solve` required to converge *and* to improve when `rtol` is
    tightened, which a stalled integrator cannot do.

### Changed

- **Removed the `optimiz-rs` dependency.** It was unused — zero references in
  any `.rs` file — and it pulled in `openblas-build`, which both contradicted
  the documented "nothing in the core depends on BLAS or LAPACK" and broke
  source installs outright (`openblas-build requires the rustls or native-tls
  feature`). The crate is now dependency-light and builds from an sdist on a
  bare machine with only a Rust toolchain.
- The `zeta` test tolerance was `1e-3` — three orders looser than the closed
  forms justify, and loose enough to pass a routine wrong in the fourth digit.
  Now `1e-14` against `pi^2/6`, `pi^4/90` and `pi^6/945`.
- Tightened the published-crate `exclude` list: notebooks, build logs, docs
  config, wheels and distribution artifacts stay in the repository.

### Fixed (packaging)

- A local-path `[patch.crates-io]` section in `Cargo.toml` was copied into the
  source distribution, where the relative path does not exist, so
  `pip install` of the sdist failed during metadata generation. Removed along
  with the dependency it patched.

### Known limitations

- `erf` / `erfc` are limited to `1.4e-7` absolute by A&S 7.1.26, which is the
  formula's construction limit rather than a port defect.
- `Gamma` is accurate to `2.7e-13` relative on `(0, 170]` and overflows past
  `x ~ 171.6`, as the type requires.
- CPU-only by design. There is no GPU backend and none is planned.

## [0.1.0] — 2026-10-18

Initial release: special functions, linear algebra, ODE and PDE solvers,
quantum mechanics, general relativity, gauge theory, topology and
electromagnetism, with an optional `PyO3` binding layer.

> **Superseded.** Eight numerical routines in this release return incorrect
> results, and the special-function suite is not reachable from Python at all.
> See the 0.2.0 entry above. Use 0.2.0 or later.

[0.2.0]: https://github.com/ThotDjehuty/pysic-rs/releases/tag/v0.2.0
[0.1.0]: https://github.com/ThotDjehuty/pysic-rs/releases/tag/v0.1.0
