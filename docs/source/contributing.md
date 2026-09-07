# Contributing

Thanks for your interest in Pysic-rs! The project values **textbook-validated physics**, fast
Rust code, and clean, reproducible documentation.

---

## Code of conduct

Be respectful and constructive. Help others level up — the library is meant to be a place to
learn physics *and* numerical methods.

---

## Development setup

```bash
git clone https://github.com/ThotDjehuty/pysic-rs.git
cd pysic-rs
pip install maturin
maturin develop          # debug build — fast iteration
```

---

## Adding a new algorithm

1. Create (or extend) the module under `src/<module>/`.
2. Implement the function following the existing style:
   - Pure Rust core operating on `&[f64]` / `&[num_complex::Complex64]`.
   - Public re-export from the module's `mod.rs`.
3. If it belongs to a documented module, add the PyO3 binding in
   `src/<module>/python_bindings.rs` and register it in `src/lib.rs`.
4. Add unit tests with closed-form reference values (see `src/special_functions` for the
   pattern). **No test → no merge.**
5. Add a documentation page under `docs/source/algorithms/` **or** extend the API page under
   `docs/source/api/`.

### Proven-physics rule

Every formula must come from a **proven, published** source (textbook / standard reference).
Experimental or research-pending models are **not** accepted. When in doubt, add the theorem
name and reference to the docstring.

---

## Running tests

```bash
cargo test                  # unit + integration
cargo test --test integration_tests
```

64 tests must pass (45 unit + 12 integration + 7 feature-independent… — the exact count is
reported by cargo; the important thing is **no failures**).

---

## Building documentation

```bash
cd docs
pip install -r requirements.txt
make html                  # or: sphinx-build -b html source build/html
open build/html/index.html
```

The docs use Sphinx + the Furo theme with `myst-parser` for the Markdown algorithm pages.
Numeric claims in docs should be reproducible from the code (they are checked by unit tests).

---

## Code style

- Rust: `rustfmt` (default settings) + conventional comments.
- Python: follow the existing `__init__.py` — clean aliases, no imports of private `_core`.
- No new external dependencies unless they're in `Cargo.toml` (keeps the wheel small and
  the build self-contained).
- Every commit compiles; run `cargo test` before pushing.

---

## Architecture

```
src/
├── lib.rs                 # crate root + PyO3 module registration
├── core.rs                # error types (PysicError)
├── constants.rs           # physical constants
├── special_functions.rs   # Gamma, Bessel, Legendre, erf, zeta, Airy, Chebyshev
├── linalg.rs              # Cholesky, LU, inverse, determinant, tensor ops
├── calculus.rs            # gradient, Hessian, quadrature, interpolation
├── ode/                   # RK4, RK45, implicit, symplectic integrators
├── pde/                   # Schrödinger, Dirac, heat, wave, Poisson, Maxwell FDTD
├── fourier/               # FFT/IFFT, spectral derivatives, PSD, Welch
├── quantum/               # Pauli, density matrices, path integrals, propagators
├── general_relativity/    # metrics, Christoffel, Riemann, geodesics, ADM
├── gauge/                 # SU(2)/SU(3), connections, Yang–Mills, instantons
├── classical/             # Hamiltonian, Lagrangian, rigid body
├── topology/              # Berry phase, Chern, winding numbers
├── em/                    # Green's functions, radiation formulas
└── casimir/               # parallel plates, sphere, cylinders, Lifshitz
```

---

## Documentation pages layout

```
docs/source/
├── index.rst                 # home page / toctree
├── getting-started.md        # prerequisites, install, first import
├── installation.md           # full build recipes, wheels, features
├── quickstart.md             # runnable examples per module
├── algorithms/*.md           # one page per module — math foundations
├── api/*.md                  # one page per module — function reference
├── theory/                   # shared mathematical foundations
├── contributing.md
└── changelog.md
```

---

## Review checklist for PRs

- [ ] `cargo test` passes (no failures)
- [ ] `rustfmt` clean
- [ ] New formula has a proven source referenced in the docstring
- [ ] Docs updated (algorithm or API page)
- [ ] No secrets/`historia/`/agent files committed

---

## Related

- [Changelog](changelog.md) — what shipped when