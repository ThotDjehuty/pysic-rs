# Pysic-rs

**High-performance mathematical physics engine in Rust with Python bindings.**

<p align="center">
  <img src="https://raw.githubusercontent.com/ThotDjehuty/pysic-rs/main/docs/source/logo_pysicrs.png"
       alt="Pysic-rs — Hofstadter butterfly logo"
       width="190" />
</p>

<div align="center">

[![Docs](https://img.shields.io/badge/docs-readthedocs-1e3a5f.svg)](https://pysic-rs.readthedocs.io/)
[![PyPI](https://img.shields.io/pypi/v/pysic-rs.svg)](https://pypi.org/project/pysic-rs/)
[![Crate](https://img.shields.io/crates/v/pysic-rs.svg)](https://crates.io/crates/pysic-rs)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/ThotDjehuty/pysic-rs/blob/main/LICENSE)
[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)
[![Rust](https://img.shields.io/badge/rust-1.70+-orange.svg)](https://www.rust-lang.org/)

**📖 [Documentation](https://pysic-rs.readthedocs.io/) · 📦 [PyPI](https://pypi.org/project/pysic-rs/) · 🦀 [crates.io](https://crates.io/crates/pysic-rs) · 🐙 [GitHub](https://github.com/ThotDjehuty/pysic-rs)**

</div>

---

## Why Pysic-rs?

- **Fast**: CPU-only Rust core with **50–100× speedup** over pure Python, and even faster inside NumPy loops thanks to PyO3 zero-overhead bindings.
- **Textbook-validated**: every solver is checked against closed-form / textbook results — graduating from a *"should be right"* to *"is right"* guarantee.
- **Fractal, on purpose**: the logo is [Hofstadter's butterfly](https://en.wikipedia.org/wiki/Hofstadter%27s_butterfly) — the spectrum of electrons in a magnetic field. pysic-rs computes that very spectrum (and a lot more).
- **Lean**: pure Rust core, zero mandatory dependencies (NumPy for Python marshalling only), optional Rayon parallelism.

The logo is the physicist's inside joke: the Harper equation that *draws* the butterfly (below) is one of the hundreds of equations pysic-rs solves in Rust.

---

## The mathematics

Pysic-rs implements **standard, textbook mathematical physics at 50–100× speed**.
Here is a taste of what ships out of the box (rendered on GitHub in MathJax; fully typeset in the [docs](https://pysic-rs.readthedocs.io/)):

**Schrödinger equation** (`pde` module — split-step Fourier solver)

$$i\hbar\frac{\partial}{\partial t}\psi(\mathbf{r},t)=\left[-\frac{\hbar^{2}}{2m}\nabla^{2}+V(\mathbf{r})\right]\psi(\mathbf{r},t)$$

**Dirac equation** (`pde` module — finite-difference test suite against the free-particle dispersion)

$$\left(i\gamma^{\mu}\partial_{\mu}-m\right)\psi=0$$

**Harper equation → Hofstadter butterfly** (`topology` module — Chern numbers, Berry phase)

$$\psi_{n+1}+\psi_{n-1}+2\cos(2\pi n\alpha)\,\psi_{n}=E\psi_{n}$$

**Schwarzschild metric** (`general_relativity` module — Christoffel symbols, geodesics)

$$ds^{2}=-\left(1-\frac{2GM}{c^{2}r}\right)c^{2}dt^{2}+\left(1-\frac{2GM}{c^{2}r}\right)^{-1}dr^{2}+r^{2}d\Omega^{2}$$

**Casimir force** (`casimir` module — parallel plates, Lifshitz theory)

$$F=-\frac{\pi^{2}\hbar c}{240\,a^{4}}A$$

**Maxwell equations** (`em` module — Green's functions, radiation)

$$\nabla\cdot\mathbf{E}=\frac{\rho}{\varepsilon_{0}},\qquad
\nabla\times\mathbf{B}-\frac{1}{c^{2}}\frac{\partial\mathbf{E}}{\partial t}=\mu_{0}\mathbf{J}$$

> **Note on rendering**: GitHub renders a *subset* of LaTeX. For fully typeset math, environments, and custom macros, open the [documentation site](https://pysic-rs.readthedocs.io/) where everything is rendered with MathJax.

---

## Modules

| Module | Description |
|--------|-------------|
| **special_functions** | Gamma, Bessel, Legendre, erf, zeta, Airy, Chebyshev |
| **linalg** | Cholesky, LU, inverse, determinant, tensor ops, Lie brackets |
| **calculus** | gradient, Hessian, Simpson/Romberg integration, interpolation |
| **ode** | RK4, RK45 (Dormand–Prince), backward Euler, Crank–Nicolson, symplectic |
| **pde** | Schrödinger (split-step), Dirac, heat, wave, Poisson, Maxwell FDTD |
| **fourier** | FFT/IFFT, spectral derivatives, Welch PSD |
| **quantum** | Pauli algebra, density matrices, path integrals, Feynman propagators |
| **general_relativity** | Schwarzschild/Kerr/FLRW metrics, Christoffel/Riemann, geodesics, ADM |
| **gauge** | SU(2)/SU(3) structure constants, Yang–Mills, connections, curvature |
| **classical** | Hamiltonian/Lagrangian, rigid body (Euler equations), rotations |
| **topology** | Berry phase, Chern numbers, winding numbers, skyrmion charge |
| **em** | Green's functions (static, retarded, Helmholtz), radiation formulas |
| **casimir** | Parallel plates, sphere, Lifshitz theory, Polder potential |

---

## Quick Start

```python
from pysicrs import schwarzschild_metric, fft, casimir_energy_parallel_plates

# Schwarzschild metric at r = 10M
g = schwarzschild_metric(r=10.0, mass=1.0)

# FFT of a sine wave
import math
x = [math.sin(2 * math.pi * k / 256) for k in range(256)]
real, imag = fft(x, [0.0] * 256)

# Casimir energy between parallel plates
E = casimir_energy_parallel_plates(
    d=1e-6, hbar=1.054571817e-34, c=2.99792458e8
)
```

---

## Installation

**From PyPI** (Python 3.8+, prebuilt `abi3` wheels):

```bash
pip install pysic-rs
```

**From crates.io** (Rust):

```bash
cargo add pysic-rs
```

**From source**:

```bash
git clone https://github.com/ThotDjehuty/pysic-rs.git
cd pysic-rs
pip install maturin
maturin develop --release
```

---

## Development

```bash
# Run tests (Rust)
cargo test

# Run tests (Python bindings)
python -m pytest tests/

# Build docs (Sphinx / Furo)
cd docs && make html

# Build wheel
maturin build --release
```

---

## License

MIT — see [LICENSE](LICENSE).

## Author

HFThot Research Lab — <contact@hfthot-lab.eu>

---

*CPU-only, textbook-validated mathematical physics — powered by Rust, bound to Python.*
