# pysic-rs

**High-performance mathematical physics engine in Rust with Python bindings.**

[![Version](https://img.shields.io/badge/version-0.1.0-blue.svg)](https://github.com/ThotDjehuty/pysic-rs/releases)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/ThotDjehuty/pysic-rs/blob/main/LICENSE)
[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)
[![Rust](https://img.shields.io/badge/rust-1.70+-orange.svg)](https://www.rust-lang.org/)

Pysic-rs provides CPU-only, textbook-validated implementations of standard mathematical physics with **50-100x speedup** over pure Python.

---

## Modules

| Module | Description |
|--------|-------------|
| **special_functions** | Gamma, Bessel, Legendre, erf, zeta, Airy, Chebyshev |
| **linalg** | Cholesky, LU, inverse, determinant, tensor ops, Lie brackets |
| **calculus** | gradient, Hessian, Simpson/Romberg integration, interpolation |
| **ode** | RK4, RK45 (Dormand-Prince), backward Euler, Crank-Nicolson, symplectic |
| **pde** | Schrödinger (split-step), Dirac, heat, wave, Poisson, Maxwell FDTD |
| **fourier** | FFT/IFFT, spectral derivatives, Welch PSD |
| **quantum** | Pauli algebra, density matrices, path integrals, Feynman propagators |
| **general_relativity** | Schwarzschild/Kerr/FLRW metrics, Christoffel/Riemann, geodesics, ADM |
| **gauge** | SU(2)/SU(3) structure constants, Yang-Mills, connections, curvature |
| **classical** | Hamiltonian/Lagrangian, rigid body (Euler equations), rotations |
| **topology** | Berry phase, Chern numbers, winding numbers, skyrmion charge |
| **em** | Green's functions (static, retarded, Helmholtz), radiation formulas |
| **casimir** | Parallel plates, sphere, Lifshitz theory, Polder potential |

---

## Quick Start

```python
from pysicrs import gamma, schwarzschild_metric, fft

# Gamma function
print(gamma(5.0))  # 24.0

# Schwarzschild metric at r=10M
g = schwarzschild_metric(r=10.0, mass=1.0)

# FFT
import math
x = [math.sin(2 * math.pi * k / 256) for k in range(256)]
real, imag = fft(x, [0.0] * 256)
```

---

## Installation

```bash
git clone https://github.com/ThotDjehuty/pysic-rs.git
cd pysic-rs
pip install maturin
maturin develop --release
```

---

## Features

- **Pure Rust core** with zero-cost abstractions
- **Python bindings** via PyO3
- **Optional parallelism** with Rayon
- **Sphinx documentation** with Furo theme
- **ReadTheDocs** integration for autodeploy

---

## Development

```bash
# Run tests
cargo test

# Build docs
cd docs && make html

# Build wheel
maturin build --release
```

---

## License

MIT — see [LICENSE](LICENSE).

## Author

HFThot Research Lab — <contact@hfthot-lab.eu>
