# Getting Started

Pysic-rs is a high-performance mathematical physics engine written in Rust with Python
bindings via PyO3. This page walks you through prerequisites, installation, a first import,
and running the test suite.

---

## Prerequisites

| Tool | Version | Purpose |
|------|---------|---------|
| **Rust toolchain** | stable ≥ 1.70 | compiles the core (via rustup) |
| **Python** | 3.8+ | bindings / usage |
| **maturin** | ≥ 1.0 | builds the Python extension |
| **num-complex, ndarray, rustfft** | crates.io | dependency resolution is automatic |

No system C compiler or BLAS/LAPACK is required — the library is pure Rust.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/ThotDjehuty/pysic-rs.git
cd pysic-rs
```

### 2. Install maturin

```bash
pip install maturin
```

### 3. Build & install the Python package (editable, into your active venv)

```bash
maturin develop --release        # release build — fast
```

For a debug, faster-to-compile cycle during development:

```bash
maturin develop                  # debug build
```

### 4. Verify the import

```python
import pysicrs
print(pysicrs.__version__)      # 0.1.0
print(pysicrs.constants()["c"]) # 2.997925e+08
```

> **Note for CI / headless builds:** every binary artifact is compiled locally; there are no
> prebuilt wheels published yet (see [Changelog](changelog.md)).

---

## Quick test

```python
from pysicrs import gamma, schwarzschild_metric, rk4_solve

print(gamma(5.0))                      # 24.0
print(schwarzschild_metric(10.0, 1.0)) # [[-0.8, ...]]
```

Run the full test suite to validate your build:

```bash
cargo test                  # 57 tests: 45 unit + 12 integration
```

---

## Jupyter note

The companion notebooks need only `pysic-rs`, NumPy, SciPy and Matplotlib.
Confirm which build you imported before trusting a number from it:

```python
# inside the notebook — always confirm the imported build
import pysicrs
assert pysicrs.__version__ >= "0.2.0", pysicrs.__version__
```

Releases before 0.2.0 contain eight incorrect numerical routines; see the
changelog.

---

## Next steps

- [Quickstart](quickstart.md) — runnable examples for each module
- [Algorithms](algorithms/special_functions.md) — mathematical foundations
- [Installation](installation.md) — detailed build recipes, conda & Rust toolchain setup