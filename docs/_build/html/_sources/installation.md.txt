# Installation

Build and install Pysic-rs on macOS, Linux, or Windows. Two options: an **editable
development install** (for active development) or a **wheel build** (for distribution).

---

## Requirements

| Tool | Version | Notes |
|------|---------|-------|
| **Rust toolchain** | ≥ 1.70 | `rustup` recommended |
| **Python** | ≥ 3.8 | 3.11 recommended (`rhftlab`) |
| **pip** | ≥ 22 | — |
| **maturin** | ≥ 1.0 | `pip install maturin` |

---

## Install from source (recommended path)

```bash
git clone https://github.com/ThotDjehuty/pysic-rs.git
cd pysic-rs

# create/activate a virtual environment
python -m venv .venv && source .venv/bin/activate   # macOS/Linux

pip install maturin
maturin develop --release
```

Verify:

```python
import pysicrs
print(pysicrs.constants()["G"])   # 6.674e-11 m³ kg⁻¹ s⁻²
```

---

## Build a Python wheel

```bash
maturin build --release
ls target/wheels/
```

Install the resulting wheel:

```bash
pip install target/wheels/pysic_rs-0.1.0-*.whl
```

---

## Cargo features

| Feature | Description | Default |
|---------|-------------|---------|
| `python` | Build the PyO3 bindings (automatically on with maturin) | no |
| `parallel` | Enable Rayon parallelism in hot loops | no |

```bash
maturin develop --release --features parallel
```

---

## Troubleshooting

**`maturin` not found**

```bash
pip install --upgrade maturin
```

**Linker / toolchain errors on macOS**

```bash
xcode-select --install          # if Command Line Tools are missing
rustup update stable
```

**Import fails with `ModuleNotFoundError: pysicrs._core`**

Rebuild the extension — the source module is Rust, so it must be compiled:

```bash
maturin develop --release
```

**Slow release builds on first run**

The first `--release` build compiles all dependencies from source; subsequent builds are
cached by cargo.

---

### Optional interop: Optimiz-rs (statistics)

For probabilistic/statistical backgrounds, Pysic-rs integrates with
[**Optimiz-rs**](https://optimiz-rs.readthedocs.io/) (MCMC, HMM, SDE, Kalman filters). The two
libraries are designed to be used side by side from the same Python environment:

```bash
pip install maturin
maturin develop --release     # optimiz-rs
maturin develop --release     # pysic-rs (in its own checkout)
```

---

## Related

- [Getting Started](getting-started.md) — the quick tour
- [Quickstart](quickstart.md) — runnable examples
- **Read the Docs configuration** — `.readthedocs.yaml` at the repo root (CI rendering, no
  Rust build on RTD)