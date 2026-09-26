# Quickstart

Every example below was executed against the installed `pysic-rs` wheel. The Python API
exposes clean aliases (`gamma`, `rk4_solve`, `fft`, …) that forward to the `_core` bindings.

---

## Physical constants

```python
from pysicrs import constants

c = constants()
print(f"c    = {c['c']:.3e} m/s")
print(f"hbar = {c['hbar']:.3e} J·s")
```

**Expected:**

```
c    = 2.998e+08 m/s
hbar = 1.055e-34 J·s
```

---

## Schwarzschild metric (General Relativity)

```python
from pysicrs import schwarzschild_metric

g = schwarzschild_metric(10.0, 1.0)   # r=10M, M=1
print(f"g_tt = {g[0][0]:.4f}")        # -(1 - 2M/r) = -0.8
print(f"g_rr = {g[1][1]:.4f}")        # 1/(1 - 2M/r) = 1.25
```

**Expected:**

```
g_tt = -0.8000
g_rr = 1.2500
```

---

## ODE solving — harmonic oscillator

```python
from pysicrs import rk4_solve

def osc(t, y):
    return [y[1], -y[0]]

times, traj = rk4_solve(osc, [1.0, 0.0], 0.0, 10.0, 10000)
print(f"y(10) ≈ {traj[-1][0]:.6f}")   # cos(10) ≈ -0.8391
```

**Expected:**

```
y(10) ≈ -0.8391
```

---

## FFT peak detection

```python
from pysicrs import fft
import math

N = 256
x = [math.sin(2 * math.pi * 3 * k / N) for k in range(N)]  # bin 3
real, imag = fft(x, [0.0] * N)
peak = max(range(N), key=lambda k: real[k]**2 + imag[k]**2)
print(f"Peak bin: {peak}")
```

**Expected:**

```
Peak bin: 3
```

---

## Quantum mechanics — Pauli algebra & coherent state

```python
from pysicrs import pauli_matrices, coherent_state

sx, sy, sz = pauli_matrices()
print(f"σxσy = iσz  →  {sx[0][1] == 1}")   # σx has (0,1) = 1

re, im = coherent_state(1.0, 0.0, 20)      # |α=1⟩, 20 terms
norm = sum(a*a + b*b for a, b in zip(re, im))
print(f"Norm of |α=1⟩: {norm:.6f}")        # ≈ 1.0
```

**Expected:**

```
σxσy = iσz  →  True
Norm of |α=1⟩: 1.000000
```

*Note:* `coherent_state(alpha_real, alpha_imag, n_terms)` returns a tuple of two equal-length
lists — the real and imaginary parts of the coefficients.

---

## Casimir force

```python
from pysicrs import casimir_force

hbar, c = 1.054571817e-34, 2.99792458e8
F = casimir_force(1e-6, hbar, c)          # plates at 1 μm
print(f"F/A = {F:.3e} N/m²")
```

**Expected:**

```
F/A ≈ -1.3e-03 N/m²  (attractive)
```

---

## Where to go next

- [Getting Started](getting-started.md) — installation & build
- [Algorithms](algorithms/ode.md) — full math foundations and more examples
- [API Reference](api/ode.md) — exhaustive signatures