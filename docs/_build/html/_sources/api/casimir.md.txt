# API Reference: Casimir Effect

The Python surface exposes Casimir routines at top level
(`from pysicrs import casimir_force`). See also
[Algorithms: Casimir Effect](../algorithms/casimir.md).

```python
import pysicrs
```

## Parallel plates

### `casimir_energy_parallel_plates`

```python
casimir_energy_parallel_plates(d: float, hbar: float, c: float) -> float
```

$\dfrac{E}{A} = -\dfrac{\pi^2\hbar c}{720 d^3}$ (per unit area).

### `casimir_force`

```python
casimir_force(d: float, hbar: float, c: float) -> float
```

$\dfrac{F}{A} = -\dfrac{\pi^2\hbar c}{240 d^4}$ (per unit area — attractive).

```python
hbar, c = 1.054571817e-34, 2.99792458e8
F = casimir_force(1e-6, hbar, c)
print(F)   # ≈ -1.3e-03 N/m²
```

## Other geometry / regularization

### `polder_potential`

```python
polder_potential(distance: float, polarizability: float, hbar: float, c: float,
                 retarded: bool) -> float
```

Retarded (distance ≫ λ) or non-retarded van der Waals / Casimir-Polder energy for an atom
near a surface.

## Rust-only (not yet bound)

| function | formula |
|----------|---------|
| `casimir_energy_sphere(a)` | $-\frac{\hbar c}{8\pi a}$ |
| `casimir_energy_zeta(d)` | $-\frac{\zeta(3)\hbar c}{8\pi d^3}$ |
| `casimir_energy_cylinders(a, d)` | $-\frac{\hbar c\,a}{24 d^2}$ |
| `casimir_energy_finite_temperature(d, T)` | $E(0)+\frac{\pi^2 k_B^4 T^4}{45\hbar^3 c^3}V$ |
| `casimir_energy_lifshitz(d, ε(iω), …)` | Lifshitz integral |
| `van_der_waals(r, c6)` | $-\frac{C_6}{r^6}$ |