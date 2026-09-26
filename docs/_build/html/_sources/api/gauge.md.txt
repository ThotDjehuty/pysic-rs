# API Reference: Gauge Theory

The Python surface exposes gauge-theory routines at top level
(`from pysicrs import su2_structure_constants`). See also
[Algorithms: Gauge Theory](../algorithms/gauge.md).

```python
import pysicrs
```

## Structure constants

### `su2_structure_constants`

```python
su2_structure_constants() -> list[list[list[float]]]
```

$f^{abc} = \varepsilon^{abc}$ (totally antisymmetric).

### `su3_structure_constants`

```python
su3_structure_constants() -> list[list[list[float]]]
```

SU(3) structure constants in the Gell-Mann basis: $f^{123}=1$,
$f^{147}=f^{246}=f^{257}=f^{345}=\tfrac12$, $f^{458}=f^{678}=\tfrac{\sqrt3}{2}$.

```python
f = su3_structure_constants()
print(f[0][1][2])   # f^{123} = 1
print(f[3][4][7])   # f^{458} = √3/2 ≈ 0.8660
```

## Instanton

### `instanton_action`

```python
instanton_action(g: float) -> float
```

Instanton action $S = 8\pi^2/g^2$ (for winding number $\pm1$ SU(2) instantons).

```python
print(instanton_action(0.5))   # ≈ 315.83
```

## Rust-only (not yet bound)

| function | purpose |
|----------|---------|
| `structure_constants("su2"/"su3")` | string-selected table |
| `gauge_connection`, `covariant_derivative` | gauge machinery |
| `field_strength` | $F_{\mu\nu}^a$ |
| `yang_mills_action`, `yang_mills_eom` | action / EOM residual |
| `gell_mann_matrices` | λ₁…λ₈ |