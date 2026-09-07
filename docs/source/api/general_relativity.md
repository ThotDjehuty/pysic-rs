# API Reference: General Relativity

The Python surface exposes GR routines at top level (`from pysicrs import schwarzschild_metric`).
See also [Algorithms: General Relativity](../algorithms/general_relativity.md).

```python
import pysicrs
```

## Metrics

### `schwarzschild_metric`

```python
schwarzschild_metric(r: float, mass: float) -> list[list[float]]
```

Schwarzschild metric in geometric units: $g_{tt} = -(1-\frac{2M}{r})$,
$g_{rr} = (1-\frac{2M}{r})^{-1}$.

```python
g = schwarzschild_metric(10.0, 1.0)
print(g[0][0])   # -0.8
print(g[1][1])   # 1.25
```

## Connection & curvature

### `christoffel_from_metric`

```python
christoffel_from_metric(coords: list[list[list[float]]], h: float = 1e-5)
    -> list[list[list[float]]]
```

Christoffel symbols $\Gamma^{\rho}_{\mu\nu}$ from a metric sampled on a grid `coords`
(indexed `[x][mu][nu]`) using central finite differences with step `h`. Returns
`[rho][mu][nu]`.

### `ricci_scalar`

```python
ricci_scalar(coords: list[list[list[float]]], h: float = 1e-5) -> float
```

$R = g^{\mu\nu}R_{\mu\nu}$. Vanishes for vacuum (e.g. Schwarzschild off the horizon).

### `einstein_tensor`

```python
einstein_tensor(coords: list[list[list[float]]], h: float = 1e-5) -> list[list[float]]
```

$G_{\mu\nu} = R_{\mu\nu} - \tfrac12 R g_{\mu\nu}$.

## ADM constraints

### `hamiltonian_constraint`

```python
hamiltonian_constraint(gamma: list[list[float]], k: list[list[float]], rho: float) -> float
```

$\mathcal H = R^{(3)} - K_{ij}K^{ij} + K^2 - 16\pi\rho$ (vanishes on constraint surface).

### `momentum_constraint`

```python
momentum_constraint(gamma: list[list[float]], k: list[list[float]], j: list[float]) -> list[float]
```

$\mathcal M_i = D^jK_{ij} - D_iK - 8\pi J_i$.

## Complete Rust surface (binding status)

Fully described in [Algorithms: General Relativity](../algorithms/general_relativity.md):

| function | Python bound? |
|----------|---------------|
| `kerr_metric`, `flrw_metric`, `minkowski_metric` | Rust (not yet bound) |
| `riemann_tensor`, `ricci_tensor` | Rust (not yet bound) |
| `geodesic_equation_rhs`, `integrate_geodesic` | Rust (not yet bound) |
| `MetricSlice`, `AdmMetric`, `evolve_metric` | Rust (not yet bound) |
| `perfect_fluid_stress_energy`, `electromagnetic_stress_energy` | Rust (not yet bound) |