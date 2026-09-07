# API Reference: Topology

The Python surface exposes topology routines at top level
(`from pysicrs import berry_phase`). See also
[Algorithms: Topology](../algorithms/topology.md).

```python
import pysicrs
```

## Berry phase & curvature

### `berry_phase`

```python
berry_phase(states_real: list[list[float]], states_imag: list[list[float]]) -> float
```

Discrete Berry phase around a closed loop of parameter points. `states_real[k]` /
`states_imag[k]` are the wavefunction components sampled along the loop.

### `berry_phase_bloch`

```python
berry_phase_bloch(theta: float) -> float
```

Berry phase accumulated by a two-level (Bloch) system along a parameter path at polar
angle `theta`.

## Chern / winding / skyrmion numbers

### `chern_number`

```python
chern_number(berry_curvature: list[list[float]], dk: float) -> float
```

$C = \frac{1}{2\pi}\int \Omega\,d^2k$ over the Brillouin zone mesh.

### `winding_number`

```python
winding_number(f1: list[float], f2: list[float]) -> float
```

$w = \frac{1}{2\pi i}\oint g'/g\,dz$ for $g = f_1 + i f_2$ sampled on a closed contour.

### `skyrmion_number`

```python
skyrmion_number(n_field: list[list[list[float]]], dx: float, dy: float) -> float
```

Degree of the map $\mathbb R^2 \to S^2$ from the spin texture `n_field` on an
`(x, y, 3)` grid.

## Rust-only (not yet bound)

| function | purpose |
|----------|---------|
| `berry_curvature` | curvature at a point |
| `chern_number_discrete` | lattice flux version |
| `winding_number_discrete` | discrete encirclement count |