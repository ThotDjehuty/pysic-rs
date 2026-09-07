# API Reference: Quantum Mechanics

The Python surface exposes quantum-mechanics routines at top level
(`from pysicrs import pauli_matrices`). See also
[Algorithms: Quantum Mechanics](../algorithms/quantum.md).

```python
import pysicrs
```

## Pauli algebra & density matrices

### `pauli_matrices`

```python
pauli_matrices() -> (sx: list[list[complex]], sy, sz)
```

The three Pauli matrices as complex 2×2.

```python
sx, sy, sz = pauli_matrices()
print(sx)   # [[0,1],[1,0]]
```

### `density_matrix`

```python
density_matrix(psi_real: list[float], psi_imag: list[float]) -> list[list[complex]]
```

Projective density matrix $\rho = |\psi\rangle\langle\psi|$.

## Hilbert space

### `coherent_state`

```python
coherent_state(alpha_real: float, alpha_imag: float, n_terms: int)
    -> (real: list[float], imag: list[float])
```

Glauber coherent state $\sum_n c_n |n\rangle$ in the Fock basis
$c_n = e^{-|\alpha|^2/2}\alpha^n/\sqrt{n!}$. Returns two equal-length lists.

```python
re, im = coherent_state(1.0, 0.0, 20)
norm = sum(a*a + b*b for a, b in zip(re, im))
print(norm)           # ≈ 1.0
```

### `number_state`

```python
number_state(n: int, dim: int) -> (real: list[float], imag: list[float])
```

Fock state $|n\rangle$ in dimension `dim`.

## Propagators

### `feynman_propagator_scalar`

```python
feynman_propagator_scalar(p_squared: float, mass: float, epsilon: float = 1e-9)
    -> complex
```

Scalar momentum-space Feynman propagator $\frac{i}{p^2 - m^2 + i\varepsilon}$. Returns
`complex` (Python `complex`, not a pair).

### `free_propagator_position`

```python
free_propagator_position(r: float, mass: float) -> float
```

Euclidean 4D position-space propagator $\frac{m}{4\pi^2 r}K_1(mr)$ (massless limit
fallback $1/(4\pi^2 r^2)$ when $mr\to 0$).

## Rust-only (not yet bound)

| function | purpose |
|----------|---------|
| `spin_operators` | $S_a = \tfrac{\hbar}{2}\sigma_a$ |
| `partial_trace`, `purity` | entanglement diagnostics |
| `inner_product`, `projector` | Hilbert-space ops |
| `rayleigh_schrodinger_1st/2nd` | perturbation theory |
| `discretized_path_integral`, `effective_potential`, `wick_rotation` | path integrals |
| `free_propagator_momentum`, `photon_propagator_feynman` | field propagators |