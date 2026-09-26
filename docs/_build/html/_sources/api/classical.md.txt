# API Reference: Classical Mechanics

The Python surface exposes classical-mechanics routines at top level
(`from pysicrs import euler_equations`). See also
[Algorithms: Classical Mechanics](../algorithms/classical.md).

```python
import pysicrs
```

## Rigid body

### `euler_equations`

```python
euler_equations(omega: list[float], moments: list[float], torque: list[float])
    -> list[float]
```

Rates $\dot{\boldsymbol\omega}$ from Euler's equations
$I_i\dot\omega_i = (I_j - I_k)\omega_j\omega_k + N_i$.

```python
d = euler_equations([1.0, 0.0, 2.0], [2.0, 2.0, 1.0], [0.0, 0.0, 0.0])
print(d)   # torque-free rates
```

### `inertia_tensor`

```python
inertia_tensor(masses: list[float], positions: list[list[float]]) -> list[list[float]]
```

$I_{ij}=\sum_a m_a(\|r_a\|^2\delta_{ij}-r_{a,i}r_{a,j})$.

### `euler_to_rotation`

```python
euler_to_rotation(phi: float, theta: float, psi: float) -> list[list[float]]
```

3D rotation matrix from ZXZ Euler angles $(\phi,\theta,\psi)$.

## Rust-only (not yet bound)

| function | purpose |
|----------|---------|
| `hamiltonian_h`, `canonical_equations`, `poisson_bracket`, `canonical_transform` | Hamiltonian picture |
| `lagrangian_l`, `euler_lagrange`, `kinetic_energy`, `potential_energy` | Lagrangian picture |
| `torque_free` | time-integrated free top |
| `angular_momentum` | $\mathbf L = I\boldsymbol\omega$ |
| `quaternion_to_rotation` | quaternion parameterization |