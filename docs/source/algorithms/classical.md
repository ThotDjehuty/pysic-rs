# Classical Mechanics

Pysic-rs implements the two classic variational formulations of mechanics — Hamiltonian and
Lagrangian — plus rigid-body dynamics with the full suite of rotation parameterizations.
Everything here is textbook (proven) mechanics.

> **Python binding status:** `euler_equations`, `inertia_tensor`, and `euler_to_rotation` are
> bound. Hamiltonian/Lagrangian picture, Poisson brackets, and `torque_free` are Rust-only.
> See the [API page](../api/classical.md).

---

## Mathematical Foundations

### Hamiltonian Mechanics

For a Hamiltonian $H(q, p, t)$ the motion is governed by **Hamilton's canonical equations**
(proven — they follow from the variational principle $\delta\int(p\dot q - H)dt=0$):

$$
\dot q_i = \frac{\partial H}{\partial p_i}, \qquad
\dot p_i = -\frac{\partial H}{\partial q_i}
$$

**Poisson bracket** of two observables:

$$
\{f, g\} = \sum_i \left(\frac{\partial f}{\partial q_i}\frac{\partial g}{\partial p_i}
- \frac{\partial f}{\partial p_i}\frac{\partial g}{\partial q_i}\right)
$$

**Proven structure:** the bracket is bilinear, antisymmetric, obeys the Jacobi identity, and
$\dfrac{df}{dt} = \{f,H\} + \dfrac{\partial f}{\partial t}$. An observable with
$\{f,H\}=0$ is a constant of motion (Noether's theorem).

**Canonical transformations** preserve the symplectic form $\Omega = dq_i\wedge dp_i$; the
necessary and sufficient condition (proven) is that the Poisson brackets are invariant,
$\{q'_i, p'_j\}_{\text{new}} = \delta_{ij}$.

### Lagrangian Mechanics

The Lagrangian $L(q, \dot q, t)$ enters the **Euler–Lagrange equation** (from
$\delta\int L\,dt=0$, Hamilton's principle):

$$
\frac{d}{dt}\frac{\partial L}{\partial \dot q_i} - \frac{\partial L}{\partial q_i} = 0
$$

**Legendre transform** $H = \sum_i p_i \dot q_i - L$ with $p_i = \partial L/\partial\dot q_i$
converts between the two pictures.

**Kinetic & potential energy** for standard potentials: $T = \frac12 m v^2$, $V = V(q)$.

### Rigid Body

**Euler's equations** for the angular velocity in the body frame (with principal moments
$I_1, I_2, I_3$):

$$
I_1 \dot\omega_1 = (I_2 - I_3)\omega_2\omega_3, \qquad
I_2 \dot\omega_2 = (I_3 - I_1)\omega_3\omega_1, \qquad
I_3 \dot\omega_3 = (I_1 - I_2)\omega_1\omega_2
$$

**Torque-free** solutions conserve kinetic energy $T=\tfrac12\omega_iI_i\omega_i$ and angular
momentum magnitude. For `I₁=I₂` (symmetric top) they are closed-form — Nutation:

$$
\omega(t) = \left(A\cos(\Omega t), A\sin(\Omega t), \omega_3\right)
$$

**Inertia tensor** for point masses $m_a$ at $\mathbf r_a$:

$$
I_{ij} = \sum_a m_a\left(\|\mathbf r_a\|^2\delta_{ij} - r_{a,i}r_{a,j}\right)
$$

**Angular momentum** $\mathbf L = I\boldsymbol\omega$; under torque $\mathbf N$,
$\dot{\mathbf L} = \mathbf N$.

---

## Routines

| Routine | Description |
|---------|-------------|
| `hamiltonian_h(T, V)` | $H = T + V$ |
| `canonical_equations(H, q, p, h)` | $\dot q, \dot p$ pairs |
| `poisson_bracket(f, g, q, p, h)` | $\{f,g\}$ |
| `canonical_transform(...)` | symplectic-coordinate check |
| `lagrangian_l(T, V)` | $L = T - V$ |
| `euler_lagrange(L, q, qdot, h)` | EL-residual |
| `kinetic_energy(m, v)`, `potential_energy(...)` | $T$, $V$ |
| `euler_equations(I, omega, torque)` | $\dot{\boldsymbol\omega}$ |
| `torque_free(I, omega0, dt, steps)` | free precession |
| `inertia_tensor(masses, positions)` | $I_{ij}$ |
| `angular_momentum(I, omega)` | $\mathbf L$ |

---

## Usage Examples

### Harmonic oscillator via canonical equations

```python
from pysicrs import canonical_equations, hamiltonian_h

H = lambda q, p: 0.5 * p**2 + 0.5 * q**2   # H = T + V
# (q, p) = (0, 1) → q̇ = 1, ṗ = 0 (rest), physically dq/dt = ∂H/∂p = p
print(canonical_equations(H, 0.0, 1.0, 1e-5))  # (1.0, ~0.0)
```

### Free symmetric top

```python
from pysicrs import torque_free, inertia_tensor

I = [2.0, 2.0, 1.0]
omega0 = [1.0, 0.0, 2.0]
traj = torque_free(I, omega0, dt=0.01, steps=1000)
print(len(traj), traj[-1])   # bounded periodic nutation
```

### Angular momentum of a rotating body

```python
from pysicrs import inertia_tensor, angular_momentum

I = inertia_tensor(masses=[1.0, 1.0],
                   positions=[[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
L = angular_momentum(I, [0.0, 0.0, 2.0])
print(L)   # about ẑ
```

---

## Advantages & Limitations

✅ Both formulations exported — match your problem's natural coordinates

✅ Symplectic integrators (ODE) preserve the bracket structure for long times

✅ Rigid-body covers free tops, precession and torque-driven rotation

❌ No constraint mechanics (holonomic constraints not solved directly)

❌ Finite-difference derivatives need a well-chosen step

❌ No field-theoretic (continuum) classical mechanics yet

---

## References

1. Goldstein, H., Poole, C. & Safko, J. (2002). *Classical Mechanics*, 3rd ed. Addison-Wesley.
2. Landau, L.D. & Lifshitz, E.M. (1976). *Mechanics*, 3rd ed. Butterworth-Heinemann.
3. Arnold, V.I. (1989). *Mathematical Methods of Classical Mechanics*, 2nd ed. Springer.
4. Hairer, E., Lubich, C. & Wanner, G. (2006). *Geometric Numerical Integration*, 2nd ed. Springer.

---

## Related Topics

- [ODE](ode.md) – symplectic integrators are the natural companion for Hamiltonians
- [Linalg](linalg.md) – inertia tensors are symmetric matrices
- [General Relativity](general_relativity.md) – geodesic framework mirrors Lagrangian views