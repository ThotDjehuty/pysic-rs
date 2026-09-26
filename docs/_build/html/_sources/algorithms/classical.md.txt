# Classical Mechanics

Pysic-rs implements the two variational formulations of mechanics — Hamiltonian and
Lagrangian — plus rigid-body dynamics with the full suite of rotation parameterisations.
Everything here is textbook (proven) mechanics.

> **Python binding status:** `euler_equations`, `inertia_tensor`, and `euler_to_rotation` are
> bound. Hamiltonian/Lagrangian picture, Poisson brackets, and `torque_free` are Rust-only.
> See the [API page](../api/classical.md).

---

## 1. Lagrangian Mechanics

### 1.1 Configuration Space & Velocity

A system with $n$ generalised coordinates $q = (q_1,\ldots,q_n)$ has configuration space
$Q$ (a manifold). The velocity is $\dot{q} = (\dot{q}_1,\ldots,\dot{q}_n)$.

### 1.2 The Lagrangian

$$
L(q,\dot{q},t) = T - V
$$

where $T = \frac{1}{2}\sum_{ij}m_{ij}(q)\dot{q}_i\dot{q}_j$ is the kinetic energy
(may be non-diagonal in curvilinear coordinates) and $V(q,t)$ is the potential.

### 1.3 Hamilton's Principle

**Theorem (principle of least action).** The physical trajectory between fixed endpoints
$q(t_1)=q_1$, $q(t_2)=q_2$ extremises the **action**:

$$
S[q] = \int_{t_1}^{t_2}L(q,\dot{q},t)\,dt
$$

*Proof.* Consider a variation $q\to q+\delta q$ with $\delta q(t_1)=\delta q(t_2)=0$:

$$
\delta S = \int_{t_1}^{t_2}\left(\frac{\partial L}{\partial q_i}\delta q_i + \frac{\partial L}{\partial\dot{q}_i}\delta\dot{q}_i\right)dt
$$

Integrate the second term by parts:

$$
\delta S = \int_{t_1}^{t_2}\left[\frac{\partial L}{\partial q_i}-\frac{d}{dt}\frac{\partial L}{\partial\dot{q}_i}\right]\delta q_i\,dt + \left[\frac{\partial L}{\partial\dot{q}_i}\delta q_i\right]_{t_1}^{t_2}
$$

The boundary term vanishes. For $\delta S=0$ for all $\delta q_i$, each bracket must vanish:
the **Euler–Lagrange equation**. $\square$

### 1.4 Euler–Lagrange Equation

$$
\boxed{\frac{d}{dt}\frac{\partial L}{\partial\dot{q}_i} - \frac{\partial L}{\partial q_i} = 0, \qquad i=1,\ldots,n}
$$

### 1.5 Noether's Theorem

**Theorem (Noether, 1918).** Every continuous symmetry of the Lagrangian corresponds to a
conserved quantity.

*Proof.* Under a one-parameter transformation $q_i\to q_i+\varepsilon\,Q_i(q,\dot{q})$
with $\delta L = 0$:

$$
0 = \frac{\partial L}{\partial q_i}Q_i + \frac{\partial L}{\partial\dot{q}_i}\dot{Q}_i
= \frac{d}{dt}\left(\frac{\partial L}{\partial\dot{q}_i}Q_i\right)
$$

Hence $J = \frac{\partial L}{\partial\dot{q}_i}Q_i$ is conserved. $\square$

**Examples:**
- Time translation → energy conservation.
- Spatial translation → momentum conservation.
- Rotation → angular momentum conservation.

---

## 2. Hamiltonian Mechanics

### 2.1 Legendre Transform

The **canonical momenta** are

$$
p_i = \frac{\partial L}{\partial\dot{q}_i}
$$

The **Hamiltonian** is the Legendre transform of $L$:

$$
H(q,p,t) = \sum_i p_i\dot{q}_i - L(q,\dot{q},t)
$$

For natural systems ($T$ quadratic in $\dot{q}$, $V$ independent of $\dot{q}$):
$H = T + V$ = total energy.

### 2.2 Hamilton's Canonical Equations

$$
\boxed{\dot{q}_i = \frac{\partial H}{\partial p_i}, \qquad \dot{p}_i = -\frac{\partial H}{\partial q_i}}
$$

*Proof.* Differentiate $H = \sum p_i\dot{q}_i - L$:

$$
\frac{\partial H}{\partial p_i} = \dot{q}_i + \sum_j p_j\frac{\partial\dot{q}_j}{\partial p_i} - \sum_j\frac{\partial L}{\partial\dot{q}_j}\frac{\partial\dot{q}_j}{\partial p_i}
= \dot{q}_i + \sum_j\left(p_j-\frac{\partial L}{\partial\dot{q}_j}\right)\frac{\partial\dot{q}_j}{\partial p_i} = \dot{q}_i
$$

since $p_j = \partial L/\partial\dot{q}_j$. Similarly for $\partial H/\partial q_i$.
$\square$

### 2.3 Poisson Bracket

$$
\{f,g\} = \sum_i\left(\frac{\partial f}{\partial q_i}\frac{\partial g}{\partial p_i} - \frac{\partial f}{\partial p_i}\frac{\partial g}{\partial q_i}\right)
$$

**Properties (proven — Lie algebra structure):**
- Bilinearity, antisymmetry: $\{f,g\}=-\{g,f\}$.
- Jacobi identity: $\{f,\{g,h\}\}+\{g,\{h,f\}\}+\{h,\{f,g\}\}=0$.
- Leibniz rule: $\{fg,h\}=f\{g,h\}+g\{f,h\}$.
- Time evolution: $\frac{df}{dt}=\{f,H\}+\frac{\partial f}{\partial t}$.

An observable with $\{f,H\}=0$ (and $\partial f/\partial t=0$) is a constant of motion.

### 2.4 Canonical Transformations

A transformation $(q,p)\to(Q,P)$ is **canonical** if it preserves the Poisson bracket
structure:

$$
\{Q_i,P_j\}_{\rm new} = \delta_{ij}, \qquad \{Q_i,Q_j\}_{\rm new}=0, \qquad \{P_i,P_j\}_{\rm new}=0
$$

Equivalently (proven): the transformation preserves the **symplectic form**
$\Omega = \sum_i dq_i\wedge dp_i$.

---

## 3. Rigid Body Dynamics

### 3.1 Inertia Tensor

For point masses $m_a$ at positions $\mathbf{r}_a$:

$$
\boxed{I_{ij} = \sum_a m_a\left(\|\mathbf{r}_a\|^2\delta_{ij} - r_{a,i}r_{a,j}\right)}
$$

**Properties:**
- Symmetric: $I_{ij}=I_{ji}$.
- Positive semi-definite.
- Diagonalisable in the principal-axis frame: $I = \mathrm{diag}(I_1,I_2,I_3)$.

### 3.2 Euler's Equations

In the body frame (rotating with angular velocity $\boldsymbol{\omega}$), with principal
moments $I_1,I_2,I_3$:

$$
\boxed{I_1\dot{\omega}_1 = (I_2-I_3)\omega_2\omega_3 + N_1}
$$
$$
I_2\dot{\omega}_2 = (I_3-I_1)\omega_3\omega_1 + N_2
$$
$$
I_3\dot{\omega}_3 = (I_1-I_2)\omega_1\omega_2 + N_3
$$

*Proof.* From $\dot{\mathbf{L}} = \mathbf{N}$ in the lab frame, and
$\dot{\mathbf{L}}_{\rm body} = \mathbf{N} - \boldsymbol{\omega}\times\mathbf{L}$ (the
transport theorem). In the principal frame, $\mathbf{L} = (I_1\omega_1,I_2\omega_2,
I_3\omega_3)$, and the cross product gives the nonlinear coupling terms. $\square$

### 3.3 Torque-Free Motion

For $\mathbf{N}=0$, both kinetic energy $T = \frac{1}{2}\sum I_i\omega_i^2$ and angular
momentum magnitude $|\mathbf{L}|$ are conserved. The tip of $\boldsymbol{\omega}$ traces
an ellipse on the **inertia ellipsoid**.

**Symmetric top** ($I_1=I_2$): closed-form solution

$$
\omega(t) = \bigl(A\cos(\Omega t),\,A\sin(\Omega t),\,\omega_3\bigr), \qquad
\Omega = \frac{I_3-I_1}{I_1}\omega_3
$$

This is **free precession** — the angular velocity vector rotates around the symmetry axis.

---

## 4. Routines

| Routine | Description |
|---------|-------------|
| `hamiltonian_h(T, V)` | $H=T+V$ |
| `canonical_equations(H, q, p, h)` | $\dot{q},\dot{p}$ pairs |
| `poisson_bracket(f, g, q, p, h)` | $\{f,g\}$ |
| `canonical_transform(...)` | symplectic-coordinate check |
| `lagrangian_l(T, V)` | $L=T-V$ |
| `euler_lagrange(L, q, qdot, h)` | EL residual |
| `kinetic_energy(m, v)`, `potential_energy(...)` | $T$, $V$ |
| `euler_equations(I, omega, torque)` | $\dot{\boldsymbol{\omega}}$ |
| `torque_free(I, omega0, dt, steps)` | free precession |
| `inertia_tensor(masses, positions)` | $I_{ij}$ |
| `angular_momentum(I, omega)` | $\mathbf{L}$ |

---

## 5. Usage Examples

### Harmonic oscillator via canonical equations

```python
from pysicrs import canonical_equations, hamiltonian_h

H = lambda q, p: 0.5*p**2 + 0.5*q**2
print(canonical_equations(H, 0.0, 1.0, 1e-5))  # (1.0, ~0.0)
```

### Free symmetric top

```python
from pysicrs import torque_free

I = [2.0, 2.0, 1.0]
omega0 = [1.0, 0.0, 2.0]
traj = torque_free(I, omega0, dt=0.01, steps=1000)
print(len(traj), traj[-1])  # bounded periodic nutation
```

### Angular momentum of a rotating body

```python
from pysicrs import inertia_tensor, angular_momentum

I = inertia_tensor(masses=[1.0, 1.0],
                   positions=[[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
L = angular_momentum(I, [0.0, 0.0, 2.0])
print(L)  # about ẑ
```

---

## 6. Advantages & Limitations

✅ Both formulations exported — match your problem's natural coordinates

✅ Symplectic integrators (ODE) preserve the bracket structure for long times

✅ Rigid-body covers free tops, precession and torque-driven rotation

❌ No constraint mechanics (holonomic constraints not solved directly)

❌ Finite-difference derivatives need a well-chosen step

❌ No field-theoretic (continuum) classical mechanics yet

---

## 7. References

1. Goldstein, H., Poole, C. & Safko, J. (2002). *Classical Mechanics*, 3rd ed. Addison-Wesley.
2. Landau, L.D. & Lifshitz, E.M. (1976). *Mechanics*, 3rd ed. Butterworth-Heinemann.
3. Arnold, V.I. (1989). *Mathematical Methods of Classical Mechanics*, 2nd ed. Springer.
4. Hairer, E., Lubich, C. & Wanner, G. (2006). *Geometric Numerical Integration*, 2nd ed. Springer.

---

## 8. Related Topics

- [ODE](ode.md) — symplectic integrators are the natural companion for Hamiltonians
- [Linear Algebra](linalg.md) — inertia tensors are symmetric matrices
- [General Relativity](general_relativity.md) — geodesic framework mirrors Lagrangian views
