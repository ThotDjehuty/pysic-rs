# General Relativity

Pysic-rs implements the standard machinery of 4D Lorentzian geometry and the Einstein
equations — exact metrics, connection/tensor calculus, geodesic evolution, the ADM
3+1 decomposition, and stress–energy tensors — using only textbook results (Wald,
Carroll, Misner–Thorne–Wheeler).

> **Python binding status:** `schwarzschild_metric`, `christoffel_from_metric`, `ricci_scalar`,
> `einstein_tensor`, `hamiltonian_constraint`, and `momentum_constraint` are bound. The other
> metrics, Riemann/Ricci tensors, geodesics, ADM evolution, and stress–energy are Rust-only.
> See the [API page](../api/general_relativity.md).

---

## 1. Manifolds & the Metric Tensor

### 1.1 Spacetime as a Manifold

General relativity describes gravity as the curvature of a 4-dimensional pseudo-Riemannian
manifold $(\mathcal{M}, g_{\mu\nu})$. The metric tensor $g_{\mu\nu}$ defines:

- **Lengths:** $ds^2 = g_{\mu\nu}\,dx^\mu\,dx^\nu$
- **Angles:** via the inner product $g(u,v) = g_{\mu\nu}u^\mu v^\nu$
- **Causal structure:** $g(u,u)<0$ (timelike), $=0$ (null), $>0$ (spacelike)

**Convention:** signature $(-,+,+,+)$ throughout (the "particle physics" or "mostly plus"
convention).

### 1.2 Levi-Civita Connection

The unique torsion-free, metric-compatible connection is the **Levi-Civita connection**,
specified by the **Koszul formula**:

$$
2g(\nabla_X Y, Z) = X\,g(Y,Z) + Y\,g(X,Z) - Z\,g(X,Y)
+ g([X,Y],Z) - g([X,Z],Y) - g([Y,Z],X)
$$

This determines the Christoffel symbols completely (see §2.1).

---

## 2. Connection & Curvature

### 2.1 Christoffel Symbols

**Definition.** The Christoffel symbols of the second kind are defined by

$$
\nabla_\mu e_\nu = \Gamma^\rho_{\mu\nu}\,e_\rho
$$

in a coordinate basis $\{e_\mu = \partial_\mu\}$.

**Theorem (Levi-Civita formula).**

$$
\boxed{\Gamma^\rho_{\mu\nu} = \frac{1}{2}g^{\rho\sigma}\left(\partial_\mu g_{\nu\sigma}
+ \partial_\nu g_{\sigma\mu} - \partial_\sigma g_{\mu\nu}\right)}
$$

*Proof.* From the Koszul formula with $X=\partial_\mu$, $Y=\partial_\nu$, $Z=\partial_\sigma$:
all Lie brackets $[\partial_\mu,\partial_\nu]=0$ in coordinate bases, and $g(\partial_\mu,\partial_\nu)
= g_{\mu\nu}$. The three terms involving metric derivatives yield the formula above after
contracting with $g^{\rho\sigma}$. $\square$

**Key properties:**
- $\Gamma^\rho_{\mu\nu} = \Gamma^\rho_{\nu\mu}$ (symmetric = torsion-free).
- NOT a tensor — transforms inhomogeneously under coordinate changes.
- Vanishes in locally inertial (freely falling) coordinates at any given point.

### 2.2 Riemann Curvature Tensor

**Definition.** The Riemann tensor measures the failure of covariant derivatives to commute:

$$
\boxed{R^\rho{}_{\sigma\mu\nu} = \partial_\mu\Gamma^\rho_{\nu\sigma}
- \partial_\nu\Gamma^\rho_{\mu\sigma}
+ \Gamma^\rho_{\mu\lambda}\Gamma^\lambda_{\nu\sigma}
- \Gamma^\rho_{\nu\lambda}\Gamma^\lambda_{\mu\sigma}}
$$

*Proof.* Apply $[\nabla_\mu,\nabla_\nu]$ to a vector $V^\rho$:

$$
[\nabla_\mu,\nabla_\nu]V^\rho = R^\rho{}_{\sigma\mu\nu}V^\sigma
+ T^\lambda_{\mu\nu}\nabla_\lambda V^\rho
$$

For torsion-free ($T=0$), the second term vanishes. Expanding the covariant derivatives
and using the symmetry of the Christoffel symbols gives the formula. $\square$

**Symmetries (proven):**

1. $R_{\rho\sigma\mu\nu} = -R_{\sigma\rho\mu\nu}$ (antisymmetric in first pair)
2. $R_{\rho\sigma\mu\nu} = -R_{\rho\sigma\nu\mu}$ (antisymmetric in second pair)
3. $R_{\rho\sigma\mu\nu} = R_{\mu\nu\rho\sigma}$ (pair symmetry)
4. $R_{\rho[\sigma\mu\nu]} = 0$ (first Bianchi identity)

These reduce the independent components from $4^4 = 256$ to **20** in 4D.

**Flatness criterion:** $R^\rho{}_{\sigma\mu\nu} = 0$ everywhere $\Longleftrightarrow$
the manifold is locally isometric to Minkowski space.

### 2.3 Ricci Tensor & Scalar

**Contraction of the Riemann tensor:**

$$
R_{\mu\nu} = R^\lambda{}_{\mu\lambda\nu}, \qquad R = g^{\mu\nu}R_{\mu\nu}
$$

The Ricci tensor is symmetric ($R_{\mu\nu}=R_{\nu\mu}$). The Ricci scalar is the
trace of curvature.

### 2.4 Einstein Tensor

$$
\boxed{G_{\mu\nu} = R_{\mu\nu} - \tfrac{1}{2}R\,g_{\mu\nu}}
$$

**Theorem (contracted Bianchi identity).** $\nabla^\mu G_{\mu\nu} = 0$.

*Proof.* The second Bianchi identity states $\nabla_{[\lambda}R_{\rho\sigma]\mu\nu}=0$.
Contracting $\lambda$ with $\mu$ and using the symmetries of the Riemann tensor yields
$\nabla^\mu G_{\mu\nu}=0$. This is the mathematical identity underlying local energy-momentum
conservation in GR. $\square$

---

## 3. Einstein Field Equations

$$
\boxed{G_{\mu\nu} = 8\pi G\,T_{\mu\nu} + \Lambda\,g_{\mu\nu}}
$$

where $\Lambda$ is the cosmological constant and $T_{\mu\nu}$ is the stress–energy tensor.
In geometric units $G=c=1$: $G_{\mu\nu} = 8\pi T_{\mu\nu} + \Lambda g_{\mu\nu}$.

---

## 4. Exact Metrics

### 4.1 Schwarzschild Metric

**Theorem (Jebsen–Birkhoff).** The unique static, spherically symmetric, vacuum solution
of Einstein's equations is the **Schwarzschild metric**:

$$
\boxed{ds^2 = -\left(1-\frac{2M}{r}\right)dt^2 + \left(1-\frac{2M}{r}\right)^{-1}dr^2
+ r^2\,d\Omega^2}
$$

where $d\Omega^2 = d\theta^2 + \sin^2\theta\,d\phi^2$.

*Proof.* Start from the general static spherically symmetric ansatz:

$$
ds^2 = -e^{2\Phi(r)}dt^2 + e^{2\Lambda(r)}dr^2 + r^2\,d\Omega^2
$$

The vacuum Einstein equations $R_{\mu\nu}=0$ give (from $R_{tt}$ and $R_{rr}$):

$$
\Phi' = -\Lambda', \qquad r\Lambda' = 1-e^{2\Lambda}
$$

These combine to $\frac{d}{dr}[r(1-e^{-2\Lambda})]=0$, giving $e^{-2\Lambda}=1-2M/r$ for
some constant $M$ (identified as the mass). The boundary condition $\Phi\to 0$ as
$r\to\infty$ fixes $\Phi = -\Lambda = \frac{1}{2}\ln(1-2M/r)$. $\square$

**Notable features:**

- **Event horizon:** $r = 2M$ (Schwarzschild radius)
- **Photon sphere:** $r = 3M$
- **Innermost stable circular orbit (ISCO):** $r = 6M$
- **Gravitational radius:** $r_g = 2GM/c^2$

### 4.2 Kerr Metric (Rotating Black Hole)

In Boyer–Lindquist coordinates with spin parameter $a = J/M$:

$$
ds^2 = -\left(1-\frac{2Mr}{\rho^2}\right)dt^2
- \frac{4Mar\sin^2\theta}{\rho^2}\,dt\,d\phi
+ \frac{\rho^2}{\Delta}dr^2 + \rho^2\,d\theta^2
+ \frac{\sin^2\theta}{\rho^2}\left[(r^2+a^2)^2-a^2\Delta\sin^2\theta\right]d\phi^2
$$

with $\rho^2 = r^2+a^2\cos^2\theta$, $\Delta = r^2-2Mr+a^2$.

The **ergosphere** lies between the outer horizon $r_+ = M+\sqrt{M^2-a^2}$ and the static
limit surface $r = M+\sqrt{M^2-a^2\cos^2\theta}$. Inside the ergosphere, no observer can
remain stationary — frame dragging is inescapable.

### 4.3 FLRW Metric (Cosmology)

The Friedmann–Lemaître–Robertson–Walker metric for a homogeneous, isotropic universe:

$$
\boxed{ds^2 = -dt^2 + a(t)^2\left[\frac{dr^2}{1-kr^2}+r^2\,d\Omega^2\right]}
$$

with $k\in\{-1,0,+1\}$ (open, flat, closed). The scale factor $a(t)$ satisfies the
Friedmann equations:

$$
\left(\frac{\dot{a}}{a}\right)^2 = \frac{8\pi G}{3}\rho - \frac{k}{a^2}+\frac{\Lambda}{3}
$$

$$
\frac{\ddot{a}}{a} = -\frac{4\pi G}{3}(\rho+3p)+\frac{\Lambda}{3}
$$

### 4.4 Minkowski Metric

$$
\eta_{\mu\nu} = \operatorname{diag}(-1,+1,+1,+1)
$$

The flat-space limit. All curvature tensors vanish identically.

---

## 5. Geodesics

### 5.1 Geodesic Equation

The geodesic equation (auto-parallel transport) is

$$
\boxed{\frac{d^2x^\mu}{d\lambda^2}+\Gamma^\mu_{\rho\sigma}\frac{dx^\rho}{d\lambda}\frac{dx^\sigma}{d\lambda} = 0}
$$

*Proof.* Geodesics extremise the proper time (or affine parameter). The Lagrangian is
$L = \frac{1}{2}g_{\mu\nu}\dot{x}^\mu\dot{x}^\nu$. The Euler–Lagrange equation gives:

$$
\frac{d}{d\lambda}(g_{\mu\nu}\dot{x}^\nu) - \frac{1}{2}\partial_\mu g_{\rho\sigma}\dot{x}^\rho\dot{x}^\sigma = 0
$$

Expanding the first term and using $\partial_\mu g_{\nu\sigma} = \Gamma_{\mu\nu\sigma}+
\Gamma_{\mu\sigma\nu}$ (metric compatibility) yields the geodesic equation. $\square$

**First-order form** (used in numerical integration):

$$
\frac{dv^\mu}{d\lambda} = -\Gamma^\mu_{\rho\sigma}\,v^\rho v^\sigma, \qquad v^\mu = \frac{dx^\mu}{d\lambda}
$$

### 5.2 Conserved Quantities

For metric symmetries (Killing vectors $\xi^\mu$), the quantity $u_\mu\xi^\mu$ is conserved
along geodesics:

- **Time translation** $\xi^\mu = (1,0,0,0)$: conservation of energy $E = (1-2M/r)\dot{t}$.
- **Rotational symmetry** $\xi^\mu = (0,0,0,1)$: conservation of angular momentum $L = r^2\sin^2\theta\,\dot{\phi}$.

### 5.3 Effective Potential (Schwarzschild)

For time-like geodesics in the equatorial plane ($\theta=\pi/2$), the radial equation
reduces to

$$
\frac{1}{2}\dot{r}^2 + V_{\rm eff}(r) = \frac{1}{2}E^2, \qquad
V_{\rm eff}(r) = \frac{1}{2}\left(1-\frac{2M}{r}\right)\left(1+\frac{L^2}{r^2}\right)
$$

This effective potential has a maximum (unstable circular orbit) at $r=3M$ (photon sphere)
and a minimum (stable circular orbit) at $r=6M$ (ISCO).

---

## 6. ADM 3+1 Decomposition

### 6.1 Foliation of Spacetime

The 4-metric is decomposed into spatial hypersurfaces $\Sigma_t$ with:

- **Lapse function** $\alpha$: rate of proper time advance between slices
- **Shift vector** $\beta^i$: spatial coordinate dragging between slices
- **3-metric** $\gamma_{ij}$: intrinsic geometry of each slice

$$
\boxed{ds^2 = -(\alpha^2-\beta_i\beta^i)dt^2 + 2\beta_i\,dx^i\,dt + \gamma_{ij}\,dx^i\,dx^j}
$$

### 6.2 Extrinsic Curvature

The extrinsic curvature of $\Sigma_t$ embedded in spacetime is

$$
K_{ij} = \frac{1}{2\alpha}\left(\dot{\gamma}_{ij}-D_i\beta_j-D_j\beta_i\right)
$$

where $D_i$ is the covariant derivative compatible with $\gamma_{ij}$.

$K_{ij}$ measures how the spatial geometry changes from one slice to the next — it is the
"time derivative of the 3-metric" in a geometrically invariant sense.

### 6.3 Constraint Equations

The Einstein equations decompose into **constraints** (on each slice) and **evolution**
(across slices).

**Hamiltonian constraint** (energy):

$$
\boxed{\mathcal{H} = \frac{1}{2\kappa}\left(R^{(3)}-K_{ij}K^{ij}+K^2\right) = 0}
$$

where $R^{(3)}$ is the 3 Ricci scalar, $K=K^i{}_i$ is the trace, and $\kappa=8\pi G$.

**Momentum constraint** (momentum):

$$
\boxed{\mathcal{M}_i = D^jK_{ij}-D_iK = 0}
$$

*Proof.* These follow from the contracted Bianchi identity applied to the foliation.
They are not evolution equations — they must be satisfied on every initial data slice.
They constrain the allowed initial data for the evolution equations. $\square$

### 6.4 BSSN Evolution

The **Baumgarte–Shapiro–Shibata–Nakamura (BSSN)** formulation reformulates the ADM
evolution equations in a way that is numerically stable. Key modifications:

1. Conformal decomposition: $\tilde{\gamma}_{ij} = e^{-4\phi}\gamma_{ij}$ with
   $\det(\tilde{\gamma})=1$.
2. Explicit trace-free extrinsic curvature: $\tilde{A}_{ij} = e^{-4\phi}(K_{ij}-\frac{1}{3}K\gamma_{ij})$.
3.引入 conformal connection functions $\tilde{\Gamma}^i = \tilde{\gamma}^{jk}\tilde{\Gamma}^i_{jk}$.

The `evolve_metric` routine applies BSSN-style evolution to the ADM variables.

---

## 7. Stress–Energy Tensors

### 7.1 Perfect Fluid

$$
T_{\mu\nu} = (\rho+P)u_\mu u_\nu + P\,g_{\mu\nu}
$$

where $\rho$ is the energy density, $P$ the pressure, and $u^\mu$ the 4-velocity
($u^\mu u_\mu = -1$).

**Conservation:** $\nabla^\mu T_{\mu\nu}=0$ gives the continuity equation and Euler
equation.

### 7.2 Electromagnetic Field

$$
T_{\mu\nu} = F_{\mu\lambda}F_\nu{}^\lambda - \tfrac{1}{4}g_{\mu\nu}F_{\rho\sigma}F^{\rho\sigma}
$$

where $F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu$ is the electromagnetic field
tensor (see [Gauge](gauge.md)).

### 7.3 Scalar Field

$$
T_{\mu\nu} = \partial_\mu\phi\,\partial_\nu\phi - g_{\mu\nu}\left(\tfrac{1}{2}g^{\rho\sigma}\partial_\rho\phi\,\partial_\sigma\phi + V(\phi)\right)
$$

---

## 8. Routines

| Routine | Returns | Derivation |
|---------|---------|------------|
| `schwarzschild_metric(r, M)` | $g_{\mu\nu}$ 4×4 | Einstein equations $R_{\mu\nu}=0$ |
| `kerr_metric(r, theta, M, a)` | $g_{\mu\nu}$ 4×4 | Rotating black hole solution |
| `flrw_metric(t, a)` | $g_{\mu\nu}$ 4×4 | Cosmological principle |
| `minkowski_metric()` | $\eta_{\mu\nu}$ | Flat spacetime |
| `christoffel_from_metric(g, coords, h)` | $\Gamma^\rho_{\mu\nu}$ | Metric compatibility |
| `riemann_tensor(g, coords, h)` | $R^\rho{}_{\sigma\mu\nu}$ | Covariant derivative commutator |
| `ricci_tensor(...)` | $R_{\mu\nu}$ | Contraction of Riemann |
| `ricci_scalar(...)` | $R$ | Trace of Ricci |
| `einstein_tensor(...)` | $G_{\mu\nu}$ | Contracted Bianchi identity |
| `geodesic_equation_rhs(g, x, v, h)` | $\ddot{x}^\mu$ | Euler–Lagrange equations |
| `integrate_geodesic(...)` | orbit $\{x^\mu(\lambda)\}$ | RK4/RK45 integration |
| `MetricSlice`, `AdmMetric`, `evolve_metric(...)` | ADM evolution | 3+1 decomposition |
| `hamiltonian_constraint(...)` | $\mathcal{H}$ | Einstein equations |
| `momentum_constraint(...)` | $\mathcal{M}_i$ | Einstein equations |
| `perfect_fluid_stress_energy(...)` | $T_{\mu\nu}$ | Energy-momentum tensor |
| `electromagnetic_stress_energy(...)` | $T^{\mu\nu}$ | Electromagnetic stress-energy |

---

## 9. Usage Examples

### Schwarzschild metric at $r=10M$

```python
from pysicrs import schwarzschild_metric, metric_signature

g = schwarzschild_metric(10.0, M=1.0)
print(g[0][0])           # -(1 - 2/10) = -0.8
print(metric_signature(g))  # (1, 3) → Lorentzian
```

### Time-like circular-orbit energy

```python
from pysicrs import schwarzschild_metric

for r in (100.0, 10.0, 6.0):
    g = schwarzschild_metric(r, 1.0)
    print(f"r={r:.1f}M: -g_tt = {-g[0][0]:.6f}")
```

### Hamiltonian constraint check on Schwarzschild

```python
from pysicrs import hamiltonian_constraint

H = hamiltonian_constraint(gamma_ij, K_ij, alpha, beta, D_i)
print(f"H ≈ {H:.2e}")   # vacuum constraint: should vanish
```

---

## 10. Numerical Notes

- Derivatives use central finite differences with step $h$ (default $\sim 10^{-5}$).
- Metrics are normalised to natural units $G=c=1$; physical values enter via stress–energy.
- Geodesic integration requires symplectic or RK4/RK45 — energy-conserving integrators
  matter for long orbits (see [ODE](ode.md)).
- BSSN evolution is numerically stable but requires the conformal Decomposition to be
  constraint-preserving.

---

## 11. Advantages & Limitations

✅ Exact metrics for all classic spacetimes (including rotating Kerr)

✅ Full connection/tensor calculus and ADM constraint machinery in one crate

✅ Directly exposes constraint equations for numerical-relativity prototyping

❌ No adaptive mesh refinement; grids are uniform

❌ No dynamical (evolved) matter; $T_{\mu\nu}$ is prescribed

❌ Derivatives are finite-difference — noisy if $h$ is too small

---

## 12. References

1. Wald, R. (1984). *General Relativity*. University of Chicago Press.
2. Carroll, S. (2004). *Spacetime and Geometry*. Addison-Wesley.
3. Misner, C., Thorne, K. & Wheeler, J. (1973). *Gravitation*. W.H. Freeman.
4. Alcubierre, M. (2008). *Introduction to 3+1 Numerical Relativity*. OUP.
5. Baumgarte, T. & Shapiro, S. (2010). *Numerical Relativity*. Cambridge.
6. Baumgarte, T. & Shapiro, S. (1998). "On the numerical integration of Einstein's field equations." *Phys. Rev. D* 59:024007.

---

## 13. Related Topics

- [Linear Algebra](linalg.md) — index raising/lowering uses $g$, $g^{-1}$
- [ODE](ode.md) — geodesic integration, ADM evolution reuse the RK suite
- [Classical](classical.md) — Lagrangian viewpoint recurs in geometric mechanics
- [Gauge](gauge.md) — same covariant-derivative & curvature pattern
- [Special Functions](special_functions.md) — $K_1$ in Euclidean propagators; spherical harmonics in perturbation theory
