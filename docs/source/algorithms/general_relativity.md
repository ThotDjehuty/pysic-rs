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

## Mathematical Foundations

### Exact Metrics

**Schwarzschild** (exterior, geometric units $G=c=1$, mass $M$):

$$
ds^2 = -\left(1-\frac{2M}{r}\right)dt^2 + \left(1-\frac{2M}{r}\right)^{-1}dr^2
+ r^2\,d\Omega^2
$$

**Derivation:** The Schwarzschild metric is the unique static, spherically symmetric, asymptotically flat vacuum solution (Jebsen–Birkhoff theorem). Starting from the general static spherically symmetric metric:

$$
ds^2 = -e^{2\Phi(r)}dt^2 + e^{2\Lambda(r)}dr^2 + r^2 d\Omega^2
$$

The Einstein equations $R_{\mu\nu} = 0$ give $\Phi' = -\Lambda'$ and $r\Lambda' = 1 - e^{2\Lambda}$. Solving with boundary condition $\Phi \to 0$ as $r \to \infty$ yields $e^{2\Lambda} = (1 - 2M/r)^{-1}$.

Notable features: event horizon at $r=2M$, ISCO at $r=6M$, photon sphere at $r=3M$.

**Kerr** (rotating, Boyer–Lindquist coordinates, spin $a=J/M$):

$$
ds^2 = -\left(1 - \frac{2Mr}{\rho^2}\right)dt^2
- \frac{4Mar\sin^2\theta}{\rho^2}\,dt\,d\phi
+ \frac{\rho^2}{\Delta}dr^2 + \rho^2 d\theta^2
+ \frac{\sin^2\theta}{\rho^2}\left[\left(r^2+a^2\right)^2 - a^2\Delta\sin^2\theta\right]d\phi^2
$$

with $\rho^2 = r^2 + a^2\cos^2\theta$, $\Delta = r^2 - 2Mr + a^2$. The ergosphere lies between
the outer horizon and the static limit surface.

**FLRW** (spatially flat, $k=0$, scaled by $a(t)$):

$$
ds^2 = -dt^2 + a(t)^2\big(dx^2 + dy^2 + dz^2\big)
$$

**Minkowski:** $\eta_{\mu\nu} = \operatorname{diag}(-1, +1, +1, +1)$.

### Connection & Curvature

**Christoffel symbols** (Levi-Civita connection):

$$
\Gamma^{\rho}_{\mu\nu} = \frac{1}{2}g^{\rho\sigma}\left(\partial_\mu g_{\nu\sigma}
+ \partial_\nu g_{\sigma\mu} - \partial_\sigma g_{\mu\nu}\right)
$$

**Derivation:** The Christoffel symbols are determined by the metric compatibility condition $\nabla_\rho g_{\mu\nu} = 0$ and the torsion-free condition. The formula follows from the definition of the covariant derivative.

**Riemann tensor:**

$$
R^{\rho}{}_{\sigma\mu\nu} = \partial_\mu\Gamma^{\rho}_{\nu\sigma}
- \partial_\nu\Gamma^{\rho}_{\mu\sigma}
+ \Gamma^{\rho}_{\mu\lambda}\Gamma^{\lambda}_{\nu\sigma}
- \Gamma^{\rho}_{\nu\lambda}\Gamma^{\lambda}_{\mu\sigma}
$$

**Derivation:** The Riemann tensor measures the failure of covariant derivatives to commute: $[\nabla_\mu, \nabla_\nu]V^\rho = R^\rho{}_{\sigma\mu\nu}V^\sigma$. The formula follows from computing this commutator.

**Ricci tensor** $R_{\mu\nu} = R^{\lambda}{}_{\mu\lambda\nu}$, **Ricci scalar** $R = g^{\mu\nu}R_{\mu\nu}$,
**Einstein tensor:**

$$
G_{\mu\nu} = R_{\mu\nu} - \tfrac12 R\,g_{\mu\nu}
$$

**Proven**: $\nabla^\mu G_{\mu\nu} = 0$ (contracted Bianchi identities), so the Einstein field
equations $G_{\mu\nu} = 8\pi G\,T_{\mu\nu}$ are automatically consistent with stress-energy
conservation.

### Geodesics

The geodesic equation (auto-parallel transport):

$$
\frac{d^2x^\mu}{d\lambda^2} + \Gamma^{\mu}_{\rho\sigma}
\frac{dx^\rho}{d\lambda}\frac{dx^\sigma}{d\lambda} = 0
$$

**Derivation:** Geodesics extremize the proper time (or affine parameter). The Euler-Lagrange equations for the Lagrangian $L = \frac{1}{2}g_{\mu\nu}\dot{x}^\mu\dot{x}^\nu$ yield the geodesic equation.

Conventionally reduced to first order with $v^\mu = dx^\mu/d\lambda$:

$$
\frac{dv^\mu}{d\lambda} = -\Gamma^{\mu}_{\rho\sigma}\,v^\rho v^\sigma
$$

Integrated with the RK4/RK45 suite; for time-like geodesics affine parameter $\lambda = \tau$.

### ADM 3+1 Decomposition

The 4-metric is foliated with lapse $\alpha$, shift vector $\beta^i$, and 3-metric $\gamma_{ij}$:

$$
ds^2 = -(\alpha^2 - \beta_i\beta^i)dt^2 + 2\beta_i\,dx^i dt + \gamma_{ij}\,dx^i dx^j
$$

**Derivation:** The ADM decomposition splits spacetime into spatial hypersurfaces. The extrinsic curvature is $K_{ij} = \frac{1}{2\alpha}(\dot{\gamma}_{ij} - D_i\beta_j - D_j\beta_i)$.

**Hamiltonian constraint** (proven — encodes energy conservation of GR):

$$
\mathcal H = \frac{1}{2\kappa}\left(R^{(3)} - K_{ij}K^{ij} + K^2\right) = 0
$$

**Momentum constraint:**

$$
\mathcal M_i = D^jK_{ij} - D_iK = 0
$$

The `evolve_metric` routine applies Baumgarte–Shapiro–Shibata–Nakamura (BSSN)-style
evolution to the ADM variables.

### Stress–Energy Tensors

**Perfect fluid:**

$$
T_{\mu\nu} = (\rho + P)u_\mu u_\nu + P\,g_{\mu\nu}
$$

**Electromagnetic (in vacuum):**

$$
T_{\mu\nu} = F_{\mu\lambda}F_\nu{}^\lambda - \tfrac14 g_{\mu\nu} F_{\rho\sigma}F^{\rho\sigma}
$$

---

## Routines

| Routine | Returns | Derivation |
|---------|---------|------------|
| `schwarzschild_metric(r, M)` | $g_{\mu\nu}$ 4×4 | Einstein equations $R_{\mu\nu}=0$ |
| `kerr_metric(r, theta, M, a)` | $g_{\mu\nu}$ 4×4 | Rotating black hole solution |
| `flrw_metric(t, a)` | $g_{\mu\nu}$ 4×4 | Cosmological principle |
| `minkowski_metric()` | $\eta_{\mu\nu}$ | Flat spacetime |
| `christoffel_from_metric(g, coords, h)` | $\Gamma^{\rho}_{\mu\nu}$ | Metric compatibility |
| `riemann_tensor(g, coords, h)` | $R^{\rho}{}_{\sigma\mu\nu}$ | Covariant derivative commutator |
| `ricci_tensor(...)`, `ricci_scalar(...)` | $R_{\mu\nu}$, $R$ | Contraction of Riemann |
| `einstein_tensor(...)` | $G_{\mu\nu}$ | Contracted Bianchi identity |
| `geodesic_equation_rhs(g, x, v, h)` | $\ddot x^\mu$ | Euler-Lagrange equations |
| `integrate_geodesic(...)` | orbit $\{x^\mu(\lambda)\}$ | RK4/RK45 integration |
| `MetricSlice`, `AdmMetric`, `evolve_metric(...)` | ADM evolution | 3+1 decomposition |
| `hamiltonian_constraint(...)`, `momentum_constraint(...)` | $\mathcal H$, $\mathcal M_i$ | Einstein equations |
| `perfect_fluid_stress_energy(...)` | $T_{\mu\nu}$ | Energy-momentum tensor |
| `electromagnetic_stress_energy(...)` | $T^{\mu\nu}$ | Electromagnetic stress-energy |

---

## Usage Examples

### Schwarzschild metric at $r=10M$

```python
from pysicrs import schwarzschild_metric, metric_signature

g = schwarzschild_metric(10.0, M=1.0)
print(g[0][0])        # -(1 - 2/10) = -0.8
print(metric_signature(g))  # (1, 3)
```

### Time-like circular-orbit speed at the ISCO

```python
from pysicrs import schwarzschild_metric

# g_tt and g_φφ give the Keplerian orbit; for comparison the Newtonian v = sqrt(M/r)
for r in (100.0, 10.0):
    g = schwarzschild_metric(r, 1.0)
    print(f"r={r}: -g_tt = { -g[0][0]:.6f}")
```

### Hamiltonian constraint check on Schwarzschild

```python
from pysicrs import hamiltonian_constraint

H = hamiltonian_constraint(gamma_ij, K_ij, alpha, beta, D_i)
print(f"H ≈ {H:.2e}")   # vacuum constraint: should vanish
```

---

## Numerical Notes

- Derivatives use central finite differences with step `h` (default ~1e-5).
- Metrics are deliberately normalized to natural units $G=c=1$ to keep formulas clean;
  physical values enter via the stress–energy tensors.
- Geodesic integration is best run with the symplectic/RK suite — energy-conserving
  integrators matter for long orbits.

---

## Advantages & Limitations

✅ Exact metrics for all the classic spacetimes (including rotating Kerr)

✅ Full connection/tensor calculus and ADM constraint machinery in one crate

✅ Directly exposes constraint equations for numerical-relativity prototyping

❌ No adaptive mesh refinement; grids are uniform

❌ No dynamical (evolved) matter; $T_{\mu\nu}$ is prescribed

❌ Derivatives are finite-difference — loud noise if `h` is too small

---

## References

1. Wald, R. (1984). *General Relativity*. University of Chicago Press.
2. Carroll, S. (2004). *Spacetime and Geometry*. Addison-Wesley.
3. Misner, C., Thorne, K. & Wheeler, J. (1973). *Gravitation*. W.H. Freeman.
4. Alcubierre, M. (2008). *Introduction to 3+1 Numerical Relativity*. OUP.

---

## Related Topics

- [Linalg](linalg.md) – index raising/lowering uses `g`, `g⁻¹` from this module
- [ODE](ode.md) – geodesic integration, ADM evolution reuse the RK suite
- [Classical](classical.md) – Lagrangian viewpoint recurs in geometric mechanics