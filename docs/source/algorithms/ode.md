# ODE Solvers

Pysic-rs provides the standard initial-value-problem (IVP) solvers used throughout physics:
explicit Runge–Kutta methods, implicit and semi-implicit schemes, and structure-preserving
*symplectic* integrators for Hamiltonian systems.

> **Python binding status:** `rk4_solve`, `rk45_solve`, `backward_euler`, and
> `leapfrog_integrate` are bound. `integrate_hamiltonian` is exposed from Python as
> `leapfrog_integrate`. See the [API page](../api/ode.md) for exact signatures; the examples
> below target the Rust API where noted.

---

## Mathematical Foundations

We solve the Cauchy problem for a first-order system

$$
\dot{y} = f(t, y), \qquad y(t_0) = y_0
$$

Any higher-order ODE is reduced to this form by introducing derivatives as new variables.
A one-step method advances $y_n \to y_{n+1}$ with step $h$.

### Runge–Kutta 4 (RK4)

The classical 4-stage method with error constant $O(h^4)$ per step:

$$
k_1 = f(t_n, y_n)
$$
$$
k_2 = f\!\left(t_n + \tfrac{h}{2}, y_n + \tfrac{h}{2} k_1\right)
$$
$$
k_3 = f\!\left(t_n + \tfrac{h}{2}, y_n + \tfrac{h}{2} k_2\right)
$$
$$
k_4 = f(t_n + h, y_n + h k_3)
$$
$$
y_{n+1} = y_n + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)
$$

**Properties (proven):** consistent of order 4 (matches the Taylor expansion through $h^4$),
explicit, zero-stable, hence **convergent** of order 4.

### Runge–Kutta 5(4) — Dormand–Prince

The Dormand–Prince pair (the method behind `scipy.integrate.solve_ivp`, `rk45`) embeds a
5th-order solution with a 4th-order error estimator in 7 stages:

$$
y_{n+1} = y_n + h\sum_i b_i k_i, \qquad \hat{y}_{n+1} = y_n + h\sum_i \hat b_i k_i
$$

The difference $\|y_{n+1} - \hat y_{n+1}\|$ estimates the local error and drives adaptive
step-size control (PI controller), targeting a user-supplied `rtol`/`atol`:

$$
h_{\text{new}} = h\left(\frac{\text{tol}}{\|e\|}\right)^{1/5}
$$

### Implicit Backward Euler

$$
y_{n+1} = y_n + h f(t_{n+1}, y_{n+1})
$$

Solved by fixed-point (Picard) iteration. **Proven A-stable** — the amplification factor is
$R(z) = (1-z)^{-1}$, whose magnitude is $< 1$ for all $\operatorname{Re} z < 0$. This makes it
unconditionally stable for stiff problems (at $O(h)$ accuracy).

### Crank–Nicolson (semi-implicit trapezoid)

$$
y_{n+1} = y_n + \frac{h}{2}\big(f(t_n, y_n) + f(t_{n+1}, y_{n+1})\big)
$$

**Proven A-stable and $O(h^2)$:** $R(z) = \dfrac{1+z/2}{1-z/2}$ lies on the unit circle in the
left half-plane. This is the ODE analogue of the trapezoidal rule used by the PDE heat solver.

### Symplectic Integrators

For separable Hamiltonians $H(q,p) = T(p) + V(q)$, symplectic integrators preserve the phase-
space volume form (Liouville), so **energy error stays bounded** over long times — unlike RK4,
whose energy error drifts monotonically.

**Leapfrog / kick-drift-kick:**

$$
p_{n+1/2} = p_n - \frac{h}{2}\nabla V(q_n), \qquad
q_{n+1} = q_n + h\, M^{-1}p_{n+1/2}, \qquad
p_{n+1} = p_{n+1/2} - \frac{h}{2}\nabla V(q_{n+1})
$$

**Velocity Verlet** is the same map written with explicit kinetic updates; both are second-
order accurate and exactly symplectic.

**Yoshida's fourth-order construction:** composing the second-order leapfrog map
$\Phi_h$ as $\Phi_{\lambda h}\,\Phi_{\mu h}\,\Phi_{\lambda h}$ with
$\lambda = (2 - 2^{1/3})^{-1}$, $\mu = 1 - 2\lambda$ cancels the leading error term,
yielding a 4th-order symplectic integrator (Yoshida 1990).

---

## Routines

| Routine | Order | Type | Best for |
|---------|-------|------|----------|
| `rk4_solve` | 4 | explicit | smooth, non-stiff systems, fixed step |
| `rk45_solve` | 5(4) adaptive | explicit | autonomous/non-autonomous, tolerance control |
| `backward_euler` | 1 | implicit A-stable | stiff systems |
| `crank_nicolson_ode` | 2 | implicit A-stable | moderately stiff, more accurate |
| `leapfrog_step` / `velocity_verlet_step` | 2 | symplectic | Hamiltonian dynamics, long-time orbits |
| `yoshida_step` | 4 | symplectic | high-accuracy Hamiltonian integration |
| `integrate_hamiltonian` | 2 | symplectic (Verlet family) | separable $H$; returns $(q_t, p_t)$ paths |

---

## Usage Examples

### Harmonic oscillator with RK4

```python
from pysicrs import rk4_solve

def osc(t, y):
    return [y[1], -y[0]]

times, traj = rk4_solve(osc, [1.0, 0.0], (0.0, 10.0), 10000)
print(traj[-1][0])        # ≈ cos(10) = -0.8391
```

### Long-time orbit with a symplectic integrator (energy stays bounded)

```python
from pysicrs import integrate_hamiltonian

# H = p²/2 + q²/2   ->  ∇V = q,  ∇T = p
grad_v = lambda q: [q[0]]
grad_t = lambda p: [p[0]]

q, p = integrate_hamiltonian([0.0], [1.0], grad_v, grad_t, 0.01, 10000, inv_mass=1.0)

E0 = 0.5          # initial energy E = p²/2 + q²/2
Ef = 0.5 * p[-1][0]**2 + 0.5 * q[-1][0]**2
print(f"energy drift after 10⁴ steps: {abs(Ef - E0):.2e}")  # tiny, bounded
```

### Stiff decay with backward Euler

```python
from pysicrs import backward_euler

def stiff(t, y):
    return [-1000.0 * y[0]]

times, y = backward_euler(stiff, [1.0], (0.0, 0.1), 100)
print(y[-1][0])   # no blow-up, stable
```

---

## Accuracy & Performance

| Method | Global error (smooth problems) | Cost/step | Stability |
|--------|-------------------------------|-----------|-----------|
| RK4 | $O(h^4)$ | 4 `f` evals | conditionally stable (non-stiff OK) |
| Dormand–Prince | adaptive to `rtol`/`atol` | 7 stages | adaptive stepping |
| Backward Euler | $O(h)$ | Newton/fixed-point | A-stable |
| Crank–Nicolson | $O(h^2)$ | implicit solve | A-stable |
| Verlet/leapfrog | $O(h^2)$ | 1 force eval | symplectic (bounded energy) |
| Yoshida | $O(h^4)$ | 3 force evals | symplectic |

---

## Advantages & Limitations

✅ Explicit methods trivial to use; symplectic integrators ideal for long Hamiltonian runs

✅ A-stable implicits handle stiff reaction kinetics

✅ Adaptivity in RK45 frees you from guessing step sizes

❌ Symplectic methods need a separable Hamiltonian

❌ RK4 has no built-in error control — choose `h` conservatively

❌ Backward Euler is only first order; use Crank–Nicolson for accuracy on stiff problems

---

## References

1. Hairer, E., Nørsett, S.P. & Wanner, G. (1993). *Solving Ordinary Differential Equations I*, 2nd ed. Springer.
2. Hairer, E., Lubich, C. & Wanner, G. (2006). *Geometric Numerical Integration*, 2nd ed. Springer.
3. Yoshida, H. (1990). "Construction of higher order symplectic integrators." *Phys. Lett. A* 150(5–7):262–268.
4. Dormand, J.R. & Prince, P.J. (1980). "A family of embedded Runge–Kutta formulae." *J. Comput. Appl. Math.* 6(1):19–26.

---

## Related Topics

- [PDE](pde.md) – the heat and wave solvers reuse Crank–Nicolson / FDTD ideas
- [Classical](classical.md) – supplies the `H`, `∇V`, `∇T` needed by symplectic integrators
- [Calculus](calculus.md) – quadrature kernels shared with integrators