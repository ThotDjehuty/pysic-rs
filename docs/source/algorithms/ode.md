# ODE Solvers

Pysic-rs provides the standard initial-value-problem (IVP) solvers used throughout physics:
explicit Runge–Kutta methods, implicit and semi-implicit schemes, and structure-preserving
*symplectic* integrators for Hamiltonian systems.

> **Python binding status:** `rk4_solve`, `rk45_solve`, `backward_euler`, and
> `leapfrog_integrate` are bound. `integrate_hamiltonian` is exposed from Python as
> `leapfrog_integrate`. See the [API page](../api/ode.md).

---

## 1. The Initial-Value Problem

We solve the Cauchy problem for a first-order system:

$$
\dot{y} = f(t,y), \qquad y(t_0) = y_0
$$

Any higher-order ODE is reduced to this form by introducing derivatives as new variables.
A **one-step method** advances $y_n\to y_{n+1}$ with step $h$.

**Convergence framework.** A method is **convergent** if the global error $\|y(t_n)-y_n\|
\to 0$ as $h\to 0$. By the Dahlquist equivalence theorem, convergence of a zero-stable
linear multistep method is equivalent to consistency + zero stability.

---

## 2. Explicit Runge–Kutta Methods

### 2.1 Runge–Kutta 4 (RK4)

The classical 4-stage method with local truncation error $\mathcal{O}(h^5)$:

$$
\begin{aligned}
k_1 &= f(t_n, y_n) \\
k_2 &= f\!\left(t_n+\tfrac{h}{2}, y_n+\tfrac{h}{2}k_1\right) \\
k_3 &= f\!\left(t_n+\tfrac{h}{2}, y_n+\tfrac{h}{2}k_2\right) \\
k_4 &= f(t_n+h, y_n+hk_3) \\
y_{n+1} &= y_n + \frac{h}{6}(k_1+2k_2+2k_3+k_4)
\end{aligned}
$$

**Properties (proven):**
- Consistent of order 4 (matches the Taylor expansion through $h^4$).
- Explicit, zero-stable, hence **convergent of order 4**.
- Global error: $\|e_n\| = \mathcal{O}(h^4)$.
- Cost: 4 function evaluations per step.

### 2.2 Runge–Kutta 5(4) — Dormand–Prince

The Dormand–Prince pair (the method behind `scipy.integrate.solve_ivp`) embeds a 5th-order
solution with a 4th-order error estimator in 7 stages:

$$
y_{n+1} = y_n + h\sum_{i=1}^7 b_i k_i, \qquad \hat{y}_{n+1} = y_n + h\sum_{i=1}^7 \hat{b}_i k_i
$$

The difference $\|y_{n+1}-\hat{y}_{n+1}\|$ estimates the local error and drives adaptive
step-size control (PI controller), targeting user-supplied `rtol`/`atol`:

$$
h_{\rm new} = h\left(\frac{\rm tol}{\|e\|}\right)^{1/5}
$$

The coefficients $\{b_i\}$ and $\{\hat{b}_i\}$ are chosen so that the leading error term
is minimised and the 4th-order solution is **FSAL** (first-same-as-last) — the last
function evaluation of one step is reused as the first of the next.

---

## 3. Implicit Methods

### 3.1 Backward Euler

$$
y_{n+1} = y_n + hf(t_{n+1}, y_{n+1})
$$

Solved by fixed-point (Picard) iteration: $y^{(0)}_{n+1}=y_n$, $y^{(k+1)}_{n+1}=y_n+
hf(t_{n+1},y^{(k)}_{n+1})$. Converges if $h\|{\partial f}/{\partial y}\| < 1$.

**Stability analysis:** the amplification factor is

$$
R(z) = \frac{1}{1-z}
$$

For $z = \lambda h$ with $\mathrm{Re}(\lambda)<0$ (stable mode): $|R(z)|<1$ for all
$\mathrm{Re}(z)<0$. Hence backward Euler is **A-stable** — unconditionally stable for
all stiff problems, at the cost of only first-order accuracy.

### 3.2 Crank–Nicolson (Implicit Trapezoid)

$$
y_{n+1} = y_n + \frac{h}{2}\bigl(f(t_n,y_n)+f(t_{n+1},y_{n+1})\bigr)
$$

**Stability:** the amplification factor is

$$
R(z) = \frac{1+z/2}{1-z/2}
$$

For $\mathrm{Re}(z)<0$: $|R(z)|<1$ (lies inside the unit circle in the left half-plane).
Hence Crank–Nicolson is **A-stable** and second-order accurate — the ODE analogue of the
trapezoidal rule used by the PDE heat solver.

---

## 4. Symplectic Integrators

### 4.1 Motivation

For Hamiltonian systems $H(q,p)=T(p)+V(q)$, standard ODE methods (RK4) suffer secular
**energy drift** — the numerical trajectory spirals inward or outward because the
symplectic structure is not preserved.

### 4.2 Leapfrog / Velocity Verlet

For separable Hamiltonians, the **kick-drift-kick** scheme:

$$
\begin{aligned}
p_{n+1/2} &= p_n - \frac{h}{2}\nabla V(q_n) \qquad\text{(half-kick)} \\
q_{n+1} &= q_n + h\,M^{-1}p_{n+1/2} \qquad\text{(drift)} \\
p_{n+1} &= p_{n+1/2} - \frac{h}{2}\nabla V(q_{n+1}) \qquad\text{(half-kick)}
\end{aligned}
$$

**Properties (proven):**
- Second-order accurate: $\mathcal{O}(h^2)$ local error.
- **Exactly symplectic:** the discrete map $\Phi_h:(q_n,p_n)\to(q_{n+1},p_{n+1})$
  preserves the symplectic 2-form: $\Phi_h^*\Omega = \Omega$.
- **Modified Hamiltonian:** the numerical method exactly conserves a modified Hamiltonian
  $\tilde{H}(q,p) = H(q,p) + \mathcal{O}(h^2)$ for exponentially long times
  $T_{\rm exp} = \mathcal{O}(e^{c/h})$.
- **No secular energy drift:** the energy error oscillates but stays bounded.

### 4.3 Yoshida's Fourth-Order Construction

Composing the second-order leapfrog map $\Phi_h$:

$$
\Phi_{\lambda h}\,\Phi_{\mu h}\,\Phi_{\lambda h}
$$

with $\lambda = (2-2^{1/3})^{-1}$, $\mu = 1-2\lambda$, cancels the leading $\mathcal{O}(h^2)$
error term, yielding a **fourth-order symplectic integrator** (Yoshida 1990). Cost: 3
force evaluations per step (vs 4 for RK4), but with symplecticity.

### 4.4 When to Use Symplectic vs RK

| Scenario | Best method |
|----------|-------------|
| Long-time Hamiltonian orbits | Symplectic (Verlet/Yoshida) |
| Smooth, non-stiff, moderate time | RK4 |
| Stiff systems | Backward Euler / Crank–Nicolson |
| Adaptive step needed | Dormand–Prince (RK45) |
| Conservation of energy critical | Symplectic |

---

## 5. Routines

| Routine | Order | Type | Best for |
|---------|-------|------|----------|
| `rk4_solve` | 4 | explicit | smooth, non-stiff, fixed step |
| `rk45_solve` | 5(4) adaptive | explicit | tolerance control |
| `backward_euler` | 1 | implicit A-stable | stiff systems |
| `crank_nicolson_ode` | 2 | implicit A-stable | moderately stiff |
| `leapfrog_step` / `velocity_verlet_step` | 2 | symplectic | Hamiltonian dynamics |
| `yoshida_step` | 4 | symplectic | high-accuracy Hamiltonian |
| `integrate_hamiltonian` | 2 | symplectic (Verlet) | separable $H$ |

---

## 6. Usage Examples

### Harmonic oscillator with RK4

```python
from pysicrs import rk4_solve

def osc(t, y):
    return [y[1], -y[0]]

times, traj = rk4_solve(osc, [1.0, 0.0], (0.0, 10.0), 10000)
print(traj[-1][0])  # ≈ cos(10) = -0.8391
```

### Long-time orbit with a symplectic integrator

```python
from pysicrs import integrate_hamiltonian

grad_v = lambda q: [q[0]]
grad_t = lambda p: [p[0]]

q, p = integrate_hamiltonian([0.0], [1.0], grad_v, grad_t, 0.01, 10000, inv_mass=1.0)

E0 = 0.5
Ef = 0.5*p[-1][0]**2 + 0.5*q[-1][0]**2
print(f"energy drift: {abs(Ef-E0):.2e}")  # tiny, bounded
```

### Stiff decay with backward Euler

```python
from pysicrs import backward_euler

def stiff(t, y):
    return [-1000.0*y[0]]

times, y = backward_euler(stiff, [1.0], (0.0, 0.1), 100)
print(y[-1][0])  # no blow-up, stable
```

---

## 7. Accuracy & Performance

| Method | Global error | Cost/step | Stability |
|--------|-------------|-----------|-----------|
| RK4 | $\mathcal{O}(h^4)$ | 4 $f$-evals | conditional |
| Dormand–Prince | adaptive to `rtol` | 7 stages | adaptive |
| Backward Euler | $\mathcal{O}(h)$ | Newton/fixed-point | A-stable |
| Crank–Nicolson | $\mathcal{O}(h^2)$ | implicit solve | A-stable |
| Verlet/leapfrog | $\mathcal{O}(h^2)$ | 1 force eval | symplectic |
| Yoshida | $\mathcal{O}(h^4)$ | 3 force evals | symplectic |

---

## 8. Advantages & Limitations

✅ Explicit methods trivial to use; symplectic integrators ideal for long Hamiltonian runs

✅ A-stable implicits handle stiff reaction kinetics

✅ Adaptivity in RK45 frees you from guessing step sizes

❌ Symplectic methods need a separable Hamiltonian

❌ RK4 has no built-in error control — choose $h$ conservatively

❌ Backward Euler is only first order; use Crank–Nicolson for accuracy

---

## 9. References

1. Hairer, E., Nørsett, S.P. & Wanner, G. (1993). *Solving Ordinary Differential Equations I*, 2nd ed. Springer.
2. Hairer, E., Lubich, C. & Wanner, G. (2006). *Geometric Numerical Integration*, 2nd ed. Springer.
3. Yoshida, H. (1990). "Construction of higher order symplectic integrators." *Phys. Lett. A* 150:262.
4. Dormand, J.R. & Prince, P.J. (1980). "A family of embedded Runge–Kutta formulae." *J. Comput. Appl. Math.* 6:19.

---

## 10. Related Topics

- [PDE](pde.md) — heat and wave solvers reuse Crank–Nicolson / FDTD ideas
- [Classical](classical.md) — supplies the $H$, $\nabla V$, $\nabla T$ for symplectic integrators
- [Calculus](calculus.md) — quadrature kernels shared with integrators
