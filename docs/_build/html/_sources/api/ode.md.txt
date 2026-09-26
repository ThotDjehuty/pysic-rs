# API Reference: ODE

The Python surface exposes ODE solvers at top level (`from pysicrs import rk4_solve`).
See also [Algorithms: ODE Solvers](../algorithms/ode.md).

```python
import pysicrs
```

## Explicit methods

### `rk4_solve`

```python
rk4_solve(f, y0: list[float], t_start: float, t_end: float, n_steps: int)
    -> (times: list[float], trajectory: list[list[float]])
```

Classical 4th-order Runge–Kutta. `f(t, y)` returns `list[float]`.

```python
from pysicrs import rk4_solve

def osc(t, y):
    return [y[1], -y[0]]

times, traj = rk4_solve(osc, [1.0, 0.0], 0.0, 10.0, 10000)
print(traj[-1][0])        # cos(10) ≈ -0.8391
```

### `rk45_solve`

```python
rk45_solve(f, y0: list[float], t_start: float, t_end: float,
           rtol: float = 1e-6, atol: float = 1e-9, max_steps: int = 100_000)
    -> OdeSolution (currently returns None; adaptive solver, structured result in progress)
```

Embedded 5(4) Dormand–Prince with adaptive stepping.

## Implicit methods

### `backward_euler`

```python
backward_euler(f, y0: list[float], t_start: float, t_end: float, n_steps: int)
    -> (times: list[float], trajectory: list[list[float]])
```

A-stable implicit Euler (Picard iteration), for stiff systems.

## Symplectic integrators

### `leapfrog_integrate`

```python
leapfrog_integrate(grad_v, grad_t, q0: list[float], p0: list[float],
                   dt: float, n_steps: int, inv_mass: float)
    -> (q: list[list[float]], p: list[list[float]])
```

Kick-drift-kick Verlet leapfrog for separable $H = T(p) + V(q)$. `grad_v(q)` returns
$\nabla V$, `grad_t(p)` returns $\nabla T$; `inv_mass` is $1/m$.

```python
from pysicrs import leapfrog_integrate

q, p = leapfrog_integrate(lambda q: [q[0]], lambda p: [p[0]],
                          [0.0], [1.0], 0.01, 10000, inv_mass=1.0)
E = 0.5*p[-1][0]**2 + 0.5*q[-1][0]**2
print(f"E ≈ {E:.6f}")     # ≈ 0.5, bounded over long times
```

## Rust-only (not yet bound)

| function | purpose |
|----------|---------|
| `crank_nicolson_ode` | semi-implicit trapezoid |
| `velocity_verlet_step` | Verlet single step |
| `yoshida_step` | 4th-order symplectic |
| `integrate_hamiltonian` | path integration of separable H |