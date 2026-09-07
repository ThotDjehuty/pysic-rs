# PDE Solvers

Pysic-rs implements the canonical numerical schemes for the core PDEs of physics: the
Schrödinger and Dirac equations (quantum), the heat and wave equations (parabolic/hyperbolic
prototypes), the Poisson equation (elliptic), and the full 3D Maxwell system on a Yee grid.

> **Python binding status:** `schrodinger_split_step_1d`, `schrodinger_eigen_1d`,
> `dirac_split_step_1d`, `heat_crank_nicolson_1d/2d`, `wave_fdtd_1d/2d`, and `poisson_fft_2d`
> are bound. `poisson_sor_2d` and the Maxwell 3D routines are Rust-only. See the
> [API page](../api/pde.md).

---

## Mathematical Foundations

### Schrödinger Equation (split-step Fourier)

The time-dependent 1D Schrödinger equation

$$
i\hbar \frac{\partial\psi}{\partial t} = -\frac{\hbar^2}{2m}\frac{\partial^2\psi}{\partial x^2} + V(x)\psi
$$

is solved with the **split-step (Fourier) method**, which exploits that the kinetic and
potential operators commute formally at each *infinitesimal* step:

$$
\psi(x, t+\Delta t) \approx e^{-iV\Delta t/2\hbar}\;
e^{-i\hat{p}^2\Delta t/2\hbar m}\;
e^{-iV\Delta t/2\hbar}\;\psi(x,t) + O(\Delta t^3)
$$

**Implementation:** FFT to momentum space (diagonal kinetic phase $\exp(-i\hbar k^2\Delta t/2m)$),
inverse FFT, potential phase applied in position space. **Proven second order in $\Delta t$**
(symmetric Strang splitting), **norm-conserving** (all three factors are unitary), and the FFT
makes it **spectrally accurate** in space.

`schrodinger_eigen_1d` computes bound-state energies/states via imaginary-time propagation —
the ground state dominates as $t\to\infty$ because each bound state decays as
$e^{-E_n \tau/\hbar}$; subtracting successive energies isolates excited states.

### Dirac Equation (split-step)

The 1D Dirac equation with a scalar potential $V(x)$ (in $c=\hbar=1$ units where convenient):

$$
i\hbar \partial_t \psi = \left(-i\hbar c\,\alpha \partial_x + mc^2\beta + V\right)\psi
$$

Solved with the same split-step idea, separating the momentum-space kinetic part (which
mixes the two spinor components) from the local potential.

### Heat Equation (Crank–Nicolson + ADI)

$$
\frac{\partial u}{\partial t} = \alpha \nabla^2 u
$$

**1D:** the Crank–Nicolson scheme (trapezoidal in time, central difference in space) is the
classical **A-stable, unconditionally stable, second-order** discretization:

$$
\frac{u^{n+1}_j - u^n_j}{\Delta t} = \frac{\alpha}{2}\left[\frac{u^{n+1}_{j-1}-2u^{n+1}_j+u^{n+1}_{j+1}}{\Delta x^2} + \frac{u^{n}_{j-1}-2u^{n}_j+u^{n}_{j+1}}{\Delta x^2}\right]
$$

This is a tridiagonal linear system per step, solved with the Thomas algorithm.

**2D:** solved with **Alternating Direction Implicit (ADI)** — sweep the $x$-direction
implicitly, then the $y$-direction implicitly. **Proven second-order accurate and
unconditionally stable for constant $α$.**

### Wave Equation (FDTD)

$$
\frac{\partial^2 u}{\partial t^2} = c^2\nabla^2 u
$$

The finite-difference time-domain (FDTD) leapfrog discretization is **explicit and
second-order**, but *conditionally* stable — the **CFL condition** (proven by von Neumann
analysis) must hold:

$$
c\,\frac{\Delta t}{\Delta x} \le 1 \qquad (\text{2D: } c\Delta t \le \frac{\Delta x}{\sqrt 2})
$$

### Poisson Equation

$$
\nabla^2 \phi = -\frac{\rho}{\varepsilon_0}
$$

Two methods:

- **FFT method** — the Laplacian is diagonal in Fourier space:
  $\hat\phi(\mathbf{k}) = \hat\rho(\mathbf{k})/(\varepsilon_0 k^2)$. Used with PBCs.

- **Successive Over-Relaxation (SOR)** — Gauss–Seidel with relaxation factor
  $\omega \in (1, 2)$ (proven optimal near $\omega \approx 2/(1+\pi/n)$ for simple domains on
  $n\times n$ grids). Converges for elliptic problems.

### Maxwell Equations (3D FDTD, Yee grid)

Maxwell's curl equations in a medium (electric conductivity $\sigma$):

$$
\frac{\partial \mathbf{E}}{\partial t} = \frac{1}{\varepsilon}\nabla\times\mathbf{H} - \frac{\sigma}{\varepsilon}\mathbf{E}, \qquad
\frac{\partial \mathbf{H}}{\partial t} = -\frac{1}{\mu}\nabla\times\mathbf{E}
$$

**Yee's staggered grid** (Yee 1966) places $E$ on cell edges and $H$ on cell faces so that
every curl is a central difference and Gauss's laws are automatically satisfied to machine
precision. The **3D CFL limit**: $c\Delta t \le \Delta x/\sqrt{3}$.

---

## Routines

| Routine | Equation | Method | Stability |
|---------|----------|--------|-----------|
| `schrodinger_split_step_1d` | TDSE | split-step FFT | unconditional (unitary) |
| `schrodinger_eigen_1d` | TISE | imaginary time | — |
| `dirac_split_step_1d` | Dirac | split-step FFT | unconditional (unitary) |
| `heat_crank_nicolson_1d/2d` | heat | C-N + ADI | unconditional |
| `wave_fdtd_1d/2d` | wave | leapfrog FDTD | CFL: $c\Delta t\le\Delta x$ |
| `poisson_fft_2d` | Poisson | spectral | PBC |
| `poisson_sor_2d` | Poisson | SOR | elliptic, monotone |
| `maxwell_step_3d` | Maxwell | Yee FDTD | CFL: $c\Delta t\le\Delta x/\sqrt3$ |
| `maxwell_max_dt(dx, c)` | — | CFL helper | returns the safe $\Delta t$ |

---

## Usage Examples

### Wave packet propagation in a harmonic trap

```python
from pysicrs import schrodinger_split_step_1d
import numpy as np

n, dx = 512, 0.05
x = (np.arange(n) - n // 2) * dx
psi0 = np.exp(-(x / 0.5)**2 / 2 + 1j * 2.0 * x)  # moving Gaussian
v = 0.5 * x**2                                   # harmonic potential

psi = schrodinger_split_step_1d(
    list(psi0.real) + list(psi0.imag),  # complex as interleaved
    list(v), dx, dt=0.01, n_steps=100,
    mass=1.0, hbar=1.0,
)
```

### Heat diffusion CFL-safe

```python
from pysicrs import heat_crank_nicolson_1d

u0 = [0.0] * 40 + [1.0] * 20 + [0.0] * 40
u = heat_crank_nicolson_1d(u0, alpha=0.05, dx=0.1, dt=0.01, n_steps=200)
```

### Maxwell 3D — safe step

```python
from pysicrs import maxwell_max_dt, maxwell_step_3d

dt = maxwell_max_dt(dx=0.01, c=1.0)   # ≤ dx/√3
# state = MaxwellState(nx, ny, nz); iterate maxwell_step_3d(state, dx, dt)
```

---

## Stability Guidelines (von Neumann analysis, proven)

| Scheme | Condition |
|--------|-----------|
| Crank–Nicolson (heat) | none — unconditionally stable |
| Leapfrog FDTD (wave) | $\dfrac{c\Delta t}{\Delta x} \le 1$ |
| Maxwell 3D Yee | $\dfrac{c\Delta t}{\Delta x} \le \dfrac1{\sqrt3}$ |
| Split-step (Schrödinger/Dirac) | none — unitary |

**Practical rule:** if a simulation blows up, first check the CFL ratio; the wave and Maxwell
solvers will absolutely require it.

---

## Advantages & Limitations

✅ Split-step gives spectral spatial accuracy + norm conservation for quantum problems

✅ Crank–Nicolson and ADI remove parabolic stability constraints entirely

✅ Yee grid conserves $\nabla\cdot B=0$ discretely (proven to round-off)

❌ Wave/Maxwell FDTD is explicitly constrained by CFL — small $\Delta t$ needed

❌ FFT-Poisson requires periodic boundaries

❌ No PML absorbing boundaries implemented yet (hard-wall boundaries default)

---

## References

1. Trefethen, L.N. (2000). *Spectral Methods in MATLAB*. SIAM. (split-step FFT)
2. Yee, K.S. (1966). "Numerical solution of initial boundary value problems involving Maxwell's equations." *IEEE Trans. AP* 14(3):302–307.
3. Press, W. et al. (2007). *Numerical Recipes*, 3rd ed. Cambs. (Crank–Nicolson, SOR)
4. Strickland, J. (2011). *FDTD lecture notes* — CFL conditions, Yee lattice.
5. Saad, Y. (2003). *Iterative Methods for Sparse Linear Systems*, 2nd ed. SIAM. (SOR convergence)

---

## Related Topics

- [Fourier](fourier.md) – supplies the FFT kernels used by the spectral solvers
- [ODE](ode.md) – method-of-lines semidiscretizations reuse the RK suite
- [Casimir](casimir.md) – the mode-sum viewpoint of vacuum energy is FDTD-compatible