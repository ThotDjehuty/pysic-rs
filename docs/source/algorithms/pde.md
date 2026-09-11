# PDE Solvers

Pysic-rs implements the canonical numerical schemes for the core PDEs of physics: the
Schrödinger and Dirac equations (quantum), the heat and wave equations (parabolic/hyperbolic
prototypes), the Poisson equation (elliptic), and the full 3D Maxwell system on a Yee grid.

> **Python binding status:** `schrodinger_split_step_1d`, `schrodinger_eigen_1d`,
> `dirac_split_step_1d`, `heat_crank_nicolson_1d/2d`, `wave_fdtd_1d/2d`, and `poisson_fft_2d`
> are bound. `poisson_sor_2d` and the Maxwell 3D routines are Rust-only. See the
> [API page](../api/pde.md).

---

## 1. Classification of Linear PDEs

| Type | Prototype | Character | Solver class |
|------|-----------|-----------|--------------|
| **Elliptic** | $\nabla^2 u = f$ | No time variable | Banded/SOR/FFT |
| **Parabolic** | $\partial_t u = \alpha\nabla^2 u$ | Diffusion, infinite speed | Crank–Nicolson, ADI |
| **Hyperbolic** | $\partial_{tt}u = c^2\nabla^2 u$ | Wave propagation, finite speed | FDTD leapfrog |
| **Dispersive** | $i\hbar\partial_t\psi = -\frac{\hbar^2}{2m}\nabla^2\psi+V\psi$ | Phase dispersion | Split-step FFT |

---

## 2. Schrödinger Equation (Split-Step Fourier)

### 2.1 Time-Dependent Schrödinger

$$
i\hbar\frac{\partial\psi}{\partial t} = -\frac{\hbar^2}{2m}\frac{\partial^2\psi}{\partial x^2} + V(x)\psi
$$

### 2.2 Split-Step Method

The evolution operator over time $\Delta t$ is

$$
e^{-i\hat{H}\Delta t/\hbar} \approx e^{-iV\Delta t/2\hbar}\,e^{-i\hat{p}^2\Delta t/(2m\hbar)}\,e^{-iV\Delta t/2\hbar} + \mathcal{O}(\Delta t^3)
$$

**Derivation (Strang splitting).** Write $\hat{H} = \hat{T}+\hat{V}$ and use the Baker–
Campbell–Hausdorff formula:

$$
e^{-i(\hat{T}+\hat{V})\Delta t/\hbar} = e^{-i\hat{V}\Delta t/2\hbar}\,e^{-i\hat{T}\Delta t/\hbar}\,e^{-i\hat{V}\Delta t/2\hbar}\,e^{\frac{1}{24}[\hat{V},[\hat{V},\hat{T}]]\Delta t^3/\hbar^3+\cdots}
$$

The symmetric (Strang) splitting is second-order: $\mathcal{O}(\Delta t^3)$ local error,
$\mathcal{O}(\Delta t^2)$ global error.

**Algorithm per step:**

1. Half-step potential phase: $\psi \leftsto \psi\,e^{-iV\Delta t/(2\hbar)}$ (position space).
2. Full-step kinetic phase: $\hat{\psi} = \text{FFT}(\psi)$; $\hat{\psi} \leftsto \hat{\psi}\,
   e^{-i\hbar k^2\Delta t/(2m)}$ (momentum space).
3. Half-step potential phase: $\psi = \text{IFFT}(\hat{\psi})$; $\psi \leftsto \psi\,
   e^{-iV\Delta t/(2\hbar)}$.

**Properties (proven):**

- **Unitary:** each exponential factor has $|e^{i\theta}|=1$, so $\|\psi\|_2$ is exactly
  conserved to machine precision.
- **Spectrally accurate in space:** the spatial derivative is exact in Fourier space.
- **Second-order in time:** the Strang splitting gives $\mathcal{O}(\Delta t^2)$ global error.

### 2.3 Imaginary-Time Propagation (Eigenstates)

`schrodinger_eigen_1d` computes bound-state energies/states by evolving in **imaginary time**
$\tau = it$:

$$
\partial_\tau\psi = -\hat{H}\psi/\hbar
$$

Each eigenstate decays as $e^{-E_n\tau/\hbar}$; the ground state ($E_0$ smallest) dominates
as $\tau\to\infty$. Excited states are isolated by **Gram–Schmidt orthogonalisation** against
converged lower states.

---

## 3. Dirac Equation (Split-Step)

### 3.1 1D Dirac Equation

$$
i\hbar\partial_t\psi = \left(-i\hbar c\,\alpha\partial_x + mc^2\beta + V\right)\psi
$$

where $\alpha = \sigma_x$ and $\beta = \sigma_z$ (see [Quantum](quantum.md) for the Pauli
algebra).

### 3.2 Split-Step Decomposition

The kinetic part $-i\hbar c\,\alpha\partial_x$ is diagonal in momentum space:

$$
\hat{T}(k) = \hbar c\,\alpha\,k + mc^2\beta
$$

The two eigenvalues are $E_\pm(k) = \pm\sqrt{(\hbar ck)^2+(mc^2)^2}$ (relativistic
dispersion). The split-step proceeds as for Schrödinger, with the momentum-space phase
being a $2\times 2$ matrix exponential.

**Norm conservation:** the split-step is unitary on each spinor component, preserving
$\int|\psi_1|^2+|\psi_2|^2\,dx$.

---

## 4. Heat Equation (Crank–Nicolson)

### 4.1 1D Discretisation

$$
\frac{\partial u}{\partial t} = \alpha\frac{\partial^2 u}{\partial x^2}
$$

The **Crank–Nicolson scheme** (trapezoidal in time, central difference in space):

$$
\frac{u_j^{n+1}-u_j^n}{\Delta t} = \frac{\alpha}{2}\left[\frac{u_{j-1}^{n+1}-2u_j^{n+1}+u_{j+1}^{n+1}}{\Delta x^2} + \frac{u_{j-1}^n-2u_j^n+u_{j+1}^n}{\Delta x^2}\right]
$$

This is a **tridiagonal linear system** per step, solved with the Thomas algorithm in
$O(N)$ operations.

### 4.2 Stability Analysis (Von Neumann)

The amplification factor is

$$
g(k) = \frac{1-\sigma\sin^2(k\Delta x/2)}{1+\sigma\sin^2(k\Delta x/2)},
\qquad \sigma = \frac{\alpha\Delta t}{\Delta x^2}
$$

Since $|g(k)|\le 1$ for all $k$ and all $\sigma>0$, the scheme is **unconditionally
stable** (A-stable). It is second-order accurate in both space and time.

### 4.3 2D: ADI Method

For 2D heat, the **Alternating Direction Implicit (ADI)** method splits each step:

1. Implicit in $x$, explicit in $y$.
2. Implicit in $y$, explicit in $x$.

Each sub-step is a tridiagonal solve (along one row/column), giving $O(N^2)$ total cost
for an $N\times N$ grid. **Proven:** ADI is unconditionally stable and second-order
accurate for constant $\alpha$.

---

## 5. Wave Equation (FDTD)

### 5.1 Leapfrog Discretisation

$$
\frac{\partial^2 u}{\partial t^2} = c^2\nabla^2 u
$$

The finite-difference time-domain (FDTD) leapfrog scheme is

$$
u_j^{n+1} = 2u_j^n - u_j^{n-1} + \left(\frac{c\Delta t}{\Delta x}\right)^2(u_{j+1}^n-2u_j^n+u_{j-1}^n)
$$

**Properties:** explicit (no linear solve), second-order in space and time, conditionally
stable.

### 5.2 CFL Condition

**Theorem (Courant–Friedrichs–Lewy).** The leapfrog scheme is stable iff

$$
\boxed{c\,\frac{\Delta t}{\Delta x}\le 1}
$$

*Proof (von Neumann analysis).* Substitute $u_j^n = g^n e^{ikj\Delta x}$. The
amplification equation is $g^2-2\lambda g+1=0$ with $\lambda = 1-2\mu^2\sin^2(k\Delta x/2)$,
$\mu = c\Delta t/\Delta x$. For $|g|\le 1$ we need $|\lambda|\le 1$, giving $\mu\le 1$.
$\square$

**2D extension:** $c\Delta t \le \Delta x/\sqrt{2}$ (the Courant number must not exceed the
ratio of grid spacing to the diagonal of a cell).

---

## 6. Poisson Equation

### 6.1 FFT Method (Spectral)

$$
\nabla^2\phi = -\frac{\rho}{\varepsilon_0}
$$

In Fourier space: $\hat{\phi}(\mathbf{k}) = \hat{\rho}(\mathbf{k})/(\varepsilon_0 k^2)$.
One inverse FFT gives $\phi$. **Requires periodic boundary conditions.**

### 6.2 SOR (Successive Over-Relaxation)

Gauss–Seidel with relaxation factor $\omega\in(1,2)$:

$$
\phi_j^{(k+1)} = (1-\omega)\phi_j^{(k)} + \frac{\omega}{4}(\phi_{j-1}^{(k+1)}+\phi_{j+1}^{(k)}+\phi_j^{(k,\text{other})}+\cdots)
$$

The optimal $\omega$ for an $n\times n$ grid with Dirichlet boundaries is

$$
\omega_{\rm opt} = \frac{2}{1+\sin(\pi/n)}
$$

Convergence is monotone for $\omega\in(1,\omega_{\rm opt})$ and oscillatory for
$\omega\in(\omega_{\rm opt},2)$.

---

## 7. Maxwell 3D (Yee FDTD)

### 7.1 Maxwell's Curl Equations

$$
\varepsilon\frac{\partial\mathbf{E}}{\partial t} = \nabla\times\mathbf{H}-\sigma\mathbf{E},
\qquad
\mu\frac{\partial\mathbf{H}}{\partial t} = -\nabla\times\mathbf{E}
$$

### 7.2 Yee's Staggered Grid (Yee 1966)

The key insight: place $E$ components on **cell edges** and $H$ components on **cell faces**
in a dual lattice. Then every curl is a **central difference** and Gauss's laws
($\nabla\cdot\mathbf{B}=0$, $\nabla\cdot\mathbf{D}=\rho$) are automatically satisfied to
machine precision.

**3D CFL limit:** $c\Delta t \le \Delta x/\sqrt{3}$.

### 7.3 Energy Conservation

The Yee scheme satisfies a discrete Poynting theorem: the total electromagnetic energy
decays only due to the resistive term $\sigma\mathbf{E}\cdot\mathbf{E}$, not due to
numerical dissipation. This is a direct consequence of the central differencing and the
staggered grid.

---

## 8. Routines

| Routine | Equation | Method | Stability |
|---------|----------|--------|-----------|
| `schrodinger_split_step_1d` | TDSE | split-step FFT | unconditional (unitary) |
| `schrodinger_eigen_1d` | TISE | imaginary time | — |
| `dirac_split_step_1d` | Dirac | split-step FFT | unconditional (unitary) |
| `heat_crank_nicolson_1d/2d` | heat | C-N + ADI | unconditional |
| `wave_fdtd_1d/2d` | wave | leapfrog FDTD | CFL: $c\Delta t\le\Delta x$ |
| `poisson_fft_2d` | Poisson | spectral | PBC |
| `poisson_sor_2d` | Poisson | SOR | elliptic, monotone |
| `maxwell_step_3d` | Maxwell | Yee FDTD | CFL: $c\Delta t\le\Delta x/\sqrt{3}$ |
| `maxwell_max_dt(dx, c)` | — | CFL helper | returns the safe $\Delta t$ |

---

## 9. Usage Examples

### Wave packet propagation in a harmonic trap

```python
from pysicrs import schrodinger_split_step_1d
import numpy as np

n, dx = 512, 0.05
x = (np.arange(n) - n//2) * dx
psi0 = np.exp(-(x/0.5)**2/2 + 1j*2.0*x)  # moving Gaussian
v = 0.5 * x**2                              # harmonic potential

psi = schrodinger_split_step_1d(
    list(psi0.real) + list(psi0.imag),
    list(v), dx, dt=0.01, n_steps=100,
    mass=1.0, hbar=1.0,
)
```

### Heat diffusion CFL-safe

```python
from pysicrs import heat_crank_nicolson_1d

u0 = [0.0]*40 + [1.0]*20 + [0.0]*40
u = heat_crank_nicolson_1d(u0, alpha=0.05, dx=0.1, dt=0.01, n_steps=200)
```

### Maxwell 3D — safe step

```python
from pysicrs import maxwell_max_dt, maxwell_step_3d

dt = maxwell_max_dt(dx=0.01, c=1.0)  # ≤ dx/√3
```

---

## 10. Stability Guidelines (von Neumann analysis, proven)

| Scheme | Condition |
|--------|-----------|
| Crank–Nicolson (heat) | none — unconditionally stable |
| Leapfrog FDTD (wave) | $c\Delta t/\Delta x \le 1$ |
| Maxwell 3D Yee | $c\Delta t/\Delta x \le 1/\sqrt{3}$ |
| Split-step (Schrödinger/Dirac) | none — unitary |

**Practical rule:** if a simulation blows up, first check the CFL ratio.

---

## 11. Advantages & Limitations

✅ Split-step gives spectral spatial accuracy + norm conservation for quantum problems

✅ Crank–Nicolson and ADI remove parabolic stability constraints entirely

✅ Yee grid conserves $\nabla\cdot B=0$ discretely (proven to round-off)

❌ Wave/Maxwell FDTD is explicitly constrained by CFL — small $\Delta t$ needed

❌ FFT-Poisson requires periodic boundaries

❌ No PML absorbing boundaries implemented yet (hard-wall boundaries default)

---

## 12. References

1. Trefethen, L.N. (2000). *Spectral Methods in MATLAB*. SIAM.
2. Yee, K.S. (1966). "Numerical solution of initial boundary value problems involving Maxwell's equations." *IEEE Trans. AP* 14(3):302–307.
3. Press, W. et al. (2007). *Numerical Recipes*, 3rd ed. Cambs.
4. Strikwerda, J.C. (2004). *Finite Difference Schemes and Partial Differential Equations*, 2nd ed. SIAM.
5. Taflove, A. & Hagness, S.C. (2005). *Computational Electrodynamics*, 3rd ed. Artech House.

---

## 13. Related Topics

- [Fourier](fourier.md) — supplies the FFT kernels used by the spectral solvers
- [ODE](ode.md) — method-of-lines semidiscretisations reuse the RK suite
- [Casimir](casimir.md) — the mode-sum viewpoint of vacuum energy is FDTD-compatible
- [Quantum](quantum.md) — the split-step Schrödinger solver is the main time-evolution engine
