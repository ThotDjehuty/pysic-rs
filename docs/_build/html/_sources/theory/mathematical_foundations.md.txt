# Mathematical Foundations

Everything in Pysic-rs rests on **proven** results from the standard literature. This page
collects the shared mathematical skeleton — the results that bind the modules together and
make the whole library consistent. Nothing here is experimental or research-pending: these
are theorems with published proofs.

---

## 1. Complex Numbers throughout

Many modules return `Complex64` (`num_complex`). The complex plane is used by:

- **PDE** — the Schrödinger/Dirac wavefunction amplitudes
- **Quantum** — coherent states, path-integral weights, propagators
- **Fourier** — the DFT is intrinsically complex
- **Gauge** — $SU(N)$ generators and field strengths
- **Topology** — winding numbers live in $\mathbb{C}$

**Convention:** all complex returns are the *canonical* complex numbers, e.g.
$e^{i\theta} = \cos\theta + i\sin\theta$, $\sqrt{-1}=i$.

---

## 2. The Central Identity of Numerical Analysis

**Taylor's theorem with remainder** — the single most important result for everything from
finite differences to RK methods:

$$
f(x+h) = f(x) + h f'(x) + \frac{h^2}{2}f''(\xi), \qquad \xi\in(x,x+h)
$$

From it, **central differencing** gives

$$
f'(x) = \frac{f(x+h) - f(x-h)}{2h} + O(h^2)
$$

and each integration / solving scheme in this library inherits its accuracy *order* from a
matching Taylor truncation.

**Implication (proven, standard):** for double precision, the optimal finite-difference step
is $h \approx \epsilon_{\text{mach}}^{1/3} \approx 10^{-5}$ (derivatives) —
see [Calculus](../algorithms/calculus.md).

---

## 3. Variational Principles (the Unifying Idea)

Classical mechanics, quantum mechanics, electrodynamics, and field theory all rest on a
variational principle: the physical path makes an action stationary.

| Domain | Action | Stationarity gives |
|--------|--------|-------------------|
| Classical | $S = \int L\,dt$ | Euler–Lagrange equations |
| Classical (phase space) | $S = \int (p\dot q - H)dt$ | Hamilton's equations |
| Quantum | $S = \int\frac12 m\dot x^2 - V\,dt$ | Feynman path integral |
| Gauge | $S = -\frac14\int F_{\mu\nu}^aF^{a\mu\nu}$ | Yang–Mills EOM |
| GR | $S = \frac{1}{16\pi G}\int R\sqrt{-g}\,d^4x$ | Einstein field equations |

**Proven** in all five cases: stationarity of $S$ (Hamilton's principle) is equivalent to
the field/path equations.

---

## 4. Three Stability Theorems You Will Meet

1. **Lax–Richtmyer equivalence (proven, 1956):** a consistent linear finite-difference
   scheme converges **iff** it is stable (von Neumann). This justifies `CFL ≤ 1` conditions
   used by the [PDE](../algorithms/pde.md) wave/Maxwell solvers.

2. **Lax–Milgram (proven):** elliptic equations (Poisson) possess unique weak solutions —
   the basis of SOR and FFT-Poisson convergence.

3. **A-stability / L-stability (Dahlquist, proven):** for stiff ODEs, multi-step methods of
   order $>2$ cannot be A-stable; implicit single-step methods (backward Euler,
   Crank–Nicolson) can. That's why the [ODE](../algorithms/ode.md) module ships the implicit
   pair.

---

## 5. Symplectic Structure

Hamiltonian flow preserves the symplectic form $dq_i \wedge dp_i$
(**Liouville's theorem** — proven). This drives:

- The [ODE](../algorithms/ode.md) symplectic integrators
- The [Classical](../algorithms/classical.md) Poisson-bracket & canonical-transform routines

**Key proven result:** discrete symplectic integrators (Verlet, Yoshida) keep the energy
error **bounded** over arbitrarily long times, whereas generic RK schemes drift linearly in
time.

---

## 6. Index Conventions

Throughout the differential-geometry modules (GR, gauge, topology) the conventions match
Wald (1984):

- Metric signature: $(-,+,+,+)$ — time-like first
- Christoffel symbol: $\Gamma^{\rho}_{\mu\nu} = \frac12 g^{\rho\sigma}(\partial_\mu g_{\nu\sigma} + \partial_\nu g_{\sigma\mu} - \partial_\sigma g_{\mu\nu})$
- Riemann: $R^{\rho}{}_{\sigma\mu\nu} = \partial_\mu\Gamma^{\rho}_{\nu\sigma} - \partial_\nu\Gamma^{\rho}_{\mu\sigma} + \Gamma^{\rho}_{\mu\lambda}\Gamma^{\lambda}_{\nu\sigma} - \Gamma^{\rho}_{\nu\lambda}\Gamma^{\lambda}_{\mu\sigma}$
- Einstein: $G_{\mu\nu} = R_{\mu\nu} - \frac12 R g_{\mu\nu}$

Natural units $G = c = \hbar = 1$ are used in GR and quantum formulas and stated explicitly
where physical units are needed.

---

## 7. The Numbers You Can Check

| Known result | Module | Value |
|--------------|--------|-------|
| $\Gamma(n) = (n-1)!$ | Special Functions | exact |
| $\zeta(2) = \pi^2/6$ | Special Functions | $1.6449340668$ |
| Schwarzschild $g_{tt} = -(1-\frac{2M}{r})$ | GR | at $r=10M$: $-0.8$ |
| Casimir $F/A = -\frac{\pi^2\hbar c}{240d^4}$ | Casimir | at $d=1\,\mu m$: $-1.3\times10^{-3}$ N/m² |
| Harmonic $x = \cos t$ | ODE (RK4) | $\cos(10) = -0.8391$ |
| $\omega^4$ dipole law | EM | verified by scaling |
| Synmplectic energy bound | ODE (Verlet) | drift $< 10^{-4}$ over $10^4$ steps |

---

## 8. Unified Real-Data Philosophy

These are **analytic** theorems — no fitting, no Monte Carlo. As such, "mode RÉEL" in this
library means *evaluating against closed-form reference values* (the tables above). The
integration/ODE/PDE results are deterministic and reproducible.

For statistical backgrounds (Bayesian inference, MCMC, SDE simulation, Kalman filtering),
Pysic-rs delegates to **Optimiz-rs** — pairing the deterministic core with the probabilistic
machinery. See the Documentation for that library for its own mathematical foundations.

---

## Recommended Reading (proven, standard)

1. Arfken, *Mathematical Methods for Physicists*, 3rd ed. — special functions, tensors, PDEs
2. Press et al., *Numerical Recipes*, 3rd ed. — the numerical-analysis backbone
3. Hairer, Lubich & Wanner, *Geometric Numerical Integration*, 2nd ed. — symplectic theory
4. Wald, *General Relativity* — index conventions
5. Jackson, *Classical Electrodynamics*, 3rd ed. — EM results
6. Sakurai, *Modern Quantum Mechanics* — the quantum postulates

---

## Related Pages

- [Special Functions](../algorithms/special_functions.md)
- [ODE Solvers](../algorithms/ode.md)
- [General Relativity](../algorithms/general_relativity.md)
- [Quantum Mechanics](../algorithms/quantum.md)
- [Casimir Effect](../algorithms/casimir.md)