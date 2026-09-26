# From Quantisation to Quantum Gravity

This page gives the mathematical formulation behind six equations that pysic-rs
solves, in the order in which they were historically forced on physics:
Schrödinger, Dirac, Maxwell/Aharonov–Bohm, the ADM constraints, and
Wheeler–DeWitt. Each section states the equation, derives the scheme pysic-rs
uses to solve it, and names the invariant that the implementation is graded
against.

```{note}
Every numerical figure quoted here is measured, not asserted. The companion
notebook
[`14_dirac_adm_wheeler_dewitt.ipynb`](https://github.com/ThotDjehuty/pysic-rs/blob/main/notebooks/14_dirac_adm_wheeler_dewitt.ipynb)
reproduces all of them, grading each solver against a closed form where one
exists and against SciPy's Airy function where one does not.
```

---

## 1. Schrödinger, and the flaw that forced Dirac

### 1.1 The equation

$$
i\hbar\,\frac{\partial\psi}{\partial t}
 = \hat H\psi
 = \left[-\frac{\hbar^2}{2m}\nabla^2 + V(\mathbf r)\right]\psi
$$

The equation is **first order in $t$** and **second order in $\mathbf x$**. Under a
Lorentz boost those derivatives transform differently, so the equation cannot be
covariant. That asymmetry is not a technical nuisance to be patched — it is what
sent Dirac looking for a first-order-in-both alternative in 1928.

### 1.2 Splitting the propagator

The formal solution is $\psi(t) = e^{-i\hat Ht/\hbar}\psi(0)$. Writing
$\hat H = \hat T + \hat V$, the two pieces are each diagonal in a different
basis: $\hat T$ in momentum space, $\hat V$ in position space. They do not
commute, so $e^{-i(\hat T+\hat V)\Delta t/\hbar} \neq e^{-i\hat T\Delta t/\hbar}e^{-i\hat V\Delta t/\hbar}$.

**Theorem (Strang splitting).** For bounded operators,

$$
e^{-i(\hat T+\hat V)\Delta t/\hbar}
 = e^{-i\hat V\Delta t/2\hbar}\,
   e^{-i\hat T\Delta t/\hbar}\,
   e^{-i\hat V\Delta t/2\hbar}
 + \mathcal O(\Delta t^3).
$$

*Proof sketch.* Expand both sides through third order in $\Delta t$ using
Baker–Campbell–Hausdorff. The first-order terms match trivially. The
second-order terms match because the symmetric arrangement cancels the
commutator $[\hat T,\hat V]$ that spoils the naive product. The leading
discrepancy is the double commutator at $\mathcal O(\Delta t^3)$ per step,
giving $\mathcal O(\Delta t^2)$ global accuracy. $\square$

Each factor is the exponential of a Hermitian operator, hence **unitary**. A
product of unitaries is unitary, so

$$
\mathcal N = \int|\psi|^2\,dx
$$

is an exact invariant of the discrete scheme — independent of $\Delta t$, and
not merely conserved to truncation order. This is the property to test: either
the norm holds to round-off, or the implementation is wrong.

### 1.3 What it is graded against

For $V=0$, a Gaussian stays Gaussian:

$$
\sigma(t) = \sigma_0\sqrt{1+\left(\frac{\hbar t}{m\sigma_0^2}\right)^{2}},
\qquad
\langle x\rangle(t) = \langle x\rangle_0 + \frac{\hbar k_0}{m}t.
$$

{func}`schrodinger_split_step_1d` reproduces this envelope to **1.6 × 10⁻¹⁴**,
with norm drift **2 × 10⁻¹⁴**.

### 1.4 Bound states

The time-independent problem is discretised on a uniform grid with Dirichlet
boundaries, giving a symmetric tridiagonal Hamiltonian:

$$
H_{ii} = \frac{\hbar^2}{m\,\Delta x^2} + V_i,
\qquad
H_{i,i\pm1} = -\frac{\hbar^2}{2m\,\Delta x^2}.
$$

For a flat potential this has an **exact** spectrum,

$$
E_k = \frac{2\hbar^2}{m\,\Delta x^2}\sin^2\!\left(\frac{k\pi}{2(N+1)}\right),
$$

which converges to $k^2\pi^2\hbar^2/2mL^2$ with $L=(N+1)\Delta x$ as
$\Delta x\to0$. {func}`schrodinger_eigen_1d` uses the implicit QL algorithm with
Wilkinson shifts and matches this to machine precision; on the harmonic
oscillator it returns $E_n=(n+\tfrac12)\hbar\omega$ to $3\times10^{-6}$, with
eigenvectors overlapping the exact Hermite functions to $10^{-9}$.

```{warning}
Count nodes above an amplitude threshold. In the classically forbidden tails
$|\psi|\sim10^{-300}$ and its sign is pure round-off, so a naive sign-change
count returns nonsense — 90 nodes for the ground state rather than 0.
```

---

## 2. Dirac: relativity forces antimatter

### 2.1 The Clifford algebra

Dirac wanted first order in both $t$ and $\mathbf x$:

$$
i\hbar\,\partial_t\psi = \left(c\,\boldsymbol\alpha\cdot\mathbf p + \beta mc^2\right)\psi .
$$

Iterating must return the relativistic dispersion $E^2 = p^2c^2+m^2c^4$.
Squaring the right-hand side gives cross terms $\alpha_i\alpha_j + \alpha_j\alpha_i$
and $\alpha_i\beta+\beta\alpha_i$, which must vanish, while the squares must be
the identity:

$$
\boxed{\{\gamma^\mu,\gamma^\nu\} = 2\eta^{\mu\nu}\mathbb 1},
\qquad
\alpha_i^2=\beta^2=\mathbb 1,
\qquad
\{\alpha_i,\beta\}=0 .
$$

Numbers commute, so no set of numbers can satisfy this. The coefficients must be
**matrices**, and $\psi$ a multi-component **spinor**. Two things follow without
being assumed:

1. **Spin ½** appears in the structure of the equation, not as a postulate.
2. The spectrum has **two branches**, $E_\pm = \pm\sqrt{p^2c^2+m^2c^4}$, separated
   by a gap $2mc^2$. Dirac read the lower branch as a new particle in 1931;
   Anderson photographed the positron in 1932.

In 1+1 dimensions the Pauli matrices supplied by {func}`pauli_matrices` realise
the algebra with $\alpha=\sigma_z$, $\beta=\sigma_x$, giving

$$
\hat H = c\,\sigma_z\,\hat p + mc^2\sigma_x + V .
$$

### 2.2 The propagator, and why unitarity is the test

pysic-rs splits the evolution as $V/2\cdot M/2\cdot K\cdot M/2\cdot V/2$. The
mass factor is a rotation in the $(\psi_L,\psi_R)$ plane:

$$
e^{-i\sigma_x\theta/2}
 = \cos\frac{\theta}{2}\,\mathbb 1 - i\sin\frac{\theta}{2}\,\sigma_x,
\qquad \theta = \frac{mc^2\Delta t}{\hbar}.
$$

```{admonition} A bug worth recording
:class: warning

An earlier version of this rotation computed `exp(i·sin(θ/2))` — a complex
exponential of the sine, of modulus 1 — where `i·sin(θ/2)` was meant. The step
matrix became approximately $\begin{pmatrix}1&-1\\-1&1\end{pmatrix}$, whose
singular values are 2 and 0. The norm then grew by a fixed factor every step
**regardless of $\Delta t$**: a Gaussian went from 1.77 to $3.6\times10^{57}$ in
100 steps. Because the growth is $\Delta t$-independent, shrinking the step size
does not help — which is exactly the signature of a broken operator rather than
an unstable one.
```

{func}`dirac_split_step_1d` now holds the norm to **1.4 × 10⁻¹²** over 5000
steps, drifting linearly with step count as accumulated round-off should.

```{important}
The spinor is passed **interleaved**: element $2i$ is $\psi_L(x_i)$ and $2i+1$
is $\psi_R(x_i)$, so the array has length $2N$. Packing it as two concatenated
blocks silently scrambles the two components and produces plausible-looking but
meaningless output.
```

### 2.3 Zitterbewegung: antimatter you can measure

A localised state is not an energy eigenstate — it superposes both branches, and
the interference between them is observable. At $p=0$ the Hamiltonian reduces to
$mc^2\sigma_x$, and the evolution of a state starting in the upper component is
exactly solvable:

$$
\psi(t) = \begin{pmatrix}\cos(mc^2t/\hbar)\\ -i\sin(mc^2t/\hbar)\end{pmatrix}
\quad\Longrightarrow\quad
|\psi_R(t)|^2 = \frac{1-\cos(2mc^2t/\hbar)}{2}.
$$

The lower-component population therefore beats at the **Compton frequency**
$\omega_{\rm zb} = 2mc^2/\hbar \approx 1.6\times10^{21}\,\mathrm{s^{-1}}$ for an
electron. For a packet of finite width each momentum component beats at
$2E_p/\hbar$, so a **broad packet in $x$** (narrow in $p$) is needed to keep them
in phase.

```{tip}
Measuring a frequency requires integrating over several periods. A first version
of this experiment ran to $t=1.2$ when the period is $\pi$ — less than half a
cycle — and its FFT returned a frequency wrong by 160% with no error raised. A
correct solver is no defence against a wrong measurement protocol.
```

---

## 3. Maxwell and Aharonov–Bohm: the potential is physical

### 3.1 Gauge freedom

Maxwell's equations are unchanged under

$$
\mathbf A \to \mathbf A + \nabla\chi,
\qquad
\phi \to \phi - \partial_t\chi ,
$$

which for a century was read as evidence that $(\phi,\mathbf A)$ are mere
bookkeeping and only $\mathbf E,\mathbf B$ are real.

### 3.2 The counterexample

Aharonov and Bohm showed in 1959 that an electron circling a solenoid picks up

$$
\boxed{\Delta\varphi_{\rm AB}
 = \frac{q}{\hbar}\oint\mathbf A\cdot d\boldsymbol\ell
 = \frac{q\Phi}{\hbar}}
$$

even though $\mathbf B = 0$ **everywhere along its path**. The line integral is
gauge invariant because $\oint\nabla\chi\cdot d\boldsymbol\ell = 0$ for
single-valued $\chi$. The potential carries information the local fields do not;
electromagnetism is a **gauge theory**, and the observable is a phase.

### 3.3 Geometric phase

This is the first known case of what Berry generalised in 1984. Transporting a
state around a closed loop in parameter space leaves

$$
\gamma = i\oint\langle\psi|\nabla_{\!R}\psi\rangle\cdot d\mathbf R ,
$$

which depends only on the geometry of the loop, not the rate of traversal. For a
spin-½ on a cone of solid angle $\Omega$, $\gamma = -\Omega/2$; an equatorial
loop subtends $\Omega = 2\pi$, so $\gamma = -\pi$ and the state returns with its
sign reversed. {func}`berry_phase` reproduces $-\pi$ to $10^{-10}$.

The associated invariant is an **integer**: {func}`winding_number` returns
$n=\pm1,\pm2,\pm3$ exactly for $e^{in\theta}$. Quantisation is what makes the AB
effect robust to disorder, and it reappears as the Chern number in the quantum
Hall effect — via the Harper equation that draws the Hofstadter butterfly in the
project logo.

```{note}
{func}`green_fn_static` returns $-1/4\pi r$, i.e. the solution of
$\nabla^2G = +\delta$. Flip the sign for the usual Coulomb convention
$\nabla^2G = -\delta$.
```

---

## 4. ADM: general relativity as a constrained system

### 4.1 The 3+1 split

Arnowitt, Deser and Misner sliced spacetime into spatial hypersurfaces labelled
by a time function. The line element becomes

$$
ds^2 = -N^2dt^2 + \gamma_{ij}\left(dx^i + N^idt\right)\left(dx^j + N^jdt\right),
$$

with $\gamma_{ij}$ the induced spatial metric, $N$ the **lapse** and $N^i$ the
**shift**. The extrinsic curvature measures how the slice bends inside
spacetime:

$$
K_{ij} = \frac{1}{2N}\left(\partial_t\gamma_{ij} - D_iN_j - D_jN_i\right).
$$

### 4.2 The constraints

Projecting Einstein's equations normal and tangential to the slice splits them
into six evolution equations and four equations containing **no time
derivative**:

$$
\boxed{\mathcal H = {}^{(3)}\!R + K^2 - K_{ij}K^{ij} - 16\pi\rho = 0}
$$

$$
\boxed{\mathcal M^i = D_j\left(K^{ij} - \gamma^{ij}K\right) - 8\pi j^i = 0}
$$

These are $\mathcal H = 2G_{\mu\nu}n^\mu n^\nu - 16\pi\rho$ and
$\mathcal M^i = -2G_{\mu\nu}n^\mu\gamma^{\nu i} - 8\pi j^i$. The Bianchi
identities guarantee they are preserved by the evolution, so it is enough to
impose them on the initial data.

The consequence is structural: **the total Hamiltonian of general relativity is
a constraint**, $\mathcal H \approx 0$. It does not generate evolution in an
external time, because there is no external time. This is what produces the
problem of time on quantisation.

### 4.3 Closed forms to grade against

On a flat slice (${}^{(3)}\!R=0$):

| Configuration | $\mathcal H$ |
|---|---|
| $K_{ij}=0$, $\rho=0$ | $0$ |
| $K_{ij}=0$, $\rho\neq0$ | $-16\pi\rho$ |
| $K_{ij}=c\,\delta_{ij}$ | $K^2-K_{ij}K^{ij}=9c^2-3c^2=6c^2$ |

{func}`hamiltonian_constraint` and {func}`momentum_constraint` reproduce all
three with **exactly zero** deviation.

```{admonition} Schwarzschild must be Ricci-flat
:class: warning

Schwarzschild is a vacuum solution, so $R \equiv 0$ at every radius and every
polar angle. An earlier {func}`ricci_scalar` returned a stable $-2/r^2$ instead,
because {func}`schwarzschild_metric` freezes $g_{\varphi\varphi}=r^2$ "assuming
$\theta=\pi/2$". The curvature routines *differentiate* the metric, so the
missing $\sin^2\theta$ meant $\partial_\theta g_{\varphi\varphi}=0$ and the
angular terms stopped cancelling. Use {func}`schwarzschild_metric_full`
whenever the metric will be differentiated; $R$ is then $\sim10^{-8}$, the
finite-difference residual.

The test that should have caught this passed a *constant* metric closure to the
routine — every derivative was zero, so $R$ vanished trivially. A test that
cannot fail is worse than no test, because it inspires confidence.
```

---

## 5. Wheeler–DeWitt: quantising the constraint

### 5.1 The equation, and the disappearance of time

Canonical quantisation promotes the constraint to an operator equation on a wave
functional of the spatial geometry:

$$
\boxed{\hat{\mathcal H}\,\Psi[\gamma_{ij}] = 0}
$$

There is **no time derivative**. The state of the universe does not evolve; it
is. This is the **problem of time**, and it remains open.

### 5.2 Minisuperspace

Freezing every degree of freedom but the scale factor $a$ reduces the functional
equation to an ODE:

$$
\frac{d^2\Psi}{da^2} + \frac{p}{a}\frac{d\Psi}{da} - U(a)\Psi = 0,
\qquad
U(a) = a^2 - \frac{\Lambda}{3}a^4 ,
$$

with $p$ an operator-ordering parameter. The potential vanishes at the **turning
point**

$$
a_t = \sqrt{3/\Lambda},
$$

which separates the classically forbidden region from the expanding universe.
Tunnelling through it is the universe-creation mechanism of Vilenkin and of
Hartle–Hawking.

### 5.3 The Airy connection

Near $a_t$,

$$
U'(a_t) = 2a_t - \frac{4\Lambda}{3}a_t^3 = -2\sqrt{3/\Lambda},
\qquad
U(a) \simeq U'(a_t)\,(a-a_t).
$$

Substituting $z = -\left[2\sqrt{3/\Lambda}\right]^{1/3}(a-a_t)$ turns the
equation into

$$
\frac{d^2\Psi}{dz^2} = z\,\Psi ,
$$

the **Airy equation**, whose decaying solution is $\mathrm{Ai}(z)$. This makes
the problem gradable: integrate the full ODE with {func}`rk45_solve` and compare
against $\mathrm{Ai}$ in the window where the linearisation holds. The agreement
is $9.2\times10^{-3}$, degrading away from $a_t$ — which is the linearisation
failing, not the integrator.

```{note}
{func}`airy_ai` is also available directly, machine-precision for $|x|\le9$ via
the Maclaurin series $\mathrm{Ai} = c_1f - c_2g$ and $1.2\times10^{-7}$ beyond
it from the four-term asymptotic expansion.
```

---

## 6. Summary

| Equation | Year | What it adds | Routine | Graded against |
|---|---|---|---|---|
| Schrödinger | 1926 | quantisation | {func}`schrodinger_split_step_1d`, {func}`schrodinger_eigen_1d` | exact Gaussian, $1.6\times10^{-14}$ |
| Dirac | 1928 | relativity, spin, antimatter | {func}`dirac_split_step_1d`, {func}`pauli_matrices` | unitarity, $1.4\times10^{-12}$ |
| Maxwell / AB | 1865 / 1959 | gauge, geometric phase | {func}`green_fn_static`, {func}`berry_phase`, {func}`winding_number` | $\gamma=-\pi$, integer winding |
| ADM | 1959 | GR as constrained evolution | {func}`hamiltonian_constraint`, {func}`ricci_scalar` | $-16\pi\rho$, $6c^2$, exact |
| Wheeler–DeWitt | 1967 | quantised constraint | {func}`rk45_solve`, {func}`airy_ai` | Airy function |

The progression is one idea pursued to its end. Schrödinger quantises a
particle. Dirac makes that relativistic and inherits antimatter.
Maxwell–Aharonov–Bohm show the phase, not the field, is fundamental. ADM
rewrites gravity as a constrained system. Wheeler–DeWitt quantises that
constraint — and finds that time has dropped out of the equation.

## Related pages

- {doc}`mathematical_foundations` — the numerical principles underneath all of this
- {doc}`../algorithms/pde` — Schrödinger and Dirac solvers
- {doc}`../algorithms/general_relativity` — metrics, curvature, ADM
- {doc}`../algorithms/topology` — Berry phase, Chern and winding numbers
- {doc}`../algorithms/special_functions` — Airy, Bessel and the rest, now plotted
