# Quantum Mechanics

Pysic-rs implements the foundational apparatus of non-relativistic quantum mechanics —
spin-1/2 algebra, density matrices, coherent states, perturbation theory, the Feynman
path integral, and free-field propagators — all at textbook (proven) level.

> **Python binding status:** `pauli_matrices`, `density_matrix`, `coherent_state`,
> `number_state`, `feynman_propagator_scalar`, and `free_propagator_position` are bound.
> `partial_trace`, `purity`, perturbation theory, and path integrals are Rust-only. See the
> [API page](../api/quantum.md).

---

## 1. Hilbert Space & State Space

### 1.1 Axioms of Quantum Mechanics

A quantum system is described by a complex Hilbert space $\mathcal{H}$ with inner product
$\langle\cdot|\cdot\rangle$ (sesquilinear, positive-definite). States are rays — unit
vectors modulo global phase:

$$
|\psi\rangle \sim e^{i\alpha}|\psi\rangle, \qquad \langle\psi|\psi\rangle = 1.
$$

**Postulate (superposition):** if $|\psi_1\rangle$ and $|\psi_2\rangle$ are
physically realisable states, so is any normalised linear combination
$\alpha|\psi_1\rangle + \beta|\psi_2\rangle$ with $|\alpha|^2+|\beta|^2=1$.

### 1.2 Projectors & the Born Rule

For a normalised state $|\psi\rangle$, the **projector** onto the ray it spans is

$$
\hat{P}_\psi = |\psi\rangle\langle\psi|, \qquad \hat{P}_\psi^2 = \hat{P}_\psi
\quad\text{(idempotent).}
$$

The probability of outcome $a$ when measuring observable $\hat{A}$ with eigenket
$|a\rangle$ is

$$
\Pr(a) = \langle\psi|\hat{P}_a|\psi\rangle = |\langle a|\psi\rangle|^2
\qquad\text{(Born rule, proven postulate).}
$$

### 1.3 Completeness Relation

The resolution of the identity $\sum_a |a\rangle\langle a| = \hat{1}$ (for discrete,
non-degenerate spectra) is the **completeness relation**. It follows from the
orthonormality $\langle a|b\rangle = \delta_{ab}$ and the fact that the eigenkets span
$\mathcal{H}$ (self-adjoint operators have a complete orthonormal eigenbasis — spectral
theorem).

---

## 2. Spin-1/2 Algebra

### 2.1 Pauli Matrices

The Pauli matrices are the $2\times 2$ generators of $\mathfrak{su}(2)$ in the
fundamental (spin-$\frac{1}{2}$) representation:

$$
\sigma_x = \begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad
\sigma_y = \begin{pmatrix}0&-i\\i&0\end{pmatrix},\qquad
\sigma_z = \begin{pmatrix}1&0\\0&-1\end{pmatrix}.
$$

**Theorem (anticommutation).** For all $a,b\in\{1,2,3\}$:

$$
\{\sigma_a,\sigma_b\} \;=\; \sigma_a\sigma_b + \sigma_b\sigma_a \;=\; 2\delta_{ab}\,I.
$$

*Proof.* Direct matrix multiplication. For $a=b$: $\sigma_x^2 = I$ (and cyclically).
For $a\neq b$: e.g.

$$
\sigma_x\sigma_y + \sigma_y\sigma_x
= \begin{pmatrix}i&0\\0&-i\end{pmatrix}
+ \begin{pmatrix}-i&0\\0&i\end{pmatrix} = 0.
$$

$\square$

**Theorem (commutation).** The commutator is

$$
[\sigma_a,\sigma_b] = 2i\,\varepsilon_{abc}\,\sigma_c,
$$

where $\varepsilon_{abc}$ is the Levi-Civita symbol. Hence $\mathfrak{su}(2)$ has
structure constants $f^{abc} = \varepsilon^{abc}$ (see [Gauge](gauge.md)).

*Proof.* From $\sigma_a\sigma_b = \tfrac{1}{2}\{\sigma_a,\sigma_b\} +
\tfrac{1}{2}[\sigma_a,\sigma_b]$, and the anticommutator is $2\delta_{ab}I$. The
commutator is therefore purely antisymmetric and must be proportional to $\varepsilon_{abc}$
by uniqueness of the invariant tensor. Computing one component, e.g. $[\sigma_x,\sigma_y]
= 2i\sigma_z$, fixes the constant. $\square$

### 2.2 Rotation Formula

**Theorem (Euler decomposition of spin rotations).** For $\hat{n} =
(\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta)$ and
$(\hat{n}\cdot\vec{\sigma})^2 = I$ (proven by direct computation):

$$
e^{i\alpha\,\hat{n}\cdot\vec{\sigma}}
= \cos\alpha\,I + i\sin\alpha\,\hat{n}\cdot\vec{\sigma}.
$$

*Proof.* Expand the exponential as a Taylor series. Since $(\hat{n}\cdot\vec{\sigma})^{2k}
= I$ for all $k\ge 1$ and $(\hat{n}\cdot\vec{\sigma})^{2k+1} = \hat{n}\cdot\vec{\sigma}$:

$$
e^{i\alpha\hat{n}\cdot\vec{\sigma}}
= \sum_{k=0}^{\infty}\frac{(i\alpha)^k}{k!}(\hat{n}\cdot\vec{\sigma})^k
= \cos\alpha\,I + i\sin\alpha\,\hat{n}\cdot\vec{\sigma}. \qquad\square
$$

### 2.3 Spin Operators

The spin operators $S_a = \frac{\hbar}{2}\sigma_a$ satisfy the canonical
$\mathfrak{su}(2)$ Lie algebra:

$$
[S_a,S_b] = i\hbar\,\varepsilon_{abc}\,S_c, \qquad
\{S_a,S_b\} = \frac{\hbar^2}{2}\delta_{ab}\,I.
$$

The Casimir operator $\vec{S}^{\,2} = S_x^2+S_y^2+S_z^2 = \frac{3\hbar^2}{4}I$ has
eigenvalue $s(s+1)\hbar^2$ with $s=\frac{1}{2}$.

---

## 3. Density Matrix Formalism

### 3.1 Pure vs Mixed States

A **pure state** $\rho = |\psi\rangle\langle\psi|$ satisfies

$$
\rho^2 = \rho \qquad\text{(idempotent)} \qquad\Longleftrightarrow\qquad
\operatorname{tr}(\rho^2) = 1.
$$

A **mixed state** $\rho = \sum_i p_i|\psi_i\rangle\langle\psi_i|$ with $p_i>0$,
$\sum_i p_i=1$ satisfies $\operatorname{tr}(\rho^2)<1$.

**Theorem (purity bound).** For any density matrix, $0<\operatorname{tr}(\rho^2)\le 1$,
with equality on the right iff $\rho$ is pure.

*Proof.* Expand $\rho = \sum_i p_i|\psi_i\rangle\langle\psi_i|$ in an orthonormal basis
$\{|j\rangle\}$:

$$
\operatorname{tr}(\rho^2) = \sum_{i,j}p_ip_j|\langle\psi_i|\psi_j\rangle|^2
\le \sum_{i,j}p_ip_j = 1,
$$

with equality iff $|\langle\psi_i|\psi_j\rangle|^2 = 1$ for all $i,j$ with $p_ip_j>0$,
i.e. all states are identical. $\square$

### 3.2 Partial Trace & Entanglement

For a bipartite system $\mathcal{H}_{AB} = \mathcal{H}_A\otimes\mathcal{H}_B$, the
**partial trace** over $B$ is defined by

$$
\rho_A = \operatorname{tr}_B\,\rho_{AB}
= \sum_{j=1}^{d_B}\langle j|_B\;\rho_{AB}\;|j\rangle_B,
$$

where $\{|j\rangle_B\}$ is any orthonormal basis of $\mathcal{H}_B$.

**Properties (proven):**

1. Linearity: $\operatorname{tr}_B(\alpha\rho_1+\beta\rho_2)=\alpha\operatorname{tr}_B\rho_1
   +\beta\operatorname{tr}_B\rho_2$.
2. Trace preservation: $\operatorname{tr}(\rho_A)=\operatorname{tr}(\rho_{AB})$.
3. Positivity: $\rho_A\ge 0$.

**Entanglement criterion:** If $\rho_{AB}$ is pure, then $\rho_A$ is pure if and only if
$\rho_{AB}$ is separable (i.e. $\rho_{AB} = \rho_A\otimes\rho_B$). For entangled states,
$\operatorname{tr}(\rho_A^2)<1$.

### 3.3 Von Neumann Entropy

The **entropy** of a density matrix is

$$
S(\rho) = -\operatorname{tr}(\rho\ln\rho) = -\sum_i\lambda_i\ln\lambda_i,
$$

where $\{\lambda_i\}$ are the eigenvalues of $\rho$. For a pure state $S=0$; for the
maximally mixed state $\rho = I/d$, $S = \ln d$.

---

## 4. Coherent States

### 4.1 Definition & Construction

The **coherent state** $|\alpha\rangle$ of the harmonic oscillator is the eigenstate of
the annihilation operator $\hat{a}$:

$$
\hat{a}|\alpha\rangle = \alpha|\alpha\rangle, \qquad
\alpha\in\mathbb{C}.
$$

**Construction via displacement operator:**

$$
|\alpha\rangle = \hat{D}(\alpha)|0\rangle, \qquad
\hat{D}(\alpha) = \exp\!\bigl(\alpha\hat{a}^\dagger - \alpha^*\hat{a}\bigr).
$$

**Theorem (Baker-Campbell-Hausdorff).** The displacement operator satisfies

$$
\hat{D}(\alpha)^\dagger\,\hat{a}\,\hat{D}(\alpha) = \hat{a}+\alpha.
$$

*Proof.* Using $[\hat{a},\hat{a}^\dagger]=1$, the BCH formula for $e^X Y e^{-X} =
Y + [X,Y] + \tfrac{1}{2}[X,[X,Y]]+\cdots$ gives $X = \alpha\hat{a}^\dagger-\alpha^*\hat{a}$,
$[X,\hat{a}] = -\alpha$, and all higher commutators vanish. Hence
$e^X\hat{a}e^{-X} = \hat{a}-\alpha$, so $\hat{D}^\dagger\hat{a}\hat{D} = \hat{a}+\alpha$.
$\square$

**Corollary:** $\hat{a}|\alpha\rangle = \hat{a}\hat{D}(\alpha)|0\rangle =
\hat{D}(\alpha)(\hat{a}+\alpha)|0\rangle = \alpha|\alpha\rangle$.

### 4.2 Fock Expansion

Expanding in the number basis:

$$
|\alpha\rangle = e^{-|\alpha|^2/2}\sum_{n=0}^{\infty}\frac{\alpha^n}{\sqrt{n!}}|n\rangle.
$$

The probability of measuring $n$ photons is the **Poisson distribution**:

$$
P(n) = |\langle n|\alpha\rangle|^2 = \frac{|\alpha|^{2n}}{n!}e^{-|\alpha|^2}.
$$

### 4.3 Minimum Uncertainty

Coherent states saturate the Heisenberg bound:

$$
\Delta x\,\Delta p = \frac{\hbar}{2}.
$$

This follows from the fact that $|\alpha\rangle$ is a Gaussian in position space with
equal variances $\Delta x^2 = \hbar/(2m\omega)$ and $\Delta p^2 = m\omega\hbar/2$.

---

## 5. Perturbation Theory

### 5.1 Rayleigh–Schrödinger (Non-Degenerate)

**Setup.** $H = H_0 + \lambda H'$ with $H_0|n^{(0)}\rangle = E_n^{(0)}|n^{(0)}\rangle$,
non-degenerate spectrum. Expand:

$$
|n\rangle = |n^{(0)}\rangle + \lambda|n^{(1)}\rangle + \lambda^2|n^{(2)}\rangle + \cdots
$$
$$
E_n = E_n^{(0)} + \lambda E_n^{(1)} + \lambda^2 E_n^{(2)} + \cdots
$$

**Intermediate normalisation:** $\langle n^{(0)}|n^{(1)}\rangle = 0$.

**Theorem (first-order energy).**

$$
\boxed{E_n^{(1)} = \langle n^{(0)}|H'|n^{(0)}\rangle}
$$

*Proof.* Substitute the expansion into $H|n\rangle = E_n|n\rangle$ and collect
$\mathcal{O}(\lambda)$:

$$
H_0|n^{(1)}\rangle + H'|n^{(0)}\rangle = E_n^{(0)}|n^{(1)}\rangle + E_n^{(1)}|n^{(0)}\rangle.
$$

Project onto $\langle n^{(0)}|$:

$$
\underbrace{\langle n^{(0)}|H_0|n^{(1)}\rangle}_{E_n^{(0)}\langle n^{(0)}|n^{(1)}\rangle=0}
+ \langle n^{(0)}|H'|n^{(0)}\rangle
= E_n^{(1)}\underbrace{\langle n^{(0)}|n^{(0)}\rangle}_{=1}. \qquad\square
$$

**Theorem (second-order energy).**

$$
\boxed{E_n^{(2)} = \sum_{m\neq n}\frac{|\langle m^{(0)}|H'|n^{(0)}\rangle|^2}{E_n^{(0)}-E_m^{(0)}}}
$$

*Proof.* Project the $\mathcal{O}(\lambda)$ equation onto $\langle m^{(0)}|$ for
$m\neq n$:

$$
\langle m^{(0)}|H'|n^{(0)}\rangle = (E_n^{(0)}-E_m^{(0)})\langle m^{(0)}|n^{(1)}\rangle
$$

$\Longrightarrow |n^{(1)}\rangle = \sum_{m\neq n}\frac{\langle m^{(0)}|H'|n^{(0)}\rangle}
{E_n^{(0)}-E_m^{(0)}}|m^{(0)}\rangle$. Collect $\mathcal{O}(\lambda^2)$ terms and
project onto $\langle n^{(0)}|$:

$$
\langle n^{(0)}|H'|n^{(1)}\rangle = E_n^{(2)}
\quad\Longrightarrow\quad
E_n^{(2)} = \sum_{m\neq n}\frac{|\langle m^{(0)}|H'|n^{(0)}\rangle|^2}{E_n^{(0)}-E_m^{(0)}}. \qquad\square
$$

**First-order state correction:**

$$
|n^{(1)}\rangle = \sum_{m\neq n}\frac{\langle m^{(0)}|H'|n^{(0)}\rangle}{E_n^{(0)}-E_m^{(0)}}|m^{(0)}\rangle
$$

**Convergence criterion** (proven): perturbation theory converges when
$\lambda|\langle m|H'|n\rangle| \ll |E_n^{(0)}-E_m^{(0)}|$ for all $m\neq n$ — the
perturbation must be small compared with the level spacing.

### 5.2 Second-Order Stark Effect (Example)

For a hydrogen atom in an electric field $H' = eEz$, the first-order shift vanishes by
parity ($\langle 1s|z|1s\rangle=0$). The second-order shift is

$$
E_{1s}^{(2)} = e^2E^2\sum_{n\neq 1}\frac{|\langle 1s|z|n\rangle|^2}{E_1^{(0)}-E_n^{(0)}} < 0,
$$

giving the quadratic Stark effect — the atom acquires an induced dipole moment
$\mathbf{p} = -\partial E^{(2)}/\partial E$.

---

## 6. Path Integral (Feynman 1948)

### 6.1 Derivation from Canonical QM

The transition amplitude between position eigenstates is

$$
K(x_f,t_f;\,x_i,t_i) = \langle x_f|e^{-i\hat{H}(t_f-t_i)/\hbar}|x_i\rangle.
$$

**Derivation.** Insert $N$ complete sets of position eigenstates:

$$
K = \int dx_1\cdots dx_{N-1}\prod_{k=0}^{N-1}\langle x_{k+1}|e^{-i\hat{H}\Delta t/\hbar}|x_k\rangle,
$$

with $\Delta t = (t_f-t_i)/N$, $x_0=x_i$, $x_N=x_f$. Use the **Trotter product
formula** (proven for self-adjoint operators with suitable domain):

$$
e^{-i(\hat{T}+\hat{V})\Delta t/\hbar}
= e^{-i\hat{V}\Delta t/\hbar}\,e^{-i\hat{T}\Delta t/\hbar} + \mathcal{O}(\Delta t^2).
$$

The free-particle matrix element is

$$
\langle x_{k+1}|e^{-i\hat{p}^2\Delta t/(2m\hbar)}|x_k\rangle
= \sqrt{\frac{m}{2\pi i\hbar\Delta t}}\exp\!\left(\frac{im(x_{k+1}-x_k)^2}{2\hbar\Delta t}\right),
$$

so

$$
K = \lim_{N\to\infty}\int\prod_{k=1}^{N-1}dx_k\;
\left(\frac{m}{2\pi i\hbar\Delta t}\right)^{N/2}
\exp\!\left(\frac{i}{\hbar}\sum_{k=0}^{N-1}
\left[\frac{m(x_{k+1}-x_k)^2}{2\Delta t} - V(x_k)\Delta t\right]\right).
$$

In the limit $N\to\infty$, the sum in the exponent becomes the classical action
$S[x] = \int_{t_i}^{t_f}L\,dt$, giving

$$
\boxed{K(x_f,t_f;\,x_i,t_i) = \int_{x(t_i)=x_i}^{x(t_f)=x_f}\mathcal{D}x\;
e^{iS[x]/\hbar}.}
$$

### 6.2 Wick Rotation & Euclidean Path Integral

The substitution $t\to -i\tau$ (**Wick rotation**) maps the Minkowski action to the
Euclidean action $S_E = \int[\tfrac{1}{2}m(\dot{x})^2+V]\,d\tau$, and the propagator
becomes

$$
K_E(x_f,\tau_f;\,x_i,\tau_i) = \int\mathcal{D}x\;e^{-S_E[x]/\hbar}.
$$

The integrand is real and exponentially suppressed away from the classical minimum — this
is the basis of **lattice field theory**, **instanton methods**, and **ground-state
projection** (imaginary-time propagation).

### 6.3 Classical Limit

By the **method of stationary phase** (steepest descent as $\hbar\to 0$):

$$
K \sim e^{iS[x_{\rm cl}]/\hbar}\cdot\sqrt{\frac{1}{2\pi i\hbar}\frac{\partial^2 S}{\partial x_f\partial x_i}}\cdot(1+\mathcal{O}(\hbar)),
$$

where $x_{\rm cl}$ solves the classical Euler–Lagrange equations. This recovers classical
mechanics in the limit $\hbar\to 0$.

---

## 7. Feynman Propagators

### 7.1 Scalar Propagator (Momentum Space)

The Feynman propagator is the Green's function of the Klein–Gordon operator:

$$
(\partial^2+m^2)\Delta_F(x-y) = -i\delta^4(x-y).
$$

Fourier transforming:

$$
\boxed{\Delta_F(p) = \frac{i}{p^2-m^2+i\varepsilon}}
$$

The $i\varepsilon$ prescription selects Feynman boundary conditions: positive frequencies
propagate forward in time, negative frequencies backward.

### 7.2 Scalar Propagator (Position Space, Euclidean)

In Euclidean 4-space, the propagator is

$$
D(x) = \frac{m}{4\pi^2|x|}\,K_1(m|x|), \qquad D(x)\big|_{m=0} = \frac{1}{4\pi^2|x|^2},
$$

where $K_1$ is the modified Bessel function of the second kind (see
[Special Functions](special_functions.md)).

**Derivation.** The Fourier transform in 4D spherical coordinates:

$$
D(x) = \int\frac{d^4p_E}{(2\pi)^4}\frac{e^{ip_E\cdot x}}{p_E^2+m^2}
= \frac{1}{(2\pi)^4}\int_0^\infty\frac{p^3\,dp}{p^2+m^2}\int d\Omega_4\,e^{ip|x|\cos\theta}.
$$

The angular integral gives $(2\pi^2)\,J_1(p|x|)/(p|x|)$, and the radial integral is the
integral representation of $K_1$. $\square$

### 7.3 Free Fermion Propagator

For a Dirac fermion:

$$
S_F(p) = \frac{i(\not{p}+m)}{p^2-m^2+i\varepsilon}, \qquad
\not{p} = \gamma^\mu p_\mu.
$$

---

## 8. Effective Potential

### 8.1 Coleman–Weinberg Formula

The one-loop effective potential in 4D scalar field theory is

$$
V_{\rm eff}(\phi) = V_0(\phi) + \frac{1}{64\pi^2}\sum_i m_i^4(\phi)\left[\ln\frac{m_i^2(\phi)}{\mu^2} - \frac{3}{2}\right],
$$

where $m_i(\phi)$ are the field-dependent masses and $\mu$ is the renormalisation scale.

**Derivation sketch.** Start from the generating functional
$W[J] = -\frac{i}{\hbar}\ln\int\mathcal{D}\phi\,e^{iS[\phi]+i\int J\phi}$. The
effective action $\Gamma[\phi]$ is the Legendre transform of $W[J]$. The one-loop
approximation gives

$$
\Gamma[\phi] = S[\phi] + \frac{i\hbar}{2}\operatorname{tr}\ln\left(\frac{\delta^2 S}{\delta\phi^2}\right).
$$

Evaluating on a constant field and taking the functional determinant via the proper-time
method (Schwinger proper time $\int_0^\infty ds/s\,e^{-im^2s}$) yields the Coleman–Weinberg
formula. $\square$

---

## 9. Routines

| Routine | Formula | Derivation |
|---------|---------|------------|
| `pauli_matrices()` | $\sigma_x,\sigma_y,\sigma_z$ | Generators of $\mathfrak{su}(2)$ |
| `spin_operators(s)` | $S_a = \frac{\hbar}{2}\sigma_a$ | Canonical commutation relations |
| `density_matrix(psi)` | $\rho = \|\psi\rangle\langle\psi\|$ | Outer product |
| `partial_trace(rho, dims, keep)` | $\operatorname{tr}_B\rho$ | Sum over basis of $B$ |
| `purity(rho)` | $\gamma = \operatorname{tr}\rho^2$ | Entanglement criterion |
| `inner_product(cr, ci, kr, ki)` | $\langle\psi|\phi\rangle$ | Sesquilinearity |
| `projector(psi)` | $P = \|\psi\rangle\langle\psi\|$ | Idempotence |
| `coherent_state(alpha, n_max)` | $\|\alpha\rangle$ | Displacement operator |
| `number_state(n, dim)` | $\|n\rangle$ | Fock basis |
| `rayleigh_schrodinger_1st(...)` | $E_n^{(1)}, \|n^{(1)}\rangle$ | Perturbation expansion |
| `rayleigh_schrodinger_2nd(...)` | $E_n^{(2)}$ | Second-order perturbation |
| `discretized_path_integral(...)` | $\langle x_f\|x_i\rangle$ via slices | Trotter decomposition |
| `effective_potential(phi, v_tree, loop)` | $V_{\rm eff}$ grid | Coleman–Weinberg |
| `wick_rotation(x, v, m, dt)` | $S_E$ | Analytic continuation $t\to -i\tau$ |
| `feynman_propagator_scalar(p2, m, eps)` | $i/(p^2-m^2+i\varepsilon)$ | Klein–Gordon Green's function |
| `free_propagator_position(r, m)` | $\frac{m}{4\pi^2 r}K_1(mr)$ | Euclidean Fourier transform |
| `free_propagator_momentum(p0, p2, m, eps)` | $i/(p_0^2-\|\vec{p}\|^2-m^2+i\varepsilon)$ | Feynman prescription |

---

## 10. Usage Examples

### Purity detects entanglement

```python
from pysicrs import density_matrix, partial_trace, purity

# Bell state |Φ⁺⟩ = (|00⟩ + |11⟩)/√2
psi = [0.70710678, 0.0, 0.0, 0.70710678]
rho = density_matrix(psi)

rho_a = partial_trace(rho, [2, 2], [0])
print(f"purity of rho_A = {purity(rho_a):.4f}  (expect 0.5 → entangled)")
```

### Coherent state of a harmonic oscillator

```python
from pysicrs import coherent_state

alpha = coherent_state(0.5 + 0.3j, 20)
print(f"number of components: {len(alpha)}  (Poisson distribution)")
```

### Second-order Stark shift

```python
from pysicrs import rayleigh_schrodinger_2nd

# Pass H0 eigenvalues, H' matrix, index n → returns E_n^(2)
```

---

## 11. Advantages & Limitations

✅ Coherent states & Pauli algebra cover spin/optical textbook problems

✅ Path integral and propagators exported directly as amplitude-value functions

✅ Warps with the pure-Rust FFT for split-step evolution (see [PDE](pde.md))

❌ No full state-vector evolution for large $N$ (no qiskit-style gate simulator)

❌ Perturbation module assumes non-degenerate unperturbed energies

❌ Propagators are the free scalar-field forms only (no interactions/couplings)

---

## 12. References

1. Sakurai, J.J. (1994). *Modern Quantum Mechanics*, rev. ed. Addison-Wesley.
2. Shankar, R. (1994). *Principles of Quantum Mechanics*, 2nd ed. Plenum.
3. Weinberg, S. (1995). *The Quantum Theory of Fields*, Vol. 1. Cambridge.
4. Feynman, R.P. (1948). "Space-time approach to non-relativistic quantum mechanics." *Rev. Mod. Phys.* 20:367.
5. Glauber, R.J. (1963). "Coherent and incoherent states of the radiation field." *Phys. Rev.* 131:2766.
6. Coleman, S. & Weinberg, E. (1973). "Radiative corrections as the origin of spontaneous symmetry breaking." *Phys. Rev. D* 7:1888.
7. Peskin, M. & Schroeder, D. (1995). *An Introduction to Quantum Field Theory*. Addison-Wesley.

---

## 13. Related Topics

- [PDE](pde.md) — the split-step Schrödinger solver drives time evolution
- [Gauge](gauge.md) — $SU(2)$ generators are the Pauli matrices here
- [Special Functions](special_functions.md) — Legendre/spherical harmonics enter angular parts
- [Topology](topology.md) — Berry phase arises from adiabatic evolution of quantum states
