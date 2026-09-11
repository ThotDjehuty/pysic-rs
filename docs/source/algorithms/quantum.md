# Quantum Mechanics

Pysic-rs covers non-relativistic quantum mechanics from the founding postulates — spin
algebra, density matrices, perturbation theory, path integrals, and free-field propagators —
using only textbook (proven) results.

> **Python binding status:** `pauli_matrices`, `density_matrix`, `coherent_state`,
> `number_state`, `feynman_propagator_scalar`, and `free_propagator_position` are bound.
> `partial_trace`, `purity`, perturbation theory, and path integrals are Rust-only. See the
> [API page](../api/quantum.md).

---

## Mathematical Foundations

### Pauli Matrices & Spin Operators

The Pauli matrices are the generators of $\mathfrak{su}(2)$ in the fundamental representation:

$$
\sigma_x = \begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
\sigma_y = \begin{pmatrix}0&-i\\i&0\end{pmatrix},\quad
\sigma_z = \begin{pmatrix}1&0\\0&-1\end{pmatrix}
$$

**Derivation of algebraic properties:**

The fundamental anticommutation relation $\{\sigma_a, \sigma_b\} = 2\delta_{ab}I$ follows from direct matrix multiplication. For $a \neq b$, the off-diagonal entries cancel:

$$
\sigma_x\sigma_y + \sigma_y\sigma_x = \begin{pmatrix}i&0\\0&-i\end{pmatrix} + \begin{pmatrix}-i&0\\0&i\end{pmatrix} = 0
$$

The commutation relation $[\sigma_a, \sigma_b] = 2i\varepsilon_{abc}\sigma_c$ is derived using the Jacobi identity and the fact that $\mathfrak{su}(2)$ is a simple Lie algebra with structure constants $\varepsilon_{abc}$.

**Rotation formula derivation:**

For $\hat{n} = (\sin\theta\cos\phi, \sin\theta\sin\phi, \cos\theta)$, we have $(\hat{n}\cdot\vec{\sigma})^2 = I$ (proven by direct computation). Therefore:

$$
e^{i\theta\hat{n}\cdot\vec{\sigma}} = \sum_{k=0}^\infty \frac{(i\theta)^k}{k!}(\hat{n}\cdot\vec{\sigma})^k = \cos\theta\,I + i\sin\theta\,\hat{n}\cdot\vec{\sigma}
$$

where we used the Taylor expansion and the idempotence property.

**Spin operators** are the doubled matrices $S_a = \frac{\hbar}{2}\sigma_a$; they obey the
canonical $\mathfrak{su}(2)$ commutation relations $[S_a,S_b]=i\hbar\varepsilon_{abc}S_c$.

### Density Matrices & Partial Trace

**Derivation of purity criterion:**

For a pure state $\rho = |\psi\rangle\langle\psi|$, we have $\rho^2 = \rho$ (idempotence), so $\operatorname{tr}(\rho^2) = \operatorname{tr}(\rho) = 1$. For a mixed state $\rho = \sum_i p_i |\psi_i\rangle\langle\psi_i|$ with $p_i > 0$, $\sum_i p_i = 1$, we compute:

$$
\operatorname{tr}(\rho^2) = \sum_{i,j} p_i p_j |\langle\psi_i|\psi_j\rangle|^2 < \sum_{i,j} p_i p_j = 1
$$

unless the states are orthogonal and the mixture is trivial.

**Partial trace derivation:**

The partial trace over system B is defined as:

$$
\rho_A = \operatorname{tr}_B\,\rho_{AB} = \sum_j \langle j|_B\ \rho_{AB}\ |j\rangle_B
$$

This operation is linear, preserves trace ($\operatorname{tr}(\rho_A) = \operatorname{tr}(\rho_{AB})$), and yields a positive semidefinite operator. The entanglement criterion follows: if $\rho_{AB}$ is a pure state, then $\rho_A$ is pure iff $\rho_{AB}$ is separable.

### Hilbert Space

- **Inner product** $\langle\psi|\phi\rangle = \sum_i \psi_i^*\,\phi_i$ (sesquilinear, positive definite).
- **Projector** $P = |\psi\rangle\langle\psi|$ (idempotent: $P^2 = P$).
- **Coherent state** of the harmonic oscillator (eigenstate of the annihilation operator,
  minimum-uncertainty, Glauber 1963):

$$
|\alpha\rangle = e^{-|\alpha|^2/2}\sum_{n=0}^\infty \frac{\alpha^n}{\sqrt{n!}}|n\rangle, \qquad
\hat a|\alpha\rangle = \alpha|\alpha\rangle, \qquad
|\alpha\rangle = \hat D(\alpha)|0\rangle,\ \hat D(\alpha)=e^{\alpha\hat a^\dagger - \alpha^*\hat a}
$$

**Derivation of coherent state properties:**

The displacement operator $\hat{D}(\alpha)$ satisfies $\hat{D}(\alpha)^\dagger \hat{a} \hat{D}(\alpha) = \hat{a} + \alpha$ (Baker-Campbell-Hausdorff formula). Applying to vacuum:

$$
\hat{a}|\alpha\rangle = \hat{a}\hat{D}(\alpha)|0\rangle = \hat{D}(\alpha)(\hat{a} + \alpha)|0\rangle = \alpha|\alpha\rangle
$$

The Poisson distribution follows from $|\langle n|\alpha\rangle|^2 = e^{-|\alpha|^2}|\alpha|^{2n}/n!$.

- **Number state** $|n\rangle$, $\hat N|n\rangle = n|n\rangle$.

### Rayleigh–Schrödinger Perturbation Theory

**Derivation of first-order energy shift:**

For $H = H_0 + \lambda H'$ with non-degenerate unperturbed spectrum $E_n^{(0)}$, $|n^{(0)}\rangle$, we expand $|n\rangle = |n^{(0)}\rangle + \lambda|n^{(1)}\rangle + \lambda^2|n^{(2)}\rangle + \cdots$ and $E_n = E_n^{(0)} + \lambda E_n^{(1)} + \lambda^2 E_n^{(2)} + \cdots$.

Substituting into $H|n\rangle = E_n|n\rangle$ and collecting $\lambda^1$ terms:

$$
H_0|n^{(1)}\rangle + H'|n^{(0)}\rangle = E_n^{(0)}|n^{(1)}\rangle + E_n^{(1)}|n^{(0)}\rangle
$$

Projecting onto $\langle n^{(0)}|$ and using $\langle n^{(0)}|n^{(1)}\rangle = 0$ (intermediate normalization):

$$
E_n^{(1)} = \langle n^{(0)}|H'|n^{(0)}\rangle
$$

**Second-order energy shift:**

Projecting the $\lambda^1$ equation onto $\langle m^{(0)}|$ for $m \neq n$:

$$
\langle m^{(0)}|H'|n^{(0)}\rangle = (E_n^{(0)} - E_m^{(0)})\langle m^{(0)}|n^{(1)}\rangle
$$

Thus $|n^{(1)}\rangle = \sum_{m\neq n} \frac{\langle m^{(0)}|H'|n^{(0)}\rangle}{E_n^{(0)} - E_m^{(0)}}|m^{(0)}\rangle$. Substituting back and projecting onto $\langle n^{(0)}|$ gives:

$$
E_n^{(2)} = \sum_{m\neq n}\frac{\left|\langle m^{(0)}|H'|n^{(0)}\rangle\right|^2}{E_n^{(0)} - E_m^{(0)}}
$$

**Proven convergence** requires the perturbation to be small compared with the level spacing: $\lambda|\langle m|H'|n\rangle| \ll |E_n^{(0)} - E_m^{(0)}|$ for all $m \neq n$.

### Path Integral (Feynman 1948)

**Derivation of the path integral representation:**

The transition amplitude is:

$$
\langle x_f|e^{-i\hat H T/\hbar}|x_i\rangle = \int \mathcal D x\; e^{\frac{i}{\hbar}\int_0^T \left[\tfrac12 m\dot x^2 - V(x)\right]dt}
$$

**Derivation:** We insert $N$ complete sets of position eigenstates:

$$
\langle x_f|e^{-i\hat H T/\hbar}|x_i\rangle = \int dx_1 \cdots dx_{N-1} \prod_{k=0}^{N-1} \langle x_{k+1}|e^{-i\hat H \Delta t/\hbar}|x_k\rangle
$$

Using the Trotter formula $e^{-i(\hat T + \hat V)\Delta t/\hbar} \approx e^{-i\hat V\Delta t/\hbar}e^{-i\hat T\Delta t/\hbar}$ and the free particle kernel:

$$
\langle x_{k+1}|e^{-i\hat p^2\Delta t/(2m\hbar)}|x_k\rangle = \sqrt{\frac{m}{2\pi i\hbar\Delta t}} \exp\left(\frac{im(x_{k+1} - x_k)^2}{2\hbar\Delta t}\right)
$$

**Discretization:** with $N$ time slices $\Delta t = T/N$ the free-particle prefactor is
$\left(\frac{m}{2\pi i\hbar\Delta t}\right)^{1/2}$ — the product of free propagators times
potential weights $\exp(-iV\Delta t/\hbar)$ converges (proven) to the amplitude as $N\to\infty$.

**Wick rotation** $t\to -i\tau$ maps this to the Euclidean (imaginary-time) path integral with
$S_E = \int\left[\tfrac12 m(\tfrac{dx}{d\tau})^2 + V\right]d\tau$, making the integrand real
and sharply peaked at the classical minimum — the basis of instanton and ground-state methods.

### Effective Potential

**Derivation of Coleman-Weinberg formula:**

The one-loop effective potential (Coleman–Weinberg) is the Legendre transform of the free
energy $\Gamma[\phi] = -W[J] + \int J\phi$, summed:

$$
V_{\text{eff}}(\phi) = V_0(\phi) + \frac{\hbar^2}{64\pi^2}\left[\ldots\right]\quad\text{(1-loop)}
$$

The derivation starts from the generating functional $W[J] = -\frac{i}{\hbar}\ln\int\mathcal{D}\phi\, e^{iS[\phi] + i\int J\phi}$. The one-loop approximation gives:

$$
V_{\text{eff}}(\phi) = V_0(\phi) + \frac{i\hbar}{2}\int\frac{d^4k}{(2\pi)^4}\ln\det\left(-\partial^2 + V''(\phi)\right)
$$

The exported `effective_potential` performs $V_{\text{tree}} + V_\text{1-loop}$ on a grid,
the textbook one-loop resummation.

### Feynman Propagators (free fields)

**Scalar (momentum space):**

$$
\Delta_F(p) = \frac{i}{p^2 - m^2 + i\varepsilon}
$$

**Derivation:** The propagator is the Green's function of the Klein-Gordon operator $(\partial^2 + m^2)\Delta_F(x-y) = -i\delta^4(x-y)$. In momentum space:

$$
(-p^2 + m^2)\tilde{\Delta}_F(p) = -i \implies \tilde{\Delta}_F(p) = \frac{i}{p^2 - m^2}
$$

The $i\varepsilon$ prescription selects the Feynman boundary conditions (positive frequencies propagate forward, negative frequencies backward).

**Position space (Euclidean 4D):**

$$
D(x) = \frac{m}{4\pi^2 |x|}\,K_1(m|x|), \qquad D(x)\big|_{m=0} = \frac{1}{4\pi^2 |x|^2}
$$

**Derivation:** The Euclidean propagator satisfies $(-\partial_E^2 + m^2)D(x) = \delta^4(x)$. Using the Fourier transform in 4D spherical coordinates and the integral representation of the modified Bessel function $K_1$:

$$
D(x) = \int\frac{d^4p_E}{(2\pi)^4}\frac{e^{ip_E\cdot x}}{p_E^2 + m^2} = \frac{m}{4\pi^2|x|}K_1(m|x|)
$$

where $K_1$ is the modified Bessel function of the second kind.

---

## Routines

| Routine | Formula | Derivation |
|---------|---------|------------|
| `pauli_matrices()` | $\sigma_x,\sigma_y,\sigma_z$ | Generators of $\mathfrak{su}(2)$ |
| `spin_operators(s)` | $S_a = \frac{\hbar}{2}\sigma_a$ | Canonical commutation relations |
| `density_matrix(psi)` | $\rho = |\psi\rangle\langle\psi|$ | Outer product |
| `partial_trace(rho, dims, keep)` | $\operatorname{tr}_B \rho$ | Sum over basis of B |
| `purity(rho)` | $\gamma = \operatorname{tr}\rho^2$ | Entanglement criterion |
| `inner_product(cre, him, kre, kim)` | $\langle\psi|\phi\rangle$ | Sesquilinearity |
| `projector(psi)` | $P = |\psi\rangle\langle\psi|$ | Idempotence |
| `coherent_state(alpha, n_max)` | $|\alpha\rangle$ | Displacement operator |
| `number_state(n, dim)` | $|n\rangle$ | Fock basis |
| `rayleigh_schrodinger_1st/2nd(...)` | $E_n^{(1)}, E_n^{(2)}, |n^{(1)}\rangle$ | Perturbation expansion |
| `discretized_path_integral(...)` | $\langle x_f|x_i\rangle$ via slices | Trotter decomposition |
| `effective_potential(phi, v_tree, loop)` | $V_{\text{eff}}$ grid | Coleman-Weinberg |
| `wick_rotation(x, v, m, dt)` | $S_E$ | Analytic continuation |
| `feynman_propagator_scalar(p², m, ε)` | $\frac{i}{p^2-m^2+i\varepsilon}$ | Klein-Gordon Green's function |
| `free_propagator_position(r, m)` | $\frac{m}{4\pi^2 r}K_1(mr)$ | Euclidean Fourier transform |
| `free_propagator_momentum(p₀, p⃗², m, ε)` | $\frac{i}{p_0^2-|\vec p|^2-m^2+i\varepsilon}$ | Feynman prescription |

---

## Usage Examples

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
# Overlap of coherent states: <α|β> = exp(-|α-β|²/2) e^{i·arg}
print(f"number of components: {len(alpha)}  (Poisson distribution)")
```

### Second-order Stark shift

```python
from pysicrs import rayleigh_schrodinger_2nd
# Pass H0 eigenvalues, H' matrix, index n → returns E_n^(2)
```

---

## Advantages & Limitations

✅ Coherent states & Pauli algebra cover spin/optical textbook problems

✅ Path integral and propagators exported directly as amplitude-value functions

✅ Warps nicely with the pure-Rust FFT for split-step evolution (see PDE)

❌ No full state-vector evolution for large $N$ (no qiskit-style gate simulator)

❌ Perturbation module assumes non-degenerate unperturbed energies

❌ Propagators are the free scalar-field forms only (no interactions/couplings)

---

## References

1. Sakurai, J.J. (1994). *Modern Quantum Mechanics*, rev. ed. Addison-Wesley.
2. Shankar, R. (1994). *Principles of Quantum Mechanics*, 2nd ed. Plenum.
3. Weinberg, S. (1995). *The Quantum Theory of Fields*, Vol. 1. Cambridge.
4. Feynman, R.P. (1948). "Space-time approach to non-relativistic quantum mechanics." *Rev. Mod. Phys.* 20:367.
5. Glauber, R.J. (1963). "Coherent and incoherent states of the radiation field." *Phys. Rev.* 131:2766.
6. Coleman, S. & Weinberg, E. (1973). "Radiative corrections as the origin of spontaneous symmetry breaking." *Phys. Rev. D* 7:1888.

---

## Related Topics

- [PDE](pde.md) – the split-step Schrödinger solver drives time evolution
- [Gauge](gauge.md) – $SU(2)$ generators are the Pauli matrices here
- [Special Functions](special_functions.md) – Legendre/spherical harmonics enter angular parts