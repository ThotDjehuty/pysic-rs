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

**Proven identities:**

- $\sigma_a\sigma_b = \delta_{ab} + i\varepsilon_{abc}\sigma_c$
- $[\sigma_a, \sigma_b] = 2i\varepsilon_{abc}\sigma_c$
- $\{\sigma_a, \sigma_b\} = 2\delta_{ab}$
- $(\hat n\cdot\vec\sigma)^2 = I$, so $e^{i\theta \hat n\cdot\vec\sigma} = \cos\theta\,I + i\sin\theta\,\hat n\cdot\vec\sigma$

**Spin operators** are the doubled matrices $S_a = \frac{\hbar}{2}\sigma_a$; they obey the
canonical $\mathfrak{su}(2)$ commutation relations $[S_a,S_b]=i\hbar\varepsilon_{abc}S_c$.

### Density Matrices & Partial Trace

A pure state has $\rho = |\psi\rangle\langle\psi|$, with $\operatorname{tr}\rho = 1$ and
**purity** $\gamma = \operatorname{tr}(\rho^2) = 1$. Mixed states have $\gamma < 1$.

**Partial trace** over system B:

$$
\rho_A = \operatorname{tr}_B\,\rho_{AB} = \sum_j \langle j|_B\ \rho_{AB}\ |j\rangle_B
$$

**Proven EE criterion:** a bipartite state $\rho_{AB}$ is *entangled* iff its reduced density
matrix $\rho_A$ is mixed (purity $<1$).

### Hilbert Space

- **Inner product** $\langle\psi|\phi\rangle = \sum_i \psi_i^*\,\phi_i$.
- **Projector** $P = |\psi\rangle\langle\psi|$ (idempotent: $P^2 = P$).
- **Coherent state** of the harmonic oscillator (eigenstate of the annihilation operator,
  minimum-uncertainty, Glauber 1963):

$$
|\alpha\rangle = e^{-|\alpha|^2/2}\sum_{n=0}^\infty \frac{\alpha^n}{\sqrt{n!}}|n\rangle, \qquad
\hat a|\alpha\rangle = \alpha|\alpha\rangle, \qquad
|\alpha\rangle = \hat D(\alpha)|0\rangle,\ \hat D(\alpha)=e^{\alpha\hat a^\dagger - \alpha^*\hat a}
$$

- **Number state** $|n\rangle$, $\hat N|n\rangle = n|n\rangle$.

### Rayleigh–Schrödinger Perturbation Theory

For $H = H_0 + \lambda H'$ with non-degenerate unperturbed spectrum $E_n^{(0)}$, $|n^{(0)}\rangle$:

$$
E_n^{(1)} = \langle n^{(0)}|H'|n^{(0)}\rangle
$$

$$
E_n^{(2)} = \sum_{m\neq n}\frac{\left|\langle m^{(0)}|H'|n^{(0)}\rangle\right|^2}{E_n^{(0)} - E_m^{(0)}}
$$

$$
|n^{(1)}\rangle = \sum_{m\neq n}\frac{\langle m^{(0)}|H'|n^{(0)}\rangle}{E_n^{(0)} - E_m^{(0)}}|m^{(0)}\rangle
$$

**Proven convergence** requires the perturbation to be small compared with the level spacing.

### Path Integral (Feynman 1948)

The transition amplitude is a sum over all paths weighted by $e^{iS/\hbar}$:

$$
\langle x_f|e^{-i\hat H T/\hbar}|x_i\rangle = \int \mathcal D x\; e^{\frac{i}{\hbar}\int_0^T
\left[\tfrac12 m\dot x^2 - V(x)\right]dt}
$$

**Discretization:** with $N$ time slices $\Delta t = T/N$ the free-particle prefactor is
$\left(\frac{m}{2\pi i\hbar\Delta t}\right)^{1/2}$ — the product of free propagators times
potential weights $\exp(-iV\Delta t/\hbar)$ converges (proven) to the amplitude as $N\to\infty$.

**Wick rotation** $t\to -i\tau$ maps this to the Euclidean (imaginary-time) path integral with
$S_E = \int\left[\tfrac12 m(\tfrac{dx}{d\tau})^2 + V\right]d\tau$, making the integrand real
and sharply peaked at the classical minimum — the basis of instanton and ground-state methods.

### Effective Potential

The one-loop effective potential (Coleman–Weinberg) is the Legendre transform of the free
energy $\Gamma[\phi] = -W[J] + \int J\phi$, summed:

$$
V_{\text{eff}}(\phi) = V_0(\phi) + \frac{\hbar^2}{64\pi^2}\left[\ldots\right]\quad\text{(1-loop)}
$$

The exported `effective_potential` performs $V_{\text{tree}} + V_\text{1-loop}$ on a grid,
the textbook one-loop resummation.

### Feynman Propagators (free fields)

**Scalar (momentum space):**

$$
\Delta_F(p) = \frac{i}{p^2 - m^2 + i\varepsilon}
$$

**Position space (Euclidean 4D):**

$$
D(x) = \frac{m}{4\pi^2 |x|}\,K_1(m|x|), \qquad D(x)\big|_{m=0} = \frac{1}{4\pi^2 |x|^2}
$$

---

## Routines

| Routine | Formula |
|---------|---------|
| `pauli_matrices()` | $\sigma_x,\sigma_y,\sigma_z$ |
| `spin_operators(s)` | $S_a = \frac{\hbar}{2}\sigma_a$ |
| `density_matrix(psi)` | $\rho = |\psi\rangle\langle\psi|$ |
| `partial_trace(rho, dims, keep)` | $\operatorname{tr}_B \rho$ |
| `purity(rho)` | $\gamma = \operatorname{tr}\rho^2$ |
| `inner_product(cre, him, kre, kim)` | $\langle\psi|\phi\rangle$ |
| `projector(psi)` | $P = |\psi\rangle\langle\psi|$ |
| `coherent_state(alpha, n_max)` | $|\alpha\rangle$ |
| `number_state(n, dim)` | $|n\rangle$ |
| `rayleigh_schrodinger_1st/2nd(...)` | $E_n^{(1)}, E_n^{(2)}, |n^{(1)}\rangle$ |
| `discretized_path_integral(...)` | $\langle x_f|x_i\rangle$ via slices |
| `effective_potential(phi, v_tree, loop)` | $V_{\text{eff}}$ grid |
| `wick_rotation(x, v, m, dt)` | $S_E$ |
| `feynman_propagator_scalar(p², m, ε)` | $\frac{i}{p^2-m^2+i\varepsilon}$ |
| `free_propagator_position(r, m)` | $\frac{m}{4\pi^2 r}K_1(mr)$ |
| `free_propagator_momentum(p₀, p⃗², m, ε)` | $\frac{i}{p_0^2-|\vec p|^2-m^2+i\varepsilon}$ |

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

---

## Related Topics

- [PDE](pde.md) – the split-step Schrödinger solver drives time evolution
- [Gauge](gauge.md) – $SU(2)$ generators are the Pauli matrices here
- [Special Functions](special_functions.md) – Legendre/spherical harmonics enter angular parts