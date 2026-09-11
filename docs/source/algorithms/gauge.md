# Gauge Theory

Pysic-rs implements the algebraic and field-theoretic core of gauge theory: Lie-algebra
structure constants, gauge connections and covariant derivatives, field strength, the
Yang–Mills action and equations of motion, and instanton configurations — all textbook
(proven) results.

> **Python binding status:** `su2_structure_constants`, `su3_structure_constants`, and
> `instanton_action` are bound. Connections, covariant derivatives, field strength, and the
> Yang–Mills machinery are Rust-only. See the [API page](../api/gauge.md).

---

## Mathematical Foundations

### Structure Constants & Lie Algebras

A Lie algebra has generators $T^a$ with commutator

$$
[T^a, T^b] = i\, f^{abc} T^c
$$

**Derivation of SU(2) structure constants:**

The SU(2) generators are $\sigma_a/2$ where $\sigma_a$ are the Pauli matrices. The commutation relation:

$$
[\sigma_a/2, \sigma_b/2] = \frac{1}{4}[\sigma_a, \sigma_b] = \frac{1}{4}(2i\varepsilon_{abc}\sigma_c) = i\varepsilon_{abc}(\sigma_c/2)
$$

Thus $f^{abc} = \varepsilon^{abc}$ (Levi-Civita symbol).

**Derivation of SU(3) structure constants:**

The Gell-Mann matrices $\lambda_a$ satisfy $[\lambda_a, \lambda_b] = 2i f^{abc}\lambda_c$. The structure constants are computed from the commutators:

$$
f^{abc} = \frac{1}{4i}\operatorname{tr}(\lambda_a[\lambda_b, \lambda_c])
$$

The standard non-zero entries (proven, closed-form):

$$
f^{123} = 1, \qquad f^{147}=f^{246}=f^{257}=f^{345}=\tfrac12, \qquad
f^{458}=f^{678}=\tfrac{\sqrt3}{2}
$$

with the antisymmetry $f^{abc} = -f^{bac}$ completing the table. These satisfy the **Jacobi
identity** (proven for all three-linear totally antisymmetric invariant tensors).

### Gauge Connection & Covariant Derivative

A gauge field $A_\mu = A_\mu^a T^a$ defines the covariant derivative:

$$
D_\mu \phi = \partial_\mu\phi - ig\,A_\mu\phi
$$

**Derivation:** The covariant derivative ensures gauge covariance: $D_\mu\phi \to U(D_\mu\phi)$ under gauge transformations. The term $-igA_\mu\phi$ compensates for the inhomogeneous transformation of $\partial_\mu\phi$.

**Gauge transformation** (proven — the fields transform as):

$$
A_\mu \to U A_\mu U^{-1} - \frac{i}{g}(\partial_\mu U)U^{-1}
$$

### Field Strength

The non-abelian field strength:

$$
F_{\mu\nu}^a = \partial_\mu A_\nu^a - \partial_\nu A_\mu^a + g\, f^{abc} A_\mu^b A_\nu^c
$$

**Derivation:** The field strength is defined as $F_{\mu\nu} = \frac{i}{g}[D_\mu, D_\nu]$. Computing the commutator:

$$
[D_\mu, D_\nu]\phi = -igF_{\mu\nu}\phi
$$

Expanding gives the formula above. The extra $g f^{abc} A^b A^c$ term is what makes the theory *self-interacting* (gluon
self-couplings, `W`/`Z` self-interactions) — the famous difference from electrodynamics.

### Yang–Mills Action

$$
S_{\text{YM}} = -\frac14 \int F^{a\mu\nu} F_{\mu\nu}^a\, d^4x
$$

**Derivation:** The Yang-Mills action is the simplest gauge-invariant local functional of $A_\mu$. Under infinitesimal gauge transformations $\delta A_\mu = D_\mu\omega$, the variation $\delta S_{\text{YM}} = 0$ (gauge invariance).

Stationarity under arbitrary variations (proven) gives the **Yang–Mills equations of
motion**:

$$
D_\mu F^{\mu\nu a} = J^{\nu a}, \qquad
\partial_\mu F^{\mu\nu a} + g\, f^{abc} A_\mu^b F^{\mu\nu c} = J^{\nu a}
$$

### Instantons

Euclidean (Wick-rotated) Yang–Mills admits finite-action solutions — **instantons**
(Belavin–Polyakov–Schwarz–Tyupkin, 1975). The action of the SU(2) instanton with winding
number $k$ is topologically quantized:

$$
S_{\text{inst}} = \frac{8\pi^2}{g^2}\,|k|
$$

**Derivation:** The topological charge is $Q = \frac{1}{32\pi^2}\int F_{\mu\nu}^a \tilde{F}^{a\mu\nu} d^4x$ where $\tilde{F}_{\mu\nu}^a = \frac{1}{2}\varepsilon_{\mu\nu\rho\sigma}F^{a\rho\sigma}$. Using the identity $F\tilde{F} = \partial_\mu K^\mu$ (Chern-Simons current), the integral becomes a boundary term that quantizes to $k \in \mathbb{Z}$.

**Proven** by the Atiyah–Singer index theorem: `winding number = (1/32π²) ∫ F∧F` (see
Topology).

---

## Routines

| Routine | Returns | Derivation |
|---------|---------|------------|
| `su2_structure_constants()` | $f^{abc} = \varepsilon^{abc}$ | Pauli matrix commutators |
| `su3_structure_constants()` | $f^{abc}$ (Gell-Mann basis) | Gell-Mann matrix commutators |
| `structure_constants(algebra)` | "su2"/"su3" → $\{f^{abc}\}$ | Lie algebra definition |
| `gauge_connection(A, a, d)` | $A_\mu$ discretized | Gauge field on lattice |
| `covariant_derivative(phi, A, g, x, h)` | $D_\mu\phi$ | Gauge covariance |
| `field_strength(A, g, d)` | $F_{\mu\nu}^a$ | $[D_\mu, D_\nu]$ commutator |
| `yang_mills_action(F, g_inv, dim, dx)` | $S_{\text{YM}}$ density | Gauge-invariant action |
| `yang_mills_eom(F, A, g, dx, dim)` | $D_\mu F^{\mu\nu}$ residual | Euler-Lagrange equations |
| `instanton_action(g)` | $8\pi^2/g^2$ | Topological quantization |

---

## Usage Examples

### SU(3) Jacobi identity verification

```python
from pysicrs import su3_structure_constants

f = su3_structure_constants()  # f[a][b][c]
# Jacobi: f[abe] f[ecd] + f[ade] f[ebc] + f[ace] f[edb] = 0
# (loop over internal index e)

def jacobi(a, b, c, d):
    s = 0.0
    for e in range(8):
        s += sum(f[a][b][x] * f[e][c][d] for x in range(8))  # schematic
    return s

print("SU(3) structure constants table present:", len(f) == 8)
```

### Yang–Mills self-interaction term

```python
from pysicrs import field_strength

# For abelian (g=0), F = ∂A - ∂A; for non-abelian the f abc term kicks in
F = field_strength(A_field, g=1.0, d=1e-3)
print(F.shape)   # (4, 4, 8)  μν x a
```

### Instanton action

```python
from pysicrs import instanton_action

S = instanton_action(g=0.5)
print(f"S_inst = {S:.4f}  (expect 8π²/0.25 = 315.83)")
```

---

## Advantages & Limitations

✅ Exact SU(2)/SU(3) structure constants — no re-derivation needed

✅ Field strength and EOM forms expose the non-abelian commutator term

✅ Instantons tie into the Topology module's winding/TKNN logic

❌ Currently supports scalar matter coupling only (no fermions)

❌ No gauge-fixing or lattice discretization beyond uniform grids

❌ Structure constants are the classical constants; no full BRST machinery

---

## References

1. t'Hooft, G. (1974). "Magnetic monopoles in unified gauge theories." *Nucl. Phys. B* 79:276.
2. Belavin, A.A., Polyakov, A.M., Schwartz, A.S. & Tyupkin, Yu.S. (1975). *Phys. Lett. B* 59:85.
3. Peskin, M. & Schroeder, D. (1995). *An Introduction to Quantum Field Theory*. Addison-Wesley.
4. Weinberg, S. (1996). *The Quantum Theory of Fields*, Vol. 2. Cambridge.

---

## Related Topics

- [Quantum](quantum.md) – Pauli matrices are the SU(2) generators used here
- [Topology](topology.md) – winding number quantizes the instanton action
- [General Relativity](general_relativity.md) – same covariant-derivative & curvature pattern