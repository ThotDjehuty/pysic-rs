# Gauge Theory

Pysic-rs implements the algebraic and field-theoretic core of gauge theory: Lie-algebra
structure constants, gauge connections and covariant derivatives, field strength, the
Yang–Mills action and equations of motion, and instanton configurations — all at
textbook (proven) level.

> **Python binding status:** `su2_structure_constants`, `su3_structure_constants`, and
> `instanton_action` are bound. Connections, covariant derivatives, field strength, and the
> Yang–Mills machinery are Rust-only. See the [API page](../api/gauge.md).

---

## 1. Lie Algebras & Structure Constants

### 1.1 Definition

A **Lie algebra** $\mathfrak{g}$ is a vector space equipped with a bilinear, antisymmetric
bracket $[\cdot,\cdot]$ satisfying the **Jacobi identity**:

$$
[X,[Y,Z]] + [Y,[Z,X]] + [Z,[X,Y]] = 0.
$$

In a basis $\{T^a\}$, the commutator is

$$
[T^a,T^b] = i\,f^{abc}\,T^c,
$$

where $f^{abc}$ are the **structure constants**, antisymmetric in all three indices.

### 1.2 SU(2) Structure Constants

The SU(2) generators in the fundamental representation are $T_a = \sigma_a/2$ where
$\sigma_a$ are the Pauli matrices (see [Quantum](quantum.md)).

**Theorem.** The structure constants of SU(2) are $f^{abc} = \varepsilon^{abc}$ (the
Levi-Civita symbol).

*Proof.* From $[\sigma_a/2,\sigma_b/2] = \frac{1}{4}[\sigma_a,\sigma_b] = \frac{1}{4}
(2i\varepsilon_{abc}\sigma_c) = i\varepsilon_{abc}(\sigma_c/2)$. Reading off: $f^{abc}
= \varepsilon^{abc}$. $\square$

**Properties:**
- $f^{123} = 1$, $f^{132} = -1$, etc. (6 sign-flipped copies of the basic triple).
- All other components vanish.
- Jacobi identity: $\varepsilon_{abe}\varepsilon_{ecd} + \varepsilon_{ade}\varepsilon_{ebc} +
  \varepsilon_{ace}\varepsilon_{edb} = 0$ (proved via the identity $\varepsilon_{abe}\varepsilon_{cde}
  = \delta_{ac}\delta_{bd} - \delta_{ad}\delta_{bc}$).

### 1.3 SU(3) Structure Constants

The Gell-Mann matrices $\lambda_a$ ($a=1,\ldots,8$) satisfy $[\lambda_a,\lambda_b] =
2if^{abc}\lambda_c$. The structure constants are extracted by

$$
f^{abc} = \frac{1}{4i}\operatorname{tr}\bigl(\lambda_a[\lambda_b,\lambda_c]\bigr).
$$

**Non-zero entries** (proven, closed form):

$$
f^{123} = 1, \qquad f^{147} = f^{246} = f^{257} = f^{345} = \frac{1}{2}, \qquad
f^{458} = f^{678} = \frac{\sqrt{3}}{2}
$$

with the complete antisymmetry $f^{abc} = -f^{bac} = -f^{acb}$ etc. filling out all
non-zero permutations. These are the same structure constants that appear in the
electroweak and strong interaction Lagrangians.

**Jacobi identity:** For SU(3), the identity
$f^{abe}f^{ecd}+f^{ade}f^{ebc}+f^{ace}f^{edb}=0$ follows from the Jacobi identity of
the Lie algebra (applied to the matrices $\lambda_a$). This is the algebraic consistency
condition that ensures gauge invariance of the Yang–Mills action.

### 1.4 Killing Form

The **Killing form** on a Lie algebra is

$$
B(X,Y) = \operatorname{tr}(\operatorname{ad}_X \operatorname{ad}_Y),
$$

where $(\operatorname{ad}_X)(Y) = [X,Y]$. For a simple Lie algebra, $B$ is
non-degenerate (Cartan's criterion). For SU($N$): $B(T^a,T^b) = 2N\delta^{ab}$ (in the
fundamental representation).

---

## 2. Gauge Connection & Covariant Derivative

### 2.1 Gauge Fields

A **gauge field** $A_\mu = A_\mu^a T^a$ is a Lie-algebra-valued 1-form on spacetime. The
**covariant derivative** acting on matter fields $\phi$ in representation $R$ is

$$
\boxed{D_\mu\phi = \partial_\mu\phi - ig\,A_\mu\phi}
$$

where $g$ is the coupling constant and the product $A_\mu\phi$ uses the representation
matrix.

### 2.2 Gauge Transformation

Under a local gauge transformation $U(x) = e^{i\omega^a(x)T^a} \in G$:

$$
\phi \to U\phi, \qquad
A_\mu \to U A_\mu U^{-1} - \frac{i}{g}(\partial_\mu U)U^{-1}
$$

*Proof.* The covariant derivative transforms covariantly: $D_\mu\phi \to U(D_\mu\phi)$,
provided $A_\mu$ transforms as above. The term $-(i/g)(\partial_\mu U)U^{-1}$ compensates
for the inhomogeneous transformation of $\partial_\mu\phi$: $\partial_\mu(U\phi) =
U\partial_\mu\phi + (\partial_\mu U)\phi$, and the gauge term in $D_\mu(U\phi)$ produces
$-ig(UA_\mu U^{-1})U\phi - i(\partial_\mu U)\phi = U(-igA_\mu\phi)$. $\square$

**Infinitesimal form** ($U \approx 1+i\omega$):

$$
\delta A_\mu = \partial_\mu\omega + ig[A_\mu,\omega] = D_\mu\omega
$$

This is an **adjoint** transformation — the gauge field transforms in the adjoint
representation of $G$.

---

## 3. Field Strength

### 3.1 Definition

The **non-abelian field strength** is defined as

$$
\boxed{F_{\mu\nu} = \frac{i}{g}[D_\mu,D_\nu]}
$$

Expanding the commutator:

$$
[D_\mu,D_\nu]\phi = (\partial_\mu\partial_\nu\phi - igA_\mu\partial_\nu\phi - ig\partial_\mu(A_\nu\phi)
- g^2 A_\mu A_\nu\phi) - (\mu\leftrightarrow\nu)
$$

$$
= -ig(\partial_\mu A_\nu - \partial_\nu A_\mu + ig[A_\mu,A_\nu])\phi
$$

Hence

$$
\boxed{F_{\mu\nu}^a = \partial_\mu A_\nu^a - \partial_\nu A_\mu^a + g\,f^{abc}A_\mu^b A_\nu^c}
$$

The extra $gf^{abc}A^bA^c$ term is the **non-abelian self-interaction** — the fundamental
difference from electrodynamics, where $F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu$
(abelian: $[A_\mu,A_\nu]=0$).

### 3.2 Bianchi Identity

$$
D_{[\lambda}F_{\mu\nu]} = 0
$$

*Proof.* From the Jacobi identity of covariant derivatives: $[D_\lambda,[D_\mu,D_\nu]] +
[\text{cyclic}] = 0$, and $[D_\mu,D_\nu] = -igF_{\mu\nu}$. $\square$

This is the **homogeneous** equation of the gauge field (analogous to $\nabla\cdot B=0$ in
electromagnetism).

### 3.3 Transformation Law

$$
F_{\mu\nu} \to U\,F_{\mu\nu}\,U^{-1}
$$

The field strength transforms **covariantly** (in the adjoint representation), confirming
that $F_{\mu\nu}^a F^{a\mu\nu}$ is gauge-invariant.

---

## 4. Yang–Mills Action & Equations of Motion

### 4.1 Yang–Mills Action

$$
\boxed{S_{\rm YM} = -\frac{1}{4}\int d^4x\,F^{a\mu\nu}F_{\mu\nu}^a = -\frac{1}{2}\int d^4x\,\operatorname{tr}(F_{\mu\nu}F^{\mu\nu})}
$$

**Gauge invariance (proven):** Under infinitesimal transformations $\delta A_\mu = D_\mu\omega$:

$$
\delta S_{\rm YM} = -\frac{1}{2}\int\operatorname{tr}(\delta F_{\mu\nu}F^{\mu\nu})
= -\frac{i}{g}\int\operatorname{tr}([D_\mu\omega]F^{\mu\nu}F_{\nu}{}^{\mu}) = 0
$$

after integration by parts and using the antisymmetry of $F_{\mu\nu}$.

### 4.2 Equations of Motion

Stationarity under $\delta A_\mu$ gives the **Yang–Mills equations**:

$$
\boxed{D_\mu F^{\mu\nu} = J^{\nu}}
$$

or in components:

$$
\partial_\mu F^{\mu\nu a} + g\,f^{abc}A_\mu^b F^{\mu\nu c} = J^{\nu a}
$$

The source current is $J^{\nu a} = g\bar{\psi}\gamma^\nu T^a\psi$ for fermionic matter.

### 4.3 Self-Interactions

The non-abelian field strength contains cubic and quartic terms in $A_\mu$:

$$
F_{\mu\nu}^a F^{a\mu\nu} = (\partial_\mu A_\nu^a - \partial_\nu A_\mu^a)^2
+ 2g f^{abc}A_\mu^b A_\nu^c(\partial^\mu A^{a\nu}-\partial^\nu A^{a\mu})
+ g^2 f^{abc}f^{ade}A_\mu^b A_\nu^c A^{\mu d} A^{\nu e}
$$

The last two terms give **3-gluon** and **4-gluon** vertices — the gluons interact with
each other (unlike photons in QED).

---

## 5. Instantons

### 5.1 Euclidean Action & Topological Charge

After Wick rotation $t\to -i\tau$, the Euclidean Yang–Mills action is

$$
S_E = \frac{1}{2g^2}\int d^4x\,\operatorname{tr}(F_{\mu\nu}F_{\mu\nu}) \ge 0
$$

The **topological charge** (Pontryagin index) is

$$
Q = \frac{1}{32\pi^2}\int d^4x\,F_{\mu\nu}^a\tilde{F}^{a\mu\nu}
$$

where $\tilde{F}_{\mu\nu}^a = \frac{1}{2}\varepsilon_{\mu\nu\rho\sigma}F^{a\rho\sigma}$ is the
Hodge dual.

### 5.2 The BPST Instanton

**Theorem (Belavin–Polyakov–Schwarz–Tyupkin, 1975).** The Euclidean YM equations admit
finite-action solutions — **instantons** — with topological charge $Q = k \in \mathbb{Z}$.
The action of the $k$-instanton is

$$
\boxed{S_{\rm inst} = \frac{8\pi^2}{g^2}|k|}
$$

*Proof.* The integrand $F\tilde{F}$ is a total derivative:

$$
F_{\mu\nu}^a\tilde{F}^{a\mu\nu} = \partial_\mu K^\mu, \qquad
K^\mu = 2\varepsilon^{\mu\nu\rho\sigma}\left(A_\nu^a\partial_\rho A_\sigma^a + \frac{g}{3}f^{abc}A_\nu^a A_\rho^b A_\sigma^c\right)
$$

(the Chern–Simons current). By Stokes' theorem, $Q$ becomes a boundary integral at
spatial infinity, which is a map $S^3\to G$. For $G = SU(2)$, this is classified by
$\pi_3(SU(2)) = \mathbb{Z}$, so $Q\in\mathbb{Z}$. The Bogomolny bound gives $S_E \ge
\frac{8\pi^2}{g^2}|Q|$, and instantons saturate this bound. $\square$

### 5.3 Physiological Significance

- **Vacuum structure:** the QCD vacuum is a superposition of topologically distinct vacua
  $|n\rangle$ labelled by the winding number $n\in\mathbb{Z}$.
- **$\theta$-vacuum:** $|\theta\rangle = \sum_n e^{in\theta}|n\rangle$, leading to the
  strong CP problem.
- **Tunneling:** instantons mediate transitions between topologically distinct vacua.
- **Chiral symmetry breaking:** instanton interactions break the $U(1)_A$ anomaly.

---

## 6. Routines

| Routine | Returns | Derivation |
|---------|---------|------------|
| `su2_structure_constants()` | $f^{abc}=\varepsilon^{abc}$ | Pauli matrix commutators |
| `su3_structure_constants()` | $f^{abc}$ (Gell-Mann basis) | Gell-Mann matrix commutators |
| `structure_constants(algebra)` | "su2"/"su3" → $\{f^{abc}\}$ | Lie algebra definition |
| `gauge_connection(A, a, d)` | $A_\mu$ discretized | Gauge field on lattice |
| `covariant_derivative(phi, A, g, x, h)` | $D_\mu\phi$ | Gauge covariance |
| `field_strength(A, g, d)` | $F_{\mu\nu}^a$ | $[D_\mu,D_\nu]$ commutator |
| `yang_mills_action(F, g_inv, dim, dx)` | $S_{\rm YM}$ density | Gauge-invariant action |
| `yang_mills_eom(F, A, g, dx, dim)` | $D_\mu F^{\mu\nu}$ residual | Euler–Lagrange equations |
| `instanton_action(g)` | $8\pi^2/g^2$ | Topological quantization |

---

## 7. Usage Examples

### SU(3) Jacobi identity verification

```python
from pysicrs import su3_structure_constants

f = su3_structure_constants()  # f[a][b][c]

def jacobi(a, b, c, d):
    s = 0.0
    for e in range(8):
        s += sum(f[a][b][x] * f[e][c][d] for x in range(8))
    return s

print("SU(3) structure constants present:", len(f) == 8)
```

### Yang–Mills self-interaction term

```python
from pysicrs import field_strength

F = field_strength(A_field, g=1.0, d=1e-3)
print(F.shape)   # (4, 4, 8)  μν × a
```

### Instanton action

```python
from pysicrs import instanton_action

S = instanton_action(g=0.5)
print(f"S_inst = {S:.4f}  (expect 8π²/0.25 ≈ 315.83)")
```

---

## 8. Advantages & Limitations

✅ Exact SU(2)/SU(3) structure constants — no re-derivation needed

✅ Field strength and EOM forms expose the non-abelian commutator term

✅ Instantons tie into the Topology module's winding/TKNN logic

❌ Currently supports scalar matter coupling only (no fermions)

❌ No gauge-fixing or lattice discretization beyond uniform grids

❌ Structure constants are the classical constants; no full BRST machinery

---

## 9. References

1. t'Hooft, G. (1974). "Magnetic monopoles in unified gauge theories." *Nucl. Phys. B* 79:276.
2. Belavin, A.A., Polyakov, A.M., Schwartz, A.S. & Tyupkin, Yu.S. (1975). *Phys. Lett. B* 59:85.
3. Peskin, M. & Schroeder, D. (1995). *An Introduction to Quantum Field Theory*. Addison-Wesley.
4. Weinberg, S. (1996). *The Quantum Theory of Fields*, Vol. 2. Cambridge.
5. Rajaraman, R. (1982). *Solitons and Instantons*. North-Holland.

---

## 10. Related Topics

- [Quantum](quantum.md) — Pauli matrices are the SU(2) generators used here
- [Topology](topology.md) — winding number quantizes the instanton action
- [General Relativity](general_relativity.md) — same covariant-derivative & curvature pattern
- [Electromagnetism](em.md) — abelian limit of Yang–Mills is Maxwell
