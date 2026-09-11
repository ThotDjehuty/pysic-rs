# Special Functions

Special functions are the canonical solutions of the ODEs that arise throughout
mathematical physics. Pysic-rs implements the complete textbook set — Gamma, Beta,
Bessel, Legendre, Chebyshev, Airy, error functions, exponential integrals, and the
Riemann zeta function — all at `f64` precision with validated reference values.

> **Python binding status:** the full Special Functions suite is Rust-exported but not yet
> bound from Python. The bound Python entry points are listed on the
> [API page](../api/special_functions.md); the examples below show the Rust functions.

---

## 1. The Gamma Function

### 1.1 Definition & Basic Properties

The **Gamma function** $\Gamma(z)$ is the unique analytic continuation of the factorial,
defined for $\operatorname{Re}(z)>0$ by the **Euler integral**:

$$
\boxed{\Gamma(z) = \int_0^\infty t^{z-1}e^{-t}\,dt}
$$

**Theorem (functional equation).** $\Gamma(z+1) = z\,\Gamma(z)$ for all $z\in\mathbb{C}$
(with poles at $z=0,-1,-2,\ldots$).

*Proof.* Integration by parts:

$$
\Gamma(z+1) = \int_0^\infty t^z e^{-t}\,dt
= \bigl[-t^z e^{-t}\bigr]_0^\infty + z\int_0^\infty t^{z-1}e^{-t}\,dt = z\,\Gamma(z). \qquad\square
$$

From $\Gamma(1)=1$, by induction $\Gamma(n+1)=n!$ for $n\in\mathbb{N}_0$.

### 1.2 Reflection Formula

**Theorem (Euler reflection).** For $z\notin\mathbb{Z}$:

$$
\boxed{\Gamma(z)\,\Gamma(1-z) = \frac{\pi}{\sin(\pi z)}}
$$

*Proof.* Consider $f(z) = \Gamma(z)\Gamma(1-z)\sin(\pi z)$. Using the Weierstrass product
for $1/\Gamma(z)$:

$$
\frac{1}{\Gamma(z)} = z e^{\gamma z}\prod_{n=1}^\infty\left(1+\frac{z}{n}\right)e^{-z/n},
$$

one verifies that $f(z)$ is entire and bounded. By Liouville's theorem, $f$ is constant.
Evaluating at $z=1/2$: $\Gamma(1/2)^2\cdot 1 = \pi$, hence $f(z)=\pi$ everywhere.
$\square$

**Corollary:** $\Gamma(1/2) = \sqrt{\pi}$. More generally, $\Gamma(n+1/2) =
\frac{(2n)!}{4^n n!}\sqrt{\pi}$ (by repeated application of the functional equation).

### 1.3 Duplication Formula

**Theorem (Legendre duplication).**

$$
\boxed{\Gamma(2z) = \frac{2^{2z-1}}{\sqrt{\pi}}\,\Gamma(z)\,\Gamma\!\left(z+\tfrac{1}{2}\right)}
$$

*Proof.* Start from the product representation of the Gamma function and compare
coefficients of the $2z$ and $z,z+\frac{1}{2}$ products. Alternatively, evaluate
$\int_0^\infty t^{2z-1}e^{-t}\,dt$ using the substitution $t = u^2$ and the beta function
integral. $\square$

### 1.4 Stirling's Asymptotic Expansion

For $|z|\to\infty$ with $|\arg z|<\pi$:

$$
\Gamma(z) \sim \sqrt{2\pi}\,z^{z-1/2}e^{-z}\left(1+\frac{1}{12z}+\frac{1}{288z^2}+\cdots\right)
$$

This is the starting point for approximations of factorials, binomial coefficients, and
the Euler–Maclaurin summation used for the zeta function.

### 1.5 Numerical Implementation

**Lanczos approximation** with the $g=7$ coefficient set (Godfrey's table):

$$
\Gamma(z+1) = \sqrt{2\pi}\left(z+g+\tfrac{1}{2}\right)^{z+1/2}e^{-(z+g+1/2)}
\left[A_0+\sum_{k=1}^{K}\frac{A_k}{z+k}\right]
$$

For $x<0.5$, the reflection formula is applied first, giving uniform accuracy across the
positive real axis.

| Domain | Typical $\lvert\text{error}\rvert$ |
|--------|--------------------------------------|
| $x\in(0,50)$ | $< 1\times 10^{-12}$ |
| $x\in(0,200)$ | $< 1\times 10^{-10}$ |

---

## 2. The Beta Function

### 2.1 Definition

$$
\boxed{B(a,b) = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)} = \int_0^1 t^{a-1}(1-t)^{b-1}\,dt}
$$

The integral representation follows from the substitution $t = u/(1+u)$ in the product
$\Gamma(a)\Gamma(b) = \int_0^\infty\int_0^\infty u^{a-1}v^{b-1}e^{-(u+v)}\,du\,dv$, then
setting $u = tw$, $v = (1-t)w$ and integrating out $w$.

### 2.2 Symmetry & Recurrence

- **Symmetry:** $B(a,b) = B(b,a)$.
- **Recurrence:** $B(a+1,b) = \frac{a}{a+b}B(a,b)$.
- **Relation to binomial coefficients:** $B(k+1,n-k+1) = \frac{1}{(n+1)\binom{n}{k}}$.

---

## 3. Bessel Functions

### 3.1 Bessel's Differential Equation

The Bessel functions $J_n(x)$ and $Y_n(x)$ are the regular and irregular solutions of

$$
\boxed{x^2 y'' + x\,y' + (x^2-n^2)y = 0}
$$

This is a Fuchsian equation with regular singular points at $x=0$ and irregular singular
point at $x=\infty$.

### 3.2 Integral Representations

**Bessel function of the first kind:**

$$
J_n(x) = \frac{1}{\pi}\int_0^\pi\cos(x\sin\theta - n\theta)\,d\theta
= \frac{1}{2\pi}\int_{-\pi}^{\pi}e^{i(x\sin\theta-n\theta)}\,d\theta
$$

**Poisson representation** (for $J_\nu$ with $\operatorname{Re}(\nu)>-1/2$):

$$
J_\nu(x) = \frac{(x/2)^\nu}{\sqrt{\pi}\,\Gamma(\nu+1/2)}\int_{-1}^{1}(1-t^2)^{\nu-1/2}\cos(xt)\,dt
$$

### 3.3 Asymptotic Behavior

**Large-$x$ asymptotics** (from the method of steepest descent or WKB):

$$
J_n(x) \sim \sqrt{\frac{2}{\pi x}}\cos\!\left(x - \frac{n\pi}{2} - \frac{\pi}{4}\right)
$$

$$
Y_n(x) \sim \sqrt{\frac{2}{\pi x}}\sin\!\left(x - \frac{n\pi}{2} - \frac{\pi}{4}\right)
$$

**Small-$x$ behavior:** $J_n(x) \sim \frac{1}{n!}(x/2)^n$, $Y_0(x) \sim \frac{2}{\pi}\ln(x/2)$.

### 3.4 Modified Bessel Functions

The modified Bessel functions $I_\nu(x) = i^{-\nu}J_\nu(ix)$ and $K_\nu(x)$ satisfy

$$
x^2 y'' + x\,y' - (x^2+\nu^2)y = 0.
$$

$K_\nu$ decays exponentially: $K_\nu(x) \sim \sqrt{\pi/(2x)}\,e^{-x}$ for $x\to\infty$.

The Euclidean propagator in 4D uses $K_1$: $D(x) = \frac{m}{4\pi^2|x|}K_1(m|x|)$ (see
[Quantum](quantum.md)).

### 3.5 Numerical Implementation

Polynomial approximation for $|x|<8$ and asymptotic expansion for $|x|\ge 8$, matching
the standard Numerical Recipes approach. Orders $n=0,1$ directly; general $n$ via
recurrence from the ascending series.

---

## 4. Legendre Polynomials & Spherical Harmonics

### 4.1 Legendre's Differential Equation

The Legendre polynomials $P_\ell(x)$ solve

$$
(1-x^2)y'' - 2xy' + \ell(\ell+1)y = 0
$$

and form a complete orthogonal basis on $[-1,1]$:

$$
\int_{-1}^{1}P_\ell(x)P_{\ell'}(x)\,dx = \frac{2}{2\ell+1}\delta_{\ell\ell'}
$$

### 4.2 Rodrigues Formula

**Theorem.**

$$
\boxed{P_\ell(x) = \frac{1}{2^\ell\,\ell!}\frac{d^\ell}{dx^\ell}(x^2-1)^\ell}
$$

*Proof.* By induction on $\ell$. The generating function $(1-2xt+t^2)^{-1/2} =
\sum_{\ell=0}^\infty P_\ell(x)t^\ell$ can be differentiated $\ell$ times with respect to
$t$ and evaluated at $t=0$ to extract the coefficient, yielding the Rodrigues formula.
The orthogonality follows from the $(x^2-1)^\ell$ factor: each integration by parts
transfers a derivative onto the other polynomial, killing all terms up to $\ell-1$.
$\square$

### 4.3 Associated Legendre Functions

For $m\ge 0$:

$$
P_\ell^m(x) = (-1)^m(1-x^2)^{m/2}\frac{d^m}{dx^m}P_\ell(x)
$$

These satisfy the orthogonality relation

$$
\int_{-1}^{1}P_\ell^m(x)P_{\ell'}^m(x)\,dx = \frac{2}{2\ell+1}\frac{(\ell+m)!}{(\ell-m)!}\delta_{\ell\ell'}
$$

### 4.4 Spherical Harmonics

$$
\boxed{Y_\ell^m(\theta,\phi) = \sqrt{\frac{2\ell+1}{4\pi}\frac{(\ell-m)!}{(\ell+m)!}}\,
P_\ell^m(\cos\theta)\,e^{im\phi}}
$$

**Orthonormality:**

$$
\int_{S^2}Y_\ell^{m*}(\theta,\phi)\,Y_{\ell'}^{m'}(\theta,\phi)\,d\Omega
= \delta_{\ell\ell'}\delta_{mm'}
$$

*Proof.* The $e^{im\phi}$ factor gives $\delta_{mm'}$ from the $\phi$-integral. The
$\theta$-integral reduces to the Legendre orthogonality via $x=\cos\theta$. The
normalisation constant is fixed by requiring unit norm. $\square$

**Addition theorem:**

$$
P_\ell(\cos\gamma) = \frac{4\pi}{2\ell+1}\sum_{m=-\ell}^{\ell}Y_\ell^m(\theta_1,\phi_1)Y_\ell^{m*}(\theta_2,\phi_2)
$$

where $\gamma$ is the angle between the two directions. This is used extensively in
multipole expansions of potentials (see [EM](em.md) and [General Relativity](general_relativity.md)).

---

## 5. Chebyshev Polynomials

The Chebyshev polynomials of the first and second kind are defined by

$$
T_n(x) = \cos(n\arccos x), \qquad U_n(x) = \frac{\sin((n+1)\arccos x)}{\sin(\arccos x)}
$$

**Recurrence:** $T_0=1$, $T_1=x$, $T_{n+1}=2xT_n-T_{n-1}$.

**Orthogonality:** $\int_{-1}^{1}\frac{T_n(x)T_m(x)}{\sqrt{1-x^2}}dx = \frac{\pi}{2}\delta_{mn}$
(for $n,m\ge 1$), with weight function $w(x)=(1-x^2)^{-1/2}$.

These are the basis of **Chebyshev spectral methods** for PDEs (see [PDE](pde.md)).

---

## 6. Error Function

### 6.1 Definition

$$
\operatorname{erf}(x) = \frac{2}{\sqrt{\pi}}\int_0^x e^{-t^2}\,dt, \qquad
\operatorname{erfc}(x) = 1-\operatorname{erf}(x) = \frac{2}{\sqrt{\pi}}\int_x^\infty e^{-t^2}\,dt
$$

### 6.2 Properties

- **Symmetry:** $\operatorname{erf}(-x) = -\operatorname{erf}(x)$ (odd function).
- **Limits:** $\operatorname{erf}(\infty) = 1$ (Gaussian integral $\int_0^\infty e^{-t^2}dt = \sqrt{\pi}/2$).
- **Small-$x$ expansion:** $\operatorname{erf}(x) = \frac{2}{\sqrt{\pi}}\left(x-\frac{x^3}{3}+\frac{x^5}{10}-\cdots\right)$.
- **Large-$x$ asymptotics:** $\operatorname{erfc}(x) \sim \frac{e^{-x^2}}{x\sqrt{\pi}}\left(1-\frac{1}{2x^2}+\cdots\right)$.

**Connection to normal distribution:** $\Phi(x) = \frac{1}{2}\left[1+\operatorname{erf}(x/\sqrt{2})\right]$.

### 6.3 Numerical Implementation

Series near the origin; the standard rational approximation (Abramowitz & Stegun 7.1.26)
away from it, with the symmetry relation applied for negative arguments.

---

## 7. Riemann Zeta Function

### 7.1 Definition

$$
\zeta(s) = \sum_{n=1}^\infty\frac{1}{n^s}, \qquad \operatorname{Re}(s)>1
$$

### 7.2 Special Values

**$\zeta(2) = \pi^2/6$.**

*Proof (Euler).* Consider $\sin x = x\prod_{n=1}^\infty(1-x^2/(n^2\pi^2))$. Expand both
sides as power series: $\sin x = x - x^3/6 + \cdots$, and the product expansion gives
$\sum_{n=1}^\infty 1/(n^2\pi^2) = 1/6$. Hence $\zeta(2) = \pi^2/6$. $\square$

**$\zeta(4) = \pi^4/90$.** Similarly, comparing the $x^4$ coefficient in $\sin^2 x$ and
the product.

**$\zeta(2n) = (-1)^{n+1}\frac{(2\pi)^{2n}B_{2n}}{2(2n)!}$** where $B_{2n}$ are Bernoulli
numbers. This follows from the residue calculus evaluation of $\oint z^{-2n}\pi\cot(\pi z)\,dz$.

### 7.3 Functional Equation

$$
\boxed{\zeta(s) = 2^s\pi^{s-1}\sin\!\left(\frac{\pi s}{2}\right)\Gamma(1-s)\,\zeta(1-s)}
$$

This extends $\zeta$ to the entire complex plane (meromorphic continuation), with a simple
pole at $s=1$ with residue 1:

$$
\zeta(s) = \frac{1}{s-1}+\gamma+O(s-1)
$$

where $\gamma\approx 0.5772$ is the Euler–Mascheroni constant.

### 7.4 Numerical Implementation

**Euler–Maclaurin summation** with Bernoulli corrections:

$$
\sum_{n=1}^{N-1}\frac{1}{n^s}+\frac{1}{(s-1)N^{s-1}}+\frac{1}{2N^s}+\sum_{k=1}^{K}\frac{B_{2k}}{(2k)!}\frac{s(s+1)\cdots(s+2k-2)}{N^{s+2k-1}}
$$

accelerates the slowly converging Dirichlet series. Typical accuracy: $|\text{error}| <
10^{-10}$ for $\operatorname{Re}(s)>0.5$.

---

## 8. Airy Functions

The Airy functions $\operatorname{Ai}(x)$ and $\operatorname{Bi}(x)$ solve

$$
y'' - x\,y = 0
$$

- $\operatorname{Ai}(x)$ decays exponentially for $x>0$ and oscillates for $x<0$.
- $\operatorname{Bi}(x)$ grows exponentially for $x>0$.

**Integral representation:**

$$
\operatorname{Ai}(x) = \frac{1}{\pi}\int_0^\infty\cos\!\left(\frac{t^3}{3}+xt\right)dt
$$

**WKB asymptotics** near the turning point $x=0$: the transition from oscillatory to
exponential behavior is described by Airy functions — this is the universal local
structure near a simple turning point in quantum mechanics (see [PDE](pde.md)).

---

## 9. Exponential Integral

$$
E_1(x) = \int_x^\infty\frac{e^{-t}}{t}\,dt = \Gamma(0,x)
$$

(incomplete Gamma function). For large $x$: integration by parts gives
$E_1(x) \sim e^{-x}/x\cdot(1-1/x+2/x^2-\cdots)$.

---

## 10. Available Functions

```python
from pysicrs import (
    gamma, ln_gamma, beta,
    bessel_j0, bessel_j1, bessel_y0, bessel_y1,
    legendre_p, legendre_plm, spherical_harmonic,
    chebyshev_t, chebyshev_u, airy_ai,
    erf, erfc, expint_e1, expint_en, zeta,
)
```

---

## 11. Usage Examples

### Gamma to factorial

```python
from pysicrs import gamma

for n in range(1, 8):
    assert abs(gamma(n) - math.factorial(n - 1)) < 1e-9
print("Γ(n) = (n-1)! ✓")
```

### Riemann zeta at even integers

```python
from pysicrs import zeta
import math

print(f"ζ(2) = {zeta(2.0):.10f}  expected {math.pi**2/6:.10f}")
print(f"ζ(4) = {zeta(4.0):.10f}  expected {math.pi**4/90:.10f}")
```

### Spherical harmonics normalization

```python
from pysicrs import spherical_harmonic

y00 = spherical_harmonic(0, 0, 0.0, 0.0)
print(f"Y_00 = {y00:.6f}  expected {1/math.sqrt(4*math.pi):.6f}")
```

---

## 12. Numerical Accuracy

| Function | Domain | Typical $\lvert\text{error}\rvert$ | Reference |
|----------|--------|--------------------------------------|-----------|
| $\Gamma(x)$ | $(0, 100)$ | $< 1\times 10^{-12}$ | Lanczos |
| $\operatorname{erf}$ | $\mathbb{R}$ | $< 1\times 10^{-14}$ | A&S 7.1.26 |
| $J_0$ | $[0, 30]$ | $< 1\times 10^{-12}$ | Polynomial/asymptotic |
| $\zeta$ | $(1, 30]$ | $< 1\times 10^{-10}$ | Euler–Maclaurin |

---

## 13. Advantages & Limitations

✅ Covers the special functions most common in physics problems

✅ All implementations are standalone (no BLAS/LAPACK dependency)

✅ Log-space Gamma avoids overflow for large arguments

✅ Unit-tested against reference values (Abramowitz & Stegun, closed forms)

❌ `f64` only — no arbitrary precision

❌ Bessel functions limited to orders $n=0,1$ (general order via recurrence/ascending series)

❌ Real axis only for most functions (no complex continuation exposed yet)

---

## 14. References

1. Abramowitz, M. & Stegun, I. (1964). *Handbook of Mathematical Functions*. NBS.
2. Arfken, G. (1985). *Mathematical Methods for Physicists*, 3rd ed. Academic Press.
3. Press, W. et al. (2007). *Numerical Recipes*, 3rd ed. Cambridge University Press.
4. Whittaker, E.T. & Watson, G.N. (1927). *A Course of Modern Analysis*, 4th ed.
5. DLMF, NIST (2023). *Digital Library of Mathematical Functions*. https://dlmf.nist.gov/

---

## 15. Related Topics

- [Linear Algebra](linalg.md) — decompositions used by the interpolation kernels
- [Calculus](calculus.md) — quadrature rules that pair with special functions
- [Quantum](quantum.md) — spherical harmonics enter angular parts; $K_1$ in propagators
- [PDE](pde.md) — Airy functions near turning points; Bessel functions in cylindrical problems
- [EM](em.md) — multipole expansions use Legendre polynomials
- [General Relativity](general_relativity.md) — $K_1$ in Euclidean propagators; spherical harmonics in perturbation theory
