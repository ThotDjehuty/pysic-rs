# Special Functions

**Special functions** are the canonical solutions of the ordinary differential equations that
dominate mathematical physics. Pysic-rs implements the textbook set — Gamma, Beta, Bessel,
Legendre, Chebyshev, Airy, error functions, the exponential integral, and the Riemann zeta
function — all at `f64` precision with validated reference values.

> **Python binding status:** the full Special Functions suite is Rust-exported but not yet
> bound from Python. The bound Python entry points are listed on the
> [API page](../api/special_functions.md); the examples below show the Rust functions.

---

## Mathematical Foundations

### Gamma Function

The Gamma function $\Gamma(z)$ is the analytic continuation of the factorial, and is defined
for $\operatorname{Re}(z)>0$ by the Euler integral:

$$
\Gamma(z) = \int_0^\infty t^{z-1} e^{-t}\,dt
$$

**Proven properties (standard references):**

- $\Gamma(n) = (n-1)!$ for $n \in \mathbb{N}^+$
- $\Gamma(1/2) = \sqrt{\pi}$
- **Functional equation**: $\Gamma(z+1) = z\,\Gamma(z)$
- **Reflection formula**: $\Gamma(z)\,\Gamma(1-z) = \dfrac{\pi}{\sin(\pi z)}$
- **Duplication formula** (Legendre): $\Gamma(z)\,\Gamma(z+\tfrac12) = 2^{1-2z}\sqrt{\pi}\,\Gamma(2z)$

**Numerical implementation:** Lanczos approximation with the $g=7$ coefficient set
(Numerical Recipes / Godfrey's table). For $x<0.5$ the reflection formula is applied first,
which gives uniform accuracy across the positive real axis.

```
Accuracy: |error| < 1e-12 for all x in (0, 50)
```

### Beta Function

$$
B(a,b) = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)} = \int_0^1 t^{a-1}(1-t)^{b-1}\,dt
$$

Computed in log space to avoid overflow for large arguments.

### Bessel Functions

$J_n(x)$ and $Y_n(x)$ are the regular and irregular solutions of Bessel's equation:

$$
x^2 y'' + x y' + (x^2 - n^2)y = 0
$$

**Proven asymptotic behavior**

- $J_0(0)=1$, $Y_0(x)\to -\infty$ as $x\to 0^+$
- As $x\to\infty$: $J_n(x) \sim \sqrt{\tfrac{2}{\pi x}}\cos\left(x - \tfrac{n\pi}{2} - \tfrac{\pi}{4}\right)$
- $Y_n(-x)$ has a logarithmic branch cut along the negative real axis

**Numerical implementation:** polynomial approximation for $|x|<8$ and the asymptotic
expansion for $|x|\ge 8$, matching the standard Numeric Recovery Recipes approach.

### Legendre Polynomials

The Legendre polynomials $P_\ell(x)$ solve Legendre's equation and form a complete orthogonal
basis on $[-1,1]$:

$$
\int_{-1}^{1} P_\ell(x)P_{\ell'}(x)\,dx = \frac{2}{2\ell+1}\,\delta_{\ell\ell'}
$$

**Associated Legendre functions** $P_\ell^m(x)$ (used for spherical harmonics and
gravitational/magnetic multipoles) extend this to $m\neq 0$.

### Spherical Harmonics

$$
Y_\ell^m(\theta,\phi) = \sqrt{\frac{(2\ell+1)}{4\pi}\frac{(\ell-m)!}{(\ell+m)!}}\,
P_\ell^m(\cos\theta)\,e^{im\phi}
$$

**Orthonormality (proven):** $\int Y_{\ell}^{m*} Y_{\ell'}^{m'}\, d\Omega = \delta_{\ell\ell'}\delta_{mm'}$

### Error Function

$$
\operatorname{erf}(x) = \frac{2}{\sqrt{\pi}}\int_0^x e^{-t^2}\,dt, \qquad
\operatorname{erfc}(x) = 1 - \operatorname{erf}(x)
$$

**Proven limits:** $\operatorname{erf}(\infty)=1$, $\operatorname{erf}(-x)=-\operatorname{erf}(x)$,
$\operatorname{erf}(x)\sim \frac{2x}{\sqrt\pi}$ as $x\to 0$.

**Numerical implementation:** series near the origin, the standard rational approximation
(Abramowitz & Stegun 7.1.26) away from it, with the symmetry relation applied for negative
arguments.

### Riemann Zeta Function

$$
\zeta(s) = \sum_{n=1}^\infty \frac{1}{n^s}, \qquad \operatorname{Re}(s)>1
$$

**Proven values:** $\zeta(2)=\pi^2/6$, $\zeta(4)=\pi^4/90$, $\zeta(3)\approx 1.2020569$ (Apéry's
constant — irrational, proven 1978). The pole at $s=1$ and the trivial zeros at negative even
integers are standard results.

**Numerical implementation:** Euler–Maclaurin summation with Bernoulli corrections, giving an
acceleration of the slowly converging Dirichlet series.

### Airy & Exponential Integral

- **Airy**: $y'' - x\,y = 0$ has the independent solutions $\operatorname{Ai}(x)$,
  $\operatorname{Bi}(x)$; $\operatorname{Ai}$ decays exponentially for $x>0$ and
  oscillates for $x<0$.
- **Exponential integral**: $E_1(x)=\displaystyle\int_x^\infty \frac{e^{-t}}{t}\,dt$
  (proven: related to the incomplete Gamma function by $E_1(x)=\Gamma(0,x)$).

---

## Available Functions

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

## Usage Examples

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

# ζ(2) = π²/6
print(f"ζ(2) = {zeta(2.0):.10f}  expected {math.pi**2/6:.10f}")
print(f"ζ(4) = {zeta(4.0):.10f}  expected {math.pi**4/90:.10f}")
```

### Spherical harmonics normalization

```python
from pysicrs import spherical_harmonic

# Y_00 = 1/sqrt(4π)
y00 = spherical_harmonic(0, 0, 0.0, 0.0)
print(f"Y_00 = {y00:.6f}  expected {1/math.sqrt(4*math.pi):.6f}")
```

---

## Numerical Accuracy

| Function | Domain | Typical |Error| · Reference |
|----------|--------|---------------------|
| $\Gamma(x)$ | $(0, 100)$ | $< 1\times10^{-12}$ | Sample values |
| $\operatorname{erf}$ | $\mathbb{R}$ | $< 1\times10^{-14}$ | A&S 7.1.26 |
| $J_0$ | $[0, 30]$ | $< 1\times10^{-12}$ | Sample values |
| $\zeta$ | $(1, 30]$ | $< 1\times10^{-10}$ | $\pi^2/6, \pi^4/90$ |

---

## Advantages & Limitations

### Advantages

✅ Covers the special functions most common in physics problems

✅ All implementations are standalone (no BLAS/LAPACK dependency)

✅ Log-space Gamma avoids overflow for large arguments

✅ Unit-tested against reference values (Abramowitz & Stegun, closed forms)

### Limitations

❌ `f64` only — no arbitrary precision

❌ Bessel functions limited to orders $n=0,1$ (general order uses recurrence/ascending series)

❌ Real axis only for most functions (no complex continuation exposed yet)

---

## References (proven standard texts)

1. Abramowitz, M. & Stegun, I. (1964). *Handbook of Mathematical Functions*. NBS.
2. Arfken, G. (1985). *Mathematical Methods for Physicists*, 3rd ed. Academic Press.
3. Press, W. et al. (2007). *Numerical Recipes*, 3rd ed. Cambridge University Press.
4. Whittaker, E.T. & Watson, G.N. (1927). *A Course of Modern Analysis*, 4th ed.

---

## Related Topics

- [Linear Algebra](linalg.md) – decompositions used by the interpolation kernels
- [Calculus](calculus.md) – quadrature rules that pair with special functions
- [Quantum](quantum.md) – spherical harmonics enter atomic wavefunctions