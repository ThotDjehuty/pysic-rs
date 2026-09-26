# Numerical Calculus

Pysic-rs provides the numerical differentiation, quadrature, and interpolation routines that
underpin every other module — from gradient-based ODE/PDE work to the finite differences used
in general relativity and gauge curvature.

> **Python binding status:** the Numerical Calculus suite is Rust-exported; no Python
> bindings have landed yet. See the [API page](../api/calculus.md).

---

## 1. Numerical Differentiation

### 1.1 Central Differences

For a smooth function $f:\mathbb{R}^n\to\mathbb{R}$, the **central difference** approximation
to the gradient is

$$
\boxed{\partial_i f(\mathbf{x}) = \frac{f(\mathbf{x}+h\mathbf{e}_i)-f(\mathbf{x}-h\mathbf{e}_i)}{2h} + \mathcal{O}(h^2)}
$$

*Proof.* Taylor expand: $f(\mathbf{x}+h\mathbf{e}_i) = f + h\partial_i f + \frac{h^2}{2}
\partial_i^2 f + \frac{h^3}{6}\partial_i^3 f + \cdots$, and similarly for $f(\mathbf{x}
-h\mathbf{e}_i)$. Subtract: the even terms cancel, the leading error is $\frac{h^2}{6}
\partial_i^3 f$. $\square$

Central differences are always preferred over forward/backward: the error constant is half,
and there is no directional bias.

### 1.2 Hessian (Second Derivatives)

**Diagonal:**

$$
\frac{\partial^2 f}{\partial x_i^2} = \frac{f(\mathbf{x}+h\mathbf{e}_i)-2f(\mathbf{x})+f(\mathbf{x}-h\mathbf{e}_i)}{h^2} + \mathcal{O}(h^2)
$$

**Off-diagonal:**

$$
\frac{\partial^2 f}{\partial x_i\partial x_j} = \frac{f(\mathbf{x}+h(\mathbf{e}_i+\mathbf{e}_j))-f(\mathbf{x}+h\mathbf{e}_i-h\mathbf{e}_j)-f(\mathbf{x}-h\mathbf{e}_i+h\mathbf{e}_j)+f(\mathbf{x}-h(\mathbf{e}_i+\mathbf{e}_j))}{4h^2}
$$

### 1.3 Optimal Step Size

For double precision ($\epsilon_{\rm mach}\approx 10^{-16}$):

- First derivative: $h\approx\epsilon_{\rm mach}^{1/3}\sim 10^{-5}$ (truncation $\sim h^2$,
  round-off $\sim \epsilon/h$; balance at $h\sim\epsilon^{1/3}$).
- Second derivative: $h\approx\epsilon_{\rm mach}^{1/4}\sim 10^{-4}$.

### 1.4 Jacobian

The Jacobian of $F:\mathbb{R}^n\to\mathbb{R}^m$ is computed column by column:

$$
J_{ij} = \frac{\partial F_i}{\partial x_j} = \frac{F_i(\mathbf{x}+h\mathbf{e}_j)-F_i(\mathbf{x}-h\mathbf{e}_j)}{2h}
$$

Cost: $2n$ function evaluations (vs $n^2$ for finite-difference Hessian).

---

## 2. Quadrature (Numerical Integration)

### 2.1 Composite Trapezoid Rule

$$
\int_a^b f(x)\,dx \approx \frac{h}{2}\left(f_0 + 2\sum_{i=1}^{N-1}f_i + f_N\right), \qquad h = \frac{b-a}{N}
$$

**Error:** $\mathcal{O}(h^2)$ (proven from the Euler–Maclaurin formula).

### 2.2 Simpson's Rule

$$
\int_a^b f(x)\,dx \approx \frac{h}{3}\left(f_0 + 4\sum_{\text{odd}}f_i + 2\sum_{\text{even}}f_i + f_N\right)
$$

**Error:** $\mathcal{O}(h^4)$ (proven — exact for polynomials of degree $\le 3$).
Requires an **even** number of subintervals.

### 2.3 Gauss–Legendre Quadrature

$$
\int_{-1}^{1}f(x)\,dx \approx \sum_{i=1}^n w_i f(x_i)
$$

where $x_i$ are the roots of the Legendre polynomial $P_n(x)$ and $w_i$ are the
Christoffel numbers.

**Theorem (exactness).** Gauss–Legendre with $n$ points is exact for polynomials of degree
$\le 2n-1$.

*Proof.* The $n$-point rule has $2n$ free parameters ($n$ nodes + $n$ weights). Matching
the moments $\int x^k\,dx$ for $k=0,\ldots,2n-1$ gives $2n$ equations in $2n$ unknowns,
which has a unique solution. The nodes must be the roots of $P_n$ (by the orthogonality
of Legendre polynomials). $\square$

### 2.4 Romberg Integration

Build the trapezoid table $R_{j,k}$ and apply Richardson extrapolation:

$$
R_{j,k} = R_{j,k-1} + \frac{R_{j,k-1}-R_{j-1,k-1}}{4^k-1}
$$

This eliminates error terms one order at a time, provably reaching $\mathcal{O}(h^{2k})$
accuracy from the basic trapezoid rule.

---

## 3. Interpolation

### 3.1 Linear Interpolation

Between two points $(x_i,y_i)$ and $(x_{i+1},y_{i+1})$:

$$
y = y_i + \frac{y_{i+1}-y_i}{x_{i+1}-x_i}(x-x_i)
$$

Error: $\mathcal{O}(h^2)$ for smooth $f$ (from the remainder of the Taylor expansion).

### 3.2 Natural Cubic Spline

The **natural cubic spline** is the unique $C^2$ interpolant that minimises $\int_a^b
(f''')^2\,dx$ (the total curvature). It solves a tridiagonal system for the second
derivatives $M_i = f''(x_i)$:

$$
h_{i-1}M_{i-1} + 2(h_{i-1}+h_i)M_i + h_i M_{i+1} = 6\left(\frac{y_{i+1}-y_i}{h_i}-\frac{y_i-y_{i-1}}{h_{i-1}}\right)
$$

with natural boundary conditions $M_0 = M_N = 0$.

**Properties (proven):**
- Second-order accurate: $\|f-S\|_\infty = \mathcal{O}(h^2)$ (for $f\in C^4$).
- The minimal-curvature property makes it the smoothest interpolant.
- Tridiagonal solve: $O(N)$ operations.

---

## 4. Routines

| Routine | Description |
|---------|-------------|
| `gradient(f, x, h)` | $\nabla f$ by central differences |
| `hessian(f, x, h)` | $\partial^2 f/\partial x_i\partial x_j$ |
| `jacobian(f, x, h)` | $J_{ij}=\partial f_i/\partial x_j$ |
| `trapz(y, x)` / `simpson(y, x)` | composite quadrature |
| `gauss_legendre(f, a, b, n)` | high-order Gaussian quadrature |
| `romberg(f, a, b, ...)` | Richardson-accelerated trapezoid |
| `lerp(x0, x1, t)` | scalar interpolation |
| `interp1d(x, y, x_query)` | piecewise linear |
| `interp_cspline(x, y, x_query)` | natural cubic spline |

---

## 5. Usage Examples

### Gradient of a quadratic form

```python
from pysicrs import gradient

f = lambda x: x[0]**2 + 3*x[1]**2
print(gradient(f, [1.0, 2.0], 1e-5))  # ≈ [2.0, 12.0]
```

### Gauss–Legendre for a well-behaved integral

```python
from pysicrs import gauss_legendre
import math

I = gauss_legendre(lambda x: math.sin(x), 0.0, math.pi, 32)
print(f"∫ sin = {I:.12f}  (exact 2.0)")
```

### Natural cubic spline smoothing

```python
from pysicrs import interp_cspline

x = [0.0, 1.0, 2.0, 3.0]
y = [0.0, 1.0, 0.0, 1.0]
print(interp_cspline(x, y, 1.5))  # ≈ 0.5
```

---

## 6. Numerical Notes

- Central differences are used exclusively (forward differences have worse error constants).
- Simpson's rule requires an even number of subintervals.
- Choose $h$ per the round-off/truncation balance; defaults ($\sim 10^{-5}$) are good for
  first derivatives.

---

## 7. Advantages & Limitations

✅ Well-understood $\mathcal{O}(h^k)$ convergence for smooth integrands

✅ Romberg gives near-machine-precision for analytic functions

✅ Splines produce smooth, differentiable interpolants

❌ High-order differentiation amplifies noise — not for noisy data

❌ Grid-based quadrature is not adaptive (no built-in error control yet)

---

## 8. References

1. Press, W. et al. (2007). *Numerical Recipes*, 3rd ed. Cambs.
2. Stoer, J. & Bulirsch, R. (2002). *Introduction to Numerical Analysis*, 3rd ed. Springer.
3. Atkinson, K. (1989). *An Introduction to Numerical Analysis*, 2nd ed. Wiley.

---

## 9. Related Topics

- [ODE](ode.md) — integrators rely on these quadrature kernels
- [General Relativity](general_relativity.md) — finite differences for Christoffel/Riemann
- [Linear Algebra](linalg.md) — spline system solved with tridiagonal (non-pivoted LU)
- [PDE](pde.md) — spatial discretisation uses these differences
