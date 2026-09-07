# Numerical Calculus

Pysic-rs provides the numerical differentiation, quadrature, and interpolation routines that
underpin every other module — from gradient-based ODE/PDE work to the finite differences used
in general relativity and gauge curvature.

> **Python binding status:** the Numerical Calculus suite is Rust-exported; no Python
> bindings have landed yet. See the [API page](../api/calculus.md) — examples below show the
> Rust functions.

---

## Mathematical Foundations

### Numerical Differentiation

**Central differences** (proven $O(h^2)$ accurate): for a smooth $f$,

$$
\partial_i f(x) = \frac{f(x + h e_i) - f(x - h e_i)}{2h} + O(h^2)
$$

The **Hessian** uses the symmetric second-difference stencil:

$$
\frac{\partial^2 f}{\partial x_i^2} = \frac{f(x+h e_i) - 2f(x) + f(x-he_i)}{h^2} + O(h^2)
$$

$$
\frac{\partial^2 f}{\partial x_i \partial x_j} = \frac{f(x+h(e_i+e_j)) - f(x+h e_i - h e_j) - f(x - h e_i + h e_j) + f(x - h(e_i+e_j))}{4h^2}
$$

**Recommended step size** for double precision: $h \approx \epsilon^{1/3} \sim 10^{-5}$ for
first derivatives, $h \approx \epsilon^{1/4} \sim 10^{-4}$ for second derivatives (round-off
vs. truncation balance).

### Quadrature

| Rule | Formula | Order |
|------|---------|-------|
| Trapezoid | $\int_a^b f \approx \frac{h}{2}\left(f_0 + 2\sum f_i + f_n\right)$ | $O(h^2)$ |
| Simpson | $\int_a^b f \approx \frac{h}{3}\left(f_0 + 4\sum_{\text{odd}} f_i + 2\sum_{\text{even}} f_i + f_n\right)$ | $O(h^4)$ |
| Gauss–Legendre | $\int_{-1}^1 f \approx \sum_i w_i f(x_i)$ ($n$-point) | exact for deg $\le 2n-1$ |
| Romberg | Richardson-accelerated trapezoid tableau | up to $O(h^{2k})$ |

**Gauss–Legendre** nodes are the roots of the Legendre polynomial $P_n$ (proven: quadrature is
exact for polynomials of degree $\le 2n-1$); the weights follow from the Christoffel numbers.

**Romberg** integration builds the trapezoid table $R_{j,k}$ and recursively applies
Richardson extrapolation — provably eliminating error terms one order at a time.

### Interpolation

**Linear interpolation** between $(x_i, y_i)$ and $(x_{i+1}, y_{i+1})$:

$$
y = y_i + \frac{y_{i+1}-y_i}{x_{i+1}-x_i}(x - x_i)
$$

**Cubic splines** solve the tridiagonal system for second-derivative continuity — the proven
minimal-curvature interpolant: among all $C^2$ interpolants, the natural cubic spline
minimizes $\int (f''')^2$.

---

## Routines

| Routine | Description |
|---------|-------------|
| `gradient(f, x, h)` | $\nabla f$ by central differences |
| `hessian(f, x, h)` | $\partial^2 f/\partial x_i \partial x_j$ |
| `jacobian(f, x, h)` | $J_{ij} = \partial f_i/\partial x_j$ |
| `trapz(y, x)` / `simpson(y, x)` | composite quadrature |
| `gauss_legendre(f, a, b, n)` | high-order Gaussian quadrature |
| `romberg(f, a, b, ...)` | Richardson-accelerated trapezoid |
| `lerp(x0, x1, t)` | scalar interpolation |
| `interp1d(x, y, x_query)` | piecewise linear |
| `interp_cspline(x, y, x_query)` | natural cubic spline |

---

## Usage Examples

### Gradient of a quadratic form

```python
from pysicrs import gradient

f = lambda x: x[0]**2 + 3*x[1]**2
# Exact gradient at (1,2) is (2, 12)
print(gradient(f, [1.0, 2.0], 1e-5))  # ≈ [2.0, 12.0]
```

### Gauss–Legendre for a well-behaved integral

```python
from pysicrs import gauss_legendre
import math

# ∫₀^π sin(x) dx = 2
I = gauss_legendre(lambda x: math.sin(x), 0.0, math.pi, 32)
print(f"∫ sin = {I:.12f}  (exact 2.0)")
```

### Natural cubic spline smoothing

```python
from pysicrs import interp_cspline

x = [0.0, 1.0, 2.0, 3.0]
y = [0.0, 1.0, 0.0, 1.0]
xq = 1.5
print(interp_cspline(x, y, xq))  # ≈ 0.5 for a symmetric spline
```

---

## Numerical Notes

- Central differences are used exclusively (forward differences have worse error constants
  and introduce a directional bias).
- Simpson's rule requires an even number of subintervals.
- Choose `h` per the round-off/truncation balance discussed above; the defaults in the code
  (~$10^{-5}$) are a good starting point for first derivatives.

---

## Advantages & Limitations

✅ Well-understood $O(h^k)$ convergence for smooth integrands

✅ Romberg gives near-machine-precision for analytic functions

✅ Splines produce smooth, differentiable interpolants

❌ High-order differentiation amplifies noise — not for noisy data

❌ Grid-based quadrature is not adaptive (no built-in error control yet)

---

## References

1. Press, W. et al. (2007). *Numerical Recipes*, 3rd ed. Cambs. — use `dea()` Richardson
   extrapolation and Gaussian quadrature chapters.
2. Stoer, J. & Bulirsch, R. (2002). *Introduction to Numerical Analysis*, 3rd ed. Springer.
3. Atkinson, K. (1989). *An Introduction to Numerical Analysis*, 2nd ed. Wiley.

---

## Related Topics

- [ODE](ode.md) – integrators rely on these quadrature kernels
- [General Relativity](general_relativity.md) – finite differences for Christoffel/Riemann
- [Linalg](linalg.md) – spline system solved with tridiagonal (non-pivoted LU)