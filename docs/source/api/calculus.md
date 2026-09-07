# API Reference: Numerical Calculus

The Python surface exposes numerical-calculus routines at top level
(`from pysicrs import <name>`). See also
[Algorithms: Numerical Calculus](../algorithms/calculus.md).

```python
import pysicrs
```

## Rust API (not yet Python-bound)

The following Rust-exported functions are described on the
[Calculus algorithms page](../algorithms/calculus.md), and are exposed to Python as soon
as their module bindings land:

| Rust function | purpose |
|---------------|---------|
| `gradient(f, x, h)` | $\nabla f$ central-difference |
| `hessian(f, x, h)` | second-derivative matrix |
| `jacobian(f, x, h)` | $J_{ij}=\partial f_i/\partial x_j$ |
| `trapz(y, x)` / `simpson(y, x)` | composite quadrature |
| `gauss_legendre(f, a, b, n)` | Gaussian quadrature |
| `romberg(f, a, b, …)` | Richardson-accelerated |
| `lerp`, `interp1d`, `interp_cspline` | linear / cubic spline |