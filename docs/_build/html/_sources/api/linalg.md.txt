# API Reference: Linear Algebra

The Python surface exposes linear-algebra routines at top level
(`from pysicrs import <name>`). Bindings are added module by module; see also
[Algorithms: Linear Algebra](../algorithms/linalg.md).

```python
import pysicrs
```

## Currently bound

```python
inertia_tensor(masses, positions) -> list[list[float]]
```

Moment-of-inertia tensor $I_{ij}=\sum_a m_a(\|r_a\|^2\delta_{ij}-r_{a,i}r_{a,j})$
for point masses. `masses` is a list of floats; `positions` a list of 3-vectors.

## Rust API (not yet Python-bound)

The following Rust-exported functions are described on the
[Linear Algebra algorithms page](../algorithms/linalg.md):

| Rust function | purpose |
|---------------|---------|
| `cholesky(A)` | $A=LL^\top$, lower $L$ |
| `solve_positive_definite(A, b)` | Cholesky-based solve |
| `lu_decompose(A)` | $PA=LU$ with pivoting |
| `inverse(A)` | Gauss–Jordan inverse |
| `determinant(A)` | LU determinant with sign |
| `trace(A)`, `norm(v)`, `cross`, `dot`, `outer` | vector/operator basics |
| `tensor_raise_lower(g, g_inv, T, …)` | index gymnastics |
| `metric_signature(g)` | $(p,q)$ signature |
| `inverse_3x3`, `det_3x3` | closed-form small matrices |
| `lie_bracket(X, Y, x, h)` | $[X,Y]$ |