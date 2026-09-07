# Linear Algebra

Pysic-rs implements the linear-algebra workhorse routines physics needs — direct factorizations,
solvers, inverses, tensor operations, and differential-geometric helpers — using `ndarray` under
the hood, with no external BLAS dependency for full reproducibility.

> **Python binding status:** `inertia_tensor` is bound from Python; the rest of the suite is
> Rust-exported. See the [API page](../api/linalg.md) for the authoritative surface.
> Examples below show the Rust functions.

---

## Mathematical Foundations

### Cholesky Decomposition

Every symmetric positive-definite matrix $A \in \mathbb{R}^{n\times n}$ admits a unique lower
triangular $L$ with positive diagonal such that:

$$
A = L L^{\mathsf T}
$$

**Algorithm** (in-place, proven stable for SPD matrices):

$$
L_{ii} = \sqrt{A_{ii} - \sum_{k<i} L_{ik}^2}, \qquad
L_{ji} = \frac{1}{L_{ii}}\left(A_{ji} - \sum_{k<i} L_{ik}L_{jk}\right)
$$

**Cost:** $\tfrac13 n^3$ flops — half of LU's. Used directly for solving $Ax=b$ via forward
and backward substitution (`solve_positive_definite`).

### LU Decomposition with Partial Pivoting

For a general square matrix $A$, LU with row pivoting produces:

$$
P A = L U
$$

where $L$ is unit lower triangular, $U$ upper triangular, $P$ a permutation matrix. Pivoting
is *required* for stability of non-SPD systems (proven: Gaussian elimination without pivoting
can be arbitrarily unstable).

### Matrix Inverse (Gauss–Jordan)

The solve-based inverse $A^{-1}$ is computed column by column using the LU factorization,
which is numerically equivalent to (and more efficient than) the classic Gauss–Jordan
elimination.

### Determinant

$$
\det A = \operatorname{sign}(P)\prod_i u_{ii}
$$

from the LU factors. Row-swap sign handling is what makes this correct — a missing sign flip
is the classic bug this implementation avoids.

### Tensor Operations & the Metric

Given the metric tensor $g_{\mu\nu}$ and its inverse $g^{\mu\nu}$ (as returned by the
`general_relativity` module), indices are raised and lowered with:

$$
V^\mu = g^{\mu\nu}V_\nu, \qquad V_\mu = g_{\mu\nu}V^\nu
$$

**Metric signature** is computed per Sylvester's law of inertia: count the eigenvalues of $g$
(normalized per row) and return the $(p,q)$ pair — the sign pattern is invariant under
congruence transformations.

### Lie Bracket (Vector Fields)

The Lie bracket of two vector fields $X, Y$ on $\mathbb{R}^n$:

$$
[X,Y] = X^j\partial_j Y - Y^j\partial_j X
$$

**Proven properties:** antisymmetry $[X,Y]=-[Y,X]$, bilinearity, and the Jacobi identity
$[X,[Y,Z]]+[Y,[Z,X]]+[Z,[X,Y]]=0$. The derivatives are computed with central finite
differences (see [Calculus](calculus.md)).

---

## Routines

| Routine | Description |
|---------|-------------|
| `cholesky(A)` | SPD → lower triangular $L$, $A=LL^\top$ |
| `solve_positive_definite(A, b)` | $Ax=b$ via Cholesky (fast, stable for SPD) |
| `lu_decompose(A)` | $PA=LU$ with partial pivoting → $(L,U,P)$ |
| `inverse(A)` | $A^{-1}$ via LU solve |
| `determinant(A)` | $\det A$ from LU with sign handling |
| `trace(A)`, `norm(v)` | matrix trace, vector norm |
| `cross(v,w)`, `dot(v,w)`, `outer(v,w)` | vector products |
| `tensor_raise_lower(g, g_inv, T, ...)` | index raising/lowering |
| `metric_signature(g)` | $(p,q)$ signature per Sylvester |
| `inverse_3x3(A)`, `det_3x3(A)` | closed-form small-matrix helpers |
| `lie_bracket(X, Y, x, h)` | $[X,Y]$ at a point |

---

## Usage Examples

### Solve a symmetric positive-definite system

```python
from pysicrs import cholesky, solve_positive_definite

A = [[4.0, 2.0],
     [2.0, 3.0]]
b = [6.0, 5.0]

x = solve_positive_definite(A, b)
print(x)          # [1.0, 1.0]

L = cholesky(A)
print(L)          # lower triangular, A = L L^T
```

### Metric signature of Minkowski space

```python
from pysicrs import metric_signature, minkowski_metric

g = minkowski_metric()
print(metric_signature(g))  # (1, 3) -> time-like, 3 space-like
```

### Lie bracket of coordinate fields

```python
from pysicrs import lie_bracket

X = lambda p: [1.0, 0.0, 0.0]   # ∂/∂x
Y = lambda p: [0.0, 1.0, 0.0]   # ∂/∂y
print(lie_bracket(X, Y, [0.0]*3, 1e-5))  # [0, 0, 0] — commute
```

---

## Numerical Notes

- All routines are pure Rust on top of `ndarray` — no optional BLAS, so results are
  deterministic across machines.
- Cholesky requires SPD input; a non-positive pivot raises `PysicError::SingularMatrix`.
- The `3x3` helpers are closed-form (Cramer / adjugate) — branch-free and fast for the
  many 3-vectors in classical mechanics.

---

## Advantages & Limitations

✅ Self-contained, deterministic, no BLAS/LAPACK installation needed

✅ Correct pivot sign handling for determinants

✅ Index raising/lowering integrates directly with GR metrics

❌ Dense solvers only — no sparse linear algebra in this module (see Optimiz-rs for ADMM/Lasso)

❌ `f64` precision; no automatic differentiation aware solvers

---

## References

1. Golub, G.H. & Van Loan, C.F. (2013). *Matrix Computations*, 4th ed. JHU Press.
2. Trefethen, L.N. & Bau, D. (1997). *Numerical Linear Algebra*. SIAM.
3. Wald, R. (1984). *General Relativity*. Chicago. (tensor index conventions)

---

## Related Topics

- [General Relativity](general_relativity.md) – sources metric tensors and inverse metrics
- [Calculus](calculus.md) – supplies the finite differences used by `lie_bracket`
- [Special Functions](special_functions.md) – analytic kernels used in solving systems