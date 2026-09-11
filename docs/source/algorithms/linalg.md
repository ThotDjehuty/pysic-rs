# Linear Algebra

Pysic-rs implements the linear-algebra workhorse routines physics needs — direct factorisations,
solvers, inverses, tensor operations, and differential-geometric helpers — using `ndarray` under
the hood, with no external BLAS dependency for full reproducibility.

> **Python binding status:** `inertia_tensor` is bound from Python; the rest of the suite is
> Rust-exported. See the [API page](../api/linalg.md) for the authoritative surface.

---

## 1. Matrix Factorisations

### 1.1 Cholesky Decomposition

**Theorem.** Every symmetric positive-definite (SPD) matrix $A\in\mathbb{R}^{n\times n}$
admits a unique lower triangular $L$ with positive diagonal entries such that

$$
\boxed{A = LL^{\mathsf{T}}}
$$

*Proof.* The existence follows from induction on $n$. Write $A = \begin{pmatrix}a_{11}&\mathbf{b}^{\mathsf{T}}\\\mathbf{b}&C\end{pmatrix}$.
Set $L_{11} = \sqrt{a_{11}}$, $\mathbf{l} = \mathbf{b}/L_{11}$. Then the $(2,2)$ block
$C - \mathbf{l}\mathbf{l}^{\mathsf{T}}$ is still SPD (by the Schur complement), so
recursion applies. Uniqueness: if $LL^{\mathsf{T}} = L'L'^{\mathsf{T}}$, then $L'^{-1}L$
is orthogonal and lower triangular, hence diagonal with unit entries, so $L=L'$. $\square$

**Algorithm (in-place):**

$$
L_{ii} = \sqrt{A_{ii}-\sum_{k<i}L_{ik}^2}, \qquad
L_{ji} = \frac{1}{L_{ii}}\left(A_{ji}-\sum_{k<i}L_{ik}L_{jk}\right)
$$

**Cost:** $\frac{1}{3}n^3$ flops — half of LU's.

### 1.2 LU Decomposition with Partial Pivoting

**Theorem.** For a general square matrix $A$, there exist a unit lower triangular $L$,
an upper triangular $U$, and a permutation matrix $P$ such that

$$
\boxed{PA = LU}
$$

*Proof.* Gaussian elimination with row pivoting (choosing the largest element in each
column as pivot) is always possible for non-singular matrices. The elimination steps
produce $L$, and the row swaps produce $P$. $\square$

**Pivoting is required** for numerical stability: Gaussian elimination without pivoting
can be arbitrarily unstable (the classic example is the Wilkinson matrix).

### 1.3 Determinant from LU

$$
\det A = \operatorname{sign}(P)\prod_{i=1}^n u_{ii}
$$

where $\operatorname{sign}(P) = (-1)^s$ with $s$ the number of row swaps. This avoids
the exponentially expensive Leibniz formula.

### 1.4 Matrix Inverse

The inverse $A^{-1}$ is computed column by column using the LU factorisation: solve
$Ax_j = e_j$ for each standard basis vector $e_j$. This is numerically equivalent to
Gauss–Jordan elimination but more efficient ($O(n^3)$ vs $O(n^3)$ with smaller constants).

---

## 2. Tensor Operations & the Metric

### 2.1 Index Raising/Lowering

Given the metric tensor $g_{\mu\nu}$ and its inverse $g^{\mu\nu}$ (as returned by the
[General Relativity](general_relativity.md) module):

$$
V^\mu = g^{\mu\nu}V_\nu, \qquad V_\mu = g_{\mu\nu}V^\nu
$$

For a rank-2 tensor: $T^{\mu}{}_{\nu} = g^{\mu\lambda}T_{\lambda\nu}$, etc.

### 2.2 Metric Signature

The **signature** of a metric $g$ is the pair $(p,q)$ counting the positive and negative
eigenvalues. By **Sylvester's law of inertia**, the signature is invariant under
congruence transformations $g\to M^{\mathsf{T}}gM$.

**Algorithm:** compute the eigenvalues of $g$ (or, more efficiently, count sign changes in
the leading principal minors — the `det_3x3` and `inverse_3x3` closed-form helpers are
used for the common 3D and 4D cases).

---

## 3. Vector Operations

### 3.1 Cross Product

In 3D:

$$
(\mathbf{v}\times\mathbf{w})_i = \varepsilon_{ijk}v_j w_k
$$

Implemented via the Levi-Civita symbol for $n=3$ only (no generalisation to other
dimensions).

### 3.2 Lie Bracket of Vector Fields

For two vector fields $X,Y$ on $\mathbb{R}^n$:

$$
[X,Y]^i = X^j\partial_j Y^i - Y^j\partial_j X^i
$$

**Properties (proven):**
- Antisymmetry: $[X,Y] = -[Y,X]$.
- Bilinearity.
- Jacobi identity: $[X,[Y,Z]] + [Y,[Z,X]] + [Z,[X,Y]] = 0$.

Derivatives are computed with central finite differences (see [Calculus](calculus.md)).

---

## 4. Routines

| Routine | Description |
|---------|-------------|
| `cholesky(A)` | SPD → lower triangular $L$, $A=LL^{\mathsf{T}}$ |
| `solve_positive_definite(A, b)` | $Ax=b$ via Cholesky |
| `lu_decompose(A)` | $PA=LU$ with partial pivoting → $(L,U,P)$ |
| `inverse(A)` | $A^{-1}$ via LU solve |
| `determinant(A)` | $\det A$ from LU with sign |
| `trace(A)`, `norm(v)` | matrix trace, vector norm |
| `cross(v,w)`, `dot(v,w)`, `outer(v,w)` | vector products |
| `tensor_raise_lower(g, g_inv, T, ...)` | index raising/lowering |
| `metric_signature(g)` | $(p,q)$ signature per Sylvester |
| `inverse_3x3(A)`, `det_3x3(A)` | closed-form small-matrix helpers |
| `lie_bracket(X, Y, x, h)` | $[X,Y]$ at a point |

---

## 5. Usage Examples

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
print(metric_signature(g))  # (1, 3) → Lorentzian
```

### Lie bracket of coordinate fields

```python
from pysicrs import lie_bracket

X = lambda p: [1.0, 0.0, 0.0]   # ∂/∂x
Y = lambda p: [0.0, 1.0, 0.0]   # ∂/∂y
print(lie_bracket(X, Y, [0.0]*3, 1e-5))  # [0, 0, 0] — commute
```

---

## 6. Numerical Notes

- All routines are pure Rust on top of `ndarray` — no optional BLAS, so results are
  deterministic across machines.
- Cholesky requires SPD input; a non-positive pivot raises `PysicError::SingularMatrix`.
- The $3\times 3$ helpers are closed-form (Cramer / adjugate) — branch-free and fast for
  the many 3-vectors in classical mechanics.

---

## 7. Advantages & Limitations

✅ Self-contained, deterministic, no BLAS/LAPACK installation needed

✅ Correct pivot sign handling for determinants

✅ Index raising/lowering integrates directly with GR metrics

❌ Dense solvers only — no sparse linear algebra (see Optimiz-rs for ADMM/Lasso)

❌ `f64` precision; no automatic-differentiation-aware solvers

---

## 8. References

1. Golub, G.H. & Van Loan, C.F. (2013). *Matrix Computations*, 4th ed. JHU Press.
2. Trefethen, L.N. & Bau, D. (1997). *Numerical Linear Algebra*. SIAM.
3. Wald, R. (1984). *General Relativity*. Chicago.

---

## 9. Related Topics

- [General Relativity](general_relativity.md) — sources metric tensors and inverse metrics
- [Calculus](calculus.md) — supplies the finite differences used by `lie_bracket`
- [Special Functions](special_functions.md) — analytic kernels used in solving systems
- [Classical](classical.md) — inertia tensors are symmetric matrices from this module
