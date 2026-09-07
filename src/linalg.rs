//! Linear algebra utilities with tensor awareness for physics.

use ndarray::{Array1, Array2};
use crate::core::{PysicError, Result};

// ============================================================================
// Decompositions
// ============================================================================

/// Cholesky decomposition A = LL^T (lower triangular).
pub fn cholesky(a: &Array2<f64>) -> Result<Array2<f64>> {
    let n = a.nrows();
    if n != a.ncols() {
        return Err(PysicError::NonSquareMatrix { rows: n, cols: a.ncols() });
    }
    let mut l = Array2::zeros((n, n));
    for i in 0..n {
        for j in 0..=i {
            let mut sum = 0.0;
            for k in 0..j {
                sum += l[[i, k]] * l[[j, k]];
            }
            let val = a[[i, j]] - sum;
            if i == j {
                if val <= 0.0 {
                    return Err(PysicError::SingularMatrix);
                }
                l[[i, j]] = val.sqrt();
            } else {
                l[[i, j]] = val / l[[j, j]];
            }
        }
    }
    Ok(l)
}

/// Solve Ax = b for symmetric positive-definite A using Cholesky.
pub fn solve_positive_definite(a: &Array2<f64>, b: &Array1<f64>) -> Result<Array1<f64>> {
    let l = cholesky(a)?;
    let n = b.len();
    let mut y = Array1::zeros(n);
    for i in 0..n {
        let mut sum = 0.0;
        for j in 0..i { sum += l[[i, j]] * y[j]; }
        y[i] = (b[i] - sum) / l[[i, i]];
    }
    let mut x = Array1::zeros(n);
    for i in (0..n).rev() {
        let mut sum = 0.0;
        for j in i + 1..n { sum += l[[j, i]] * x[j]; }
        x[i] = (y[i] - sum) / l[[i, i]];
    }
    Ok(x)
}

/// LU decomposition with partial pivoting. Returns (L, U, P).
pub fn lu_decompose(a: &Array2<f64>) -> Result<(Array2<f64>, Array2<f64>, Array2<f64>)> {
    let n = a.nrows();
    if n != a.ncols() {
        return Err(PysicError::NonSquareMatrix { rows: n, cols: a.ncols() });
    }
    let mut lu = a.clone();
    let mut piv: Vec<usize> = (0..n).collect();

    for k in 0..n {
        // Find pivot
        let mut max_val = 0.0;
        let mut max_idx = k;
        for i in k..n {
            if lu[[i, k]].abs() > max_val {
                max_val = lu[[i, k]].abs();
                max_idx = i;
            }
        }
        if max_val < 1e-15 {
            return Err(PysicError::SingularMatrix);
        }
        if max_idx != k {
            piv.swap(k, max_idx);
            for j in 0..n {
                let tmp = lu[[k, j]];
                lu[[k, j]] = lu[[max_idx, j]];
                lu[[max_idx, j]] = tmp;
            }
        }
        for i in k + 1..n {
            lu[[i, k]] /= lu[[k, k]];
            for j in k + 1..n {
                lu[[i, j]] -= lu[[i, k]] * lu[[k, j]];
            }
        }
    }

    let mut l = Array2::eye(n);
    let mut u = Array2::zeros((n, n));
    let mut p = Array2::zeros((n, n));
    for i in 0..n {
        p[[i, piv[i]]] = 1.0;
        for j in 0..i { l[[i, j]] = lu[[i, j]]; }
        for j in i..n { u[[i, j]] = lu[[i, j]]; }
    }

    Ok((l, u, p))
}

/// Matrix inverse via Gauss-Jordan elimination.
pub fn inverse(a: &Array2<f64>) -> Result<Array2<f64>> {
    let n = a.nrows();
    if n != a.ncols() {
        return Err(PysicError::NonSquareMatrix { rows: n, cols: a.ncols() });
    }
    let mut aug = Array2::<f64>::zeros((n, 2 * n));
    for i in 0..n {
        for j in 0..n {
            aug[[i, j]] = a[[i, j]];
        }
        aug[[i, n + i]] = 1.0;
    }

    for col in 0..n {
        let mut max_row = col;
        for row in col..n {
            if aug[[row, col]].abs() > aug[[max_row, col]].abs() {
                max_row = row;
            }
        }
        if max_row != col {
            for j in 0..2 * n {
                let tmp = aug[[col, j]];
                aug[[col, j]] = aug[[max_row, j]];
                aug[[max_row, j]] = tmp;
            }
        }

        let pivot = aug[[col, col]];
        if pivot.abs() < 1e-15 {
            return Err(PysicError::SingularMatrix);
        }

        for j in 0..2 * n {
            aug[[col, j]] /= pivot;
        }

        for row in 0..n {
            if row == col { continue; }
            let factor = aug[[row, col]];
            for j in 0..2 * n {
                aug[[row, j]] -= factor * aug[[col, j]];
            }
        }
    }

    let mut inv = Array2::zeros((n, n));
    for i in 0..n {
        for j in 0..n {
            inv[[i, j]] = aug[[i, n + j]];
        }
    }
    Ok(inv)
}

// ============================================================================
// Matrix properties
// ============================================================================

/// Trace of a square matrix.
pub fn trace(m: &Array2<f64>) -> Result<f64> {
    let n = m.nrows();
    if n != m.ncols() {
        return Err(PysicError::NonSquareMatrix { rows: n, cols: m.ncols() });
    }
    Ok((0..n).map(|i| m[[i, i]]).sum())
}

/// Determinant via LU.
pub fn determinant(a: &Array2<f64>) -> Result<f64> {
    let n = a.nrows();
    if n != a.ncols() {
        return Err(PysicError::NonSquareMatrix { rows: n, cols: a.ncols() });
    }
    let mut det = 1.0;
    let mut lu = a.clone();

    for k in 0..n {
        let mut max_val = 0.0;
        let mut max_idx = k;
        for i in k..n {
            if lu[[i, k]].abs() > max_val {
                max_val = lu[[i, k]].abs();
                max_idx = i;
            }
        }
        if max_val < 1e-15 {
            return Ok(0.0);
        }
        if max_idx != k {
            for j in 0..n {
                let tmp = lu[[k, j]];
                lu[[k, j]] = lu[[max_idx, j]];
                lu[[max_idx, j]] = tmp;
            }
            det = -det;
        }
        det *= lu[[k, k]];
        for i in k + 1..n {
            lu[[i, k]] /= lu[[k, k]];
            for j in k + 1..n {
                lu[[i, j]] -= lu[[i, k]] * lu[[k, j]];
            }
        }
    }
    Ok(det)
}

/// Frobenius norm.
pub fn matrix_norm_frobenius(m: &Array2<f64>) -> f64 {
    m.iter().map(|&x| x * x).sum::<f64>().sqrt()
}

/// Vector L2 norm.
pub fn vector_norm(v: &Array1<f64>) -> f64 {
    v.iter().map(|&x| x * x).sum::<f64>().sqrt()
}

/// Normalize vector to unit length.
pub fn normalize(v: &Array1<f64>) -> Result<Array1<f64>> {
    let n = vector_norm(v);
    if n < 1e-15 {
        return Err(PysicError::InvalidInput("Cannot normalize zero vector".into()));
    }
    Ok(v / n)
}

/// Cross product of two 3-vectors.
pub fn cross_product_3d(a: &Array1<f64>, b: &Array1<f64>) -> Result<Array1<f64>> {
    if a.len() != 3 || b.len() != 3 {
        return Err(PysicError::InvalidInput("Cross product requires 3-vectors".into()));
    }
    Ok(Array1::from_vec(vec![
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]))
}

/// Dot product of two vectors.
pub fn dot_product(a: &Array1<f64>, b: &Array1<f64>) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Outer product of two vectors: A_ij = a_i * b_j.
pub fn outer_product(a: &Array1<f64>, b: &Array1<f64>) -> Array2<f64> {
    let mut result = Array2::zeros((a.len(), b.len()));
    for i in 0..a.len() {
        for j in 0..b.len() {
            result[[i, j]] = a[i] * b[j];
        }
    }
    result
}

// ============================================================================
// Tensor operations for physics
// ============================================================================

/// Raise an index using the inverse metric: v^i = g^{ij} v_j.
pub fn tensor_raise_lower(metric_inv: &Array2<f64>, covector: &Array1<f64>) -> Result<Array1<f64>> {
    let n = metric_inv.nrows();
    if n != covector.len() {
        return Err(PysicError::DimensionMismatch { expected: n, actual: covector.len() });
    }
    let mut result = Array1::zeros(n);
    for i in 0..n {
        for j in 0..n {
            result[i] += metric_inv[[i, j]] * covector[j];
        }
    }
    Ok(result)
}

/// Determine metric signature from a diagonal matrix.
/// Returns (num_positive, num_negative, num_zero).
pub fn metric_signature(g: &Array2<f64>) -> (usize, usize, usize) {
    let n = g.nrows();
    let mut pos = 0;
    let mut neg = 0;
    let mut zero = 0;
    for i in 0..n {
        let d = g[[i, i]];
        if d > 0.0 { pos += 1; }
        else if d < 0.0 { neg += 1; }
        else { zero += 1; }
    }
    (pos, neg, zero)
}

/// 3x3 matrix inverse (direct formula).
pub fn inverse_3x3(m: &Array2<f64>) -> Result<Array2<f64>> {
    if m.nrows() != 3 || m.ncols() != 3 {
        return Err(PysicError::NonSquareMatrix { rows: m.nrows(), cols: m.ncols() });
    }
    let det = m[[0,0]] * (m[[1,1]] * m[[2,2]] - m[[1,2]] * m[[2,1]])
            - m[[0,1]] * (m[[1,0]] * m[[2,2]] - m[[1,2]] * m[[2,0]])
            + m[[0,2]] * (m[[1,0]] * m[[2,1]] - m[[1,1]] * m[[2,0]]);

    if det.abs() < 1e-15 {
        return Err(PysicError::SingularMatrix);
    }

    let inv_det = 1.0 / det;
    Ok(Array2::from_shape_fn((3, 3), |(i, j)| {
        let (r0, r1) = match i { 0 => (1,2), 1 => (0,2), _ => (0,1) };
        let (c0, c1) = match j { 0 => (1,2), 1 => (0,2), _ => (0,1) };
        let sign = if (i + j) % 2 == 0 { 1.0 } else { -1.0 };
        sign * inv_det * (m[[r0,c0]] * m[[r1,c1]] - m[[r0,c1]] * m[[r1,c0]])
    }))
}

/// Determinant of a 3x3 matrix.
pub fn det_3x3(m: &Array2<f64>) -> f64 {
    m[[0,0]] * (m[[1,1]] * m[[2,2]] - m[[1,2]] * m[[2,1]])
  - m[[0,1]] * (m[[1,0]] * m[[2,2]] - m[[1,2]] * m[[2,0]])
  + m[[0,2]] * (m[[1,0]] * m[[2,1]] - m[[1,1]] * m[[2,0]])
}

/// Lie bracket [X, Y] for vector fields (discrete: finite differences).
pub fn lie_bracket(
    x_field: &dyn Fn(&[f64]) -> Vec<f64>,
    y_field: &dyn Fn(&[f64]) -> Vec<f64>,
    point: &[f64],
    h: f64,
) -> Vec<f64> {
    let n = point.len();
    let mut result = vec![0.0; n];

    // [X, Y]^i = X^j ∂_j Y^i - Y^j ∂_j X^i
    for k in 0..n {
        // Numerical derivatives of Y along X
        let mut d_y_dx = vec![0.0; n];
        let mut d_x_dy = vec![0.0; n];

        for j in 0..n {
            // ∂_j Y^i
            let mut x_plus = point.to_vec();
            let mut x_minus = point.to_vec();
            x_plus[j] += h;
            x_minus[j] -= h;
            let y_plus = y_field(&x_plus);
            let y_minus = y_field(&x_minus);
            let x_plus2 = x_field(&x_plus);
            let x_minus2 = x_field(&x_minus);

            for i in 0..n {
                d_y_dx[i] += x_field(point)[j] * (y_plus[i] - y_minus[i]) / (2.0 * h);
                d_x_dy[i] += y_field(point)[j] * (x_plus2[i] - x_minus2[i]) / (2.0 * h);
            }
        }

        result[k] = d_y_dx[k] - d_x_dy[k];
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn test_cholesky() {
        let a = array![[4.0, 2.0], [2.0, 3.0]];
        let l = cholesky(&a).unwrap();
        // Verify L * L^T ≈ A
        let lt = l.t();
        let recon = l.dot(&lt);
        assert!((recon[[0,0]] - 4.0).abs() < 1e-10);
        assert!((recon[[0,1]] - 2.0).abs() < 1e-10);
        assert!((recon[[1,0]] - 2.0).abs() < 1e-10);
        assert!((recon[[1,1]] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_determinant() {
        let a = array![[1.0, 2.0], [3.0, 4.0]];
        assert!((determinant(&a).unwrap() - (-2.0)).abs() < 1e-10);
    }

    #[test]
    fn test_inverse() {
        let a = array![[1.0, 2.0], [3.0, 4.0]];
        let inv = inverse(&a).unwrap();
        let product = a.dot(&inv);
        assert!((product[[0,0]] - 1.0).abs() < 1e-10);
        assert!((product[[0,1]] - 0.0).abs() < 1e-10);
        assert!((product[[1,0]] - 0.0).abs() < 1e-10);
        assert!((product[[1,1]] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_cross_product() {
        let a = array![1.0, 0.0, 0.0];
        let b = array![0.0, 1.0, 0.0];
        let c = cross_product_3d(&a, &b).unwrap();
        assert!((c[0] - 0.0).abs() < 1e-10);
        assert!((c[1] - 0.0).abs() < 1e-10);
        assert!((c[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_metric_signature() {
        let mut diag = Array2::zeros((4, 4));
        diag[[0,0]] = -1.0;
        diag[[1,1]] = 1.0;
        diag[[2,2]] = 1.0;
        diag[[3,3]] = 1.0;
        let sig = metric_signature(&diag);
        assert_eq!(sig, (3, 1, 0)); // (+,+,+,-) = Lorentzian
    }
}
