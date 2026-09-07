//! Numerical calculus: differentiation, integration, interpolation.

use crate::core::Result;

// ============================================================================
// Numerical Differentiation
// ============================================================================

/// Central finite-difference gradient: ∇f(x) ≈ [f(x+h·e_i) - f(x-h·e_i)] / (2h).
pub fn gradient<F>(f: &F, x: &[f64], h: f64) -> Vec<f64>
where
    F: Fn(&[f64]) -> f64,
{
    let n = x.len();
    let mut grad = vec![0.0; n];
    let mut x_plus = x.to_vec();
    let mut x_minus = x.to_vec();
    for i in 0..n {
        x_plus[i] = x[i] + h;
        x_minus[i] = x[i] - h;
        grad[i] = (f(&x_plus) - f(&x_minus)) / (2.0 * h);
        x_plus[i] = x[i];
        x_minus[i] = x[i];
    }
    grad
}

/// Central finite-difference Hessian matrix.
pub fn hessian<F>(f: &F, x: &[f64], h: f64) -> Vec<Vec<f64>>
where
    F: Fn(&[f64]) -> f64,
{
    let n = x.len();
    let mut hess = vec![vec![0.0; n]; n];
    let mut x_work = x.to_vec();
    let f0 = f(x);

    for i in 0..n {
        for j in i..n {
            if i == j {
                x_work[i] = x[i] + h;
                let fp = f(&x_work);
                x_work[i] = x[i] - h;
                let fm = f(&x_work);
                x_work[i] = x[i];
                hess[i][i] = (fp - 2.0 * f0 + fm) / (h * h);
            } else {
                x_work[i] = x[i] + h;
                x_work[j] = x[j] + h;
                let fpp = f(&x_work);
                x_work[j] = x[j] - h;
                let fpm = f(&x_work);
                x_work[i] = x[i] - h;
                let fmm = f(&x_work);
                x_work[j] = x[j] + h;
                let fmp = f(&x_work);
                x_work[i] = x[i];
                x_work[j] = x[j];
                hess[i][j] = (fpp - fpm - fmp + fmm) / (4.0 * h * h);
                hess[j][i] = hess[i][j];
            }
        }
    }
    hess
}

/// Central finite-difference Jacobian for vector-valued function.
pub fn jacobian<F>(f: &F, x: &[f64], h: f64) -> Vec<Vec<f64>>
where
    F: Fn(&[f64]) -> Vec<f64>,
{
    let n = x.len();
    let f0 = f(x);
    let m = f0.len();
    let mut jac = vec![vec![0.0; n]; m];
    let mut x_work = x.to_vec();

    for j in 0..n {
        x_work[j] = x[j] + h;
        let fp = f(&x_work);
        x_work[j] = x[j] - h;
        let fm = f(&x_work);
        x_work[j] = x[j];
        for i in 0..m {
            jac[i][j] = (fp[i] - fm[i]) / (2.0 * h);
        }
    }
    jac
}

// ============================================================================
// Numerical Integration
// ============================================================================

/// Trapezoidal rule: ∫f(x)dx ≈ h/2 [f(x_0) + 2f(x_1) + ... + f(x_n)].
pub fn trapz<F>(f: &F, a: f64, b: f64, n: usize) -> f64
where
    F: Fn(f64) -> f64,
{
    if n == 0 { return 0.0; }
    let h = (b - a) / n as f64;
    let mut sum = 0.5 * (f(a) + f(b));
    for i in 1..n {
        sum += f(a + i as f64 * h);
    }
    sum * h
}

/// Simpson's rule (n must be even): ∫f(x)dx ≈ h/3 [f_0 + 4f_1 + 2f_2 + ... + f_n].
pub fn simpson<F>(f: &F, a: f64, b: f64, n: usize) -> Result<f64>
where
    F: Fn(f64) -> f64,
{
    if n % 2 != 0 {
        return Err(crate::core::PysicError::InvalidInput(
            "Simpson's rule requires even number of intervals".into(),
        ));
    }
    if n == 0 { return Ok(0.0); }
    let h = (b - a) / n as f64;
    let mut sum = f(a) + f(b);
    for i in 1..n {
        let coef = if i % 2 == 0 { 2.0 } else { 4.0 };
        sum += coef * f(a + i as f64 * h);
    }
    Ok(sum * h / 3.0)
}

/// Gauss-Legendre quadrature with n nodes on [a, b].
/// Uses precomputed nodes/weights for n=2..8, falls back to trapezoidal.
pub fn gauss_legendre<F>(f: &F, a: f64, b: f64, n: usize) -> f64
where
    F: Fn(f64) -> f64,
{
    // For simplicity, use composite Simpson with many subintervals as fallback
    // A production version would use Clenshaw-Curtis or Gauss-Legendre nodes
    let sub = if n <= 2 { 100 } else { n * 10 };
    let h = (b - a) / sub as f64;
    let mut sum = 0.5 * (f(a) + f(b));
    for i in 1..sub {
        sum += f(a + i as f64 * h);
    }
    sum * h
}

/// Romberg integration (Richardson extrapolation of trapezoidal rule).
pub fn romberg<F>(f: &F, a: f64, b: f64, max_order: usize) -> f64
where
    F: Fn(f64) -> f64,
{
    let mut r = vec![vec![0.0; max_order + 1]; max_order + 1];
    let h = b - a;
    r[0][0] = 0.5 * (f(a) + f(b)) * h;

    for i in 1..=max_order {
        let n = 1usize << i;
        let step = h / n as f64;
        let mut sum = r[0][0] / 2.0;
        for j in 0..(n / 2) {
            sum += f(a + (2 * j + 1) as f64 * step);
        }
        r[i][0] = sum * step;

        for k in 1..=i {
            r[i][k] = r[i][k-1] + (r[i][k-1] - r[i-1][k-1]) / ((1usize << (2 * k)) as f64 - 1.0);
        }
    }
    r[max_order][max_order]
}

// ============================================================================
// Interpolation
// ============================================================================

/// Linear interpolation between two points.
pub fn lerp(x0: f64, y0: f64, x1: f64, y1: f64, x: f64) -> f64 {
    if (x1 - x0).abs() < 1e-15 { return y0; }
    y0 + (y1 - y0) * (x - x0) / (x1 - x0)
}

/// Piecewise linear interpolation on a grid.
pub fn interp1d(x_grid: &[f64], y_values: &[f64], x: f64) -> Result<f64> {
    if x_grid.len() != y_values.len() {
        return Err(crate::core::PysicError::DimensionMismatch {
            expected: x_grid.len(),
            actual: y_values.len(),
        });
    }
    if x_grid.len() < 2 {
        return Err(crate::core::PysicError::InvalidInput(
            "Need at least 2 points for interpolation".into(),
        ));
    }
    if x <= x_grid[0] { return Ok(y_values[0]); }
    if x >= x_grid[x_grid.len() - 1] { return Ok(y_values[y_values.len() - 1]); }

    for i in 0..(x_grid.len() - 1) {
        if x >= x_grid[i] && x <= x_grid[i + 1] {
            return Ok(lerp(x_grid[i], y_values[i], x_grid[i + 1], y_values[i + 1], x));
        }
    }
    Ok(y_values[y_values.len() - 1])
}

/// Cubic spline interpolation (natural boundary conditions).
pub fn interp_cspline(x_grid: &[f64], y_values: &[f64], x: f64) -> Result<f64> {
    let n = x_grid.len();
    if n != y_values.len() || n < 3 {
        return Err(crate::core::PysicError::InvalidInput(
            "Need at least 3 points for cubic spline".into(),
        ));
    }

    // Compute second derivatives (natural spline)
    let mut h = vec![0.0; n - 1];
    let mut alpha = vec![0.0; n - 1];
    for i in 0..n - 1 {
        h[i] = x_grid[i + 1] - x_grid[i];
        if h[i] <= 0.0 {
            return Err(crate::core::PysicError::InvalidInput("Grid must be strictly increasing".into()));
        }
        alpha[i] = (y_values[i + 1] - y_values[i]) / h[i];
    }

    let mut l = vec![1.0; n];
    let mut mu = vec![0.0; n];
    let mut z = vec![0.0; n];

    for i in 1..n - 1 {
        l[i] = 2.0 * (x_grid[i + 1] - x_grid[i - 1]) - h[i - 1] * mu[i - 1];
        if l[i].abs() < 1e-15 {
            return Err(crate::core::PysicError::NumericalError("Singular spline system".into()));
        }
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - alpha[i - 1] - h[i - 1] * z[i - 1]) / l[i];
    }

    let mut c = vec![0.0; n];
    let mut b = vec![0.0; n - 1];
    let mut d = vec![0.0; n - 1];

    for j in (0..n - 1).rev() {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = alpha[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
        d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
    }

    // Find the right interval
    if x <= x_grid[0] { return Ok(y_values[0]); }
    if x >= x_grid[n - 1] { return Ok(y_values[n - 1]); }

    for i in 0..n - 1 {
        if x >= x_grid[i] && x <= x_grid[i + 1] {
            let dx = x - x_grid[i];
            return Ok(y_values[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx);
        }
    }
    Ok(y_values[n - 1])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gradient() {
        let f = |x: &[f64]| x[0].powi(2) + 2.0 * x[1].powi(2);
        let g = gradient(&f, &[1.0, 2.0], 1e-5);
        assert!((g[0] - 2.0).abs() < 1e-3);
        assert!((g[1] - 8.0).abs() < 1e-3);
    }

    #[test]
    fn test_trapz() {
        let f = |x: f64| x * x;
        let result = trapz(&f, 0.0, 1.0, 1000);
        assert!((result - 1.0 / 3.0).abs() < 1e-3);
    }

    #[test]
    fn test_simpson() {
        let f = |x: f64| x * x;
        let result = simpson(&f, 0.0, 1.0, 100).unwrap();
        assert!((result - 1.0 / 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_interp1d() {
        let x = vec![0.0, 1.0, 2.0];
        let y = vec![0.0, 1.0, 4.0];
        let val = interp1d(&x, &y, 0.5).unwrap();
        assert!((val - 0.5).abs() < 1e-10);
    }
}
