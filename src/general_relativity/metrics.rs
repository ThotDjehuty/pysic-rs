//! Standard spacetime metrics: Schwarzschild, Kerr, FLRW, Minkowski.

use ndarray::Array2;

/// Minkowski metric: η_μν = diag(-1, +1, +1, +1).
pub fn minkowski_metric() -> Array2<f64> {
    let mut g = Array2::zeros((4, 4));
    g[[0, 0]] = -1.0;
    g[[1, 1]] = 1.0;
    g[[2, 2]] = 1.0;
    g[[3, 3]] = 1.0;
    g
}

/// Inverse Minkowski metric η^μν = diag(-1, +1, +1, +1).
pub fn minkowski_inverse() -> Array2<f64> {
    minkowski_metric()
}

/// Schwarzschild metric in Schwarzschild coordinates (t, r, θ, φ):
/// ds² = -(1 - 2GM/rc²) dt² + (1 - 2GM/rc²)^{-1} dr² + r² dΩ²
///
/// `m` is the mass in geometric units (G = c = 1).
pub fn schwarzschild_metric(r: f64, mass: f64) -> Array2<f64> {
    let rs = 2.0 * mass; // Schwarzschild radius
    let f = 1.0 - rs / r;

    let mut g = Array2::zeros((4, 4));
    g[[0, 0]] = -f;
    g[[1, 1]] = 1.0 / f;
    g[[2, 2]] = r * r;
    g[[3, 3]] = r * r; // assuming θ = π/2
    g
}

/// Inverse Schwarzschild metric.
pub fn schwarzschild_inverse(r: f64, mass: f64) -> Array2<f64> {
    let rs = 2.0 * mass;
    let f = 1.0 - rs / r;

    let mut g_inv = Array2::zeros((4, 4));
    g_inv[[0, 0]] = -1.0 / f;
    g_inv[[1, 1]] = f;
    g_inv[[2, 2]] = 1.0 / (r * r);
    g_inv[[3, 3]] = 1.0 / (r * r);
    g_inv
}

/// Kerr metric in Boyer-Lindquist coordinates (t, r, θ, φ):
/// ds² = -(1 - 2GMr/Σ) dt² - (4GMr a sin²θ/Σ) dt dφ + Σ/Δ dr² + Σ dθ²
///       + (r² + a² + 2GMr a² sin²θ/Σ) sin²θ dφ²
///
/// `m` = mass, `a` = spin parameter (a = J/Mc).
pub fn kerr_metric(r: f64, theta: f64, mass: f64, a: f64) -> Array2<f64> {
    let rs = 2.0 * mass;
    let sigma = r * r + a * a * theta.cos().powi(2);
    let delta = r * r - rs * r + a * a;
    let w = (2.0 * mass * r * a * theta.sin().powi(2)) / sigma;

    let mut g = Array2::zeros((4, 4));
    g[[0, 0]] = -(1.0 - rs * r / sigma);
    g[[0, 3]] = -w;
    g[[3, 0]] = -w;
    g[[1, 1]] = sigma / delta;
    g[[2, 2]] = sigma;
    g[[3, 3]] = (r * r + a * a + rs * r * a * a * theta.sin().powi(2) / sigma) * theta.sin().powi(2);
    g
}

/// FLRW metric for a flat universe:
/// ds² = -dt² + a(t)² (dx² + dy² + dz²)
///
/// `a` is the scale factor.
pub fn flrw_metric(a: f64) -> Array2<f64> {
    let mut g = Array2::zeros((4, 4));
    g[[0, 0]] = -1.0;
    g[[1, 1]] = a * a;
    g[[2, 2]] = a * a;
    g[[3, 3]] = a * a;
    g
}

/// Inverse FLRW metric.
pub fn flrw_inverse(a: f64) -> Array2<f64> {
    let a2_inv = 1.0 / (a * a);
    let mut g_inv = Array2::zeros((4, 4));
    g_inv[[0, 0]] = -1.0;
    g_inv[[1, 1]] = a2_inv;
    g_inv[[2, 2]] = a2_inv;
    g_inv[[3, 3]] = a2_inv;
    g_inv
}

/// Kerr-Newman metric (charged, rotating black hole).
pub fn kerr_newman_metric(r: f64, theta: f64, mass: f64, a: f64, charge: f64) -> Array2<f64> {
    let rs = 2.0 * mass;
    let r_q = charge * charge; // In geometric units
    let sigma = r * r + a * a * theta.cos().powi(2);
    let delta = r * r - rs * r + a * a + r_q;
    let w = (rs * r * a * theta.sin().powi(2)) / sigma;

    let mut g = Array2::zeros((4, 4));
    g[[0, 0]] = -(1.0 - rs * r / sigma);
    g[[0, 3]] = -w;
    g[[3, 0]] = -w;
    g[[1, 1]] = sigma / delta;
    g[[2, 2]] = sigma;
    g[[3, 3]] = (r * r + a * a + rs * r * a * a * theta.sin().powi(2) / sigma) * theta.sin().powi(2);
    g
}
