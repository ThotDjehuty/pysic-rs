//! Christoffel symbols, Riemann tensor, Ricci tensor, Einstein tensor.

use ndarray::Array2;

/// Compute Christoffel symbols Γ^λ_μν from the metric tensor.
///
/// Γ^λ_μν = ½ g^{λρ} (∂_μ g_{νρ} + ∂_ν g_{μρ} - ∂_ρ g_{μν})
///
/// Uses finite differences to compute metric derivatives.
pub fn christoffel_from_metric(
    metric: &dyn Fn(&[f64]) -> Array2<f64>,
    coords: &[f64],
    h: f64,
) -> Vec<Vec<Vec<f64>>> {
    let n = coords.len();
    let g = metric(coords);

    // Compute inverse metric
    let g_inv = invert_metric(&g);

    // Compute metric derivatives
    let mut dg = vec![vec![vec![0.0; n]; n]; n]; // dg[i][j][k] = ∂_i g_{jk}
    for i in 0..n {
        let mut coords_plus = coords.to_vec();
        let mut coords_minus = coords.to_vec();
        coords_plus[i] += h;
        coords_minus[i] -= h;
        let g_plus = metric(&coords_plus);
        let g_minus = metric(&coords_minus);
        for j in 0..n {
            for k in 0..n {
                dg[i][j][k] = (g_plus[[j, k]] - g_minus[[j, k]]) / (2.0 * h);
            }
        }
    }

    // Christoffel symbols
    let mut gamma = vec![vec![vec![0.0; n]; n]; n];
    for lam in 0..n {
        for mu in 0..n {
            for nu in 0..n {
                let mut sum = 0.0;
                for rho in 0..n {
                    sum += 0.5 * g_inv[[lam, rho]]
                        * (dg[mu][nu][rho] + dg[nu][mu][rho] - dg[rho][mu][nu]);
                }
                gamma[lam][mu][nu] = sum;
            }
        }
    }

    gamma
}

/// Riemann tensor R^ρ_σμν from Christoffel symbols and their derivatives.
pub fn riemann_tensor(
    metric: &dyn Fn(&[f64]) -> Array2<f64>,
    coords: &[f64],
    h: f64,
) -> Vec<Vec<Vec<Vec<f64>>>> {
    let n = coords.len();
    let gamma = christoffel_from_metric(metric, coords, h);

    // Compute derivatives of Christoffel symbols
    let mut dgamma = vec![vec![vec![vec![0.0; n]; n]; n]; n]; // dgamma[i][lam][mu][nu]
    for i in 0..n {
        let mut coords_plus = coords.to_vec();
        let mut coords_minus = coords.to_vec();
        coords_plus[i] += h;
        coords_minus[i] -= h;
        let gamma_plus = christoffel_from_metric(metric, &coords_plus, h);
        let gamma_minus = christoffel_from_metric(metric, &coords_minus, h);
        for lam in 0..n {
            for mu in 0..n {
                for nu in 0..n {
                    dgamma[i][lam][mu][nu] = (gamma_plus[lam][mu][nu] - gamma_minus[lam][mu][nu]) / (2.0 * h);
                }
            }
        }
    }

    // R^ρ_σμν = ∂_μ Γ^ρ_νσ - ∂_ν Γ^ρ_μσ + Γ^ρ_μλ Γ^λ_νσ - Γ^ρ_νλ Γ^λ_μσ
    let mut riemann = vec![vec![vec![vec![0.0; n]; n]; n]; n];
    for rho in 0..n {
        for sigma in 0..n {
            for mu in 0..n {
                for nu in 0..n {
                    let mut val = dgamma[mu][rho][nu][sigma] - dgamma[nu][rho][mu][sigma];
                    for lam in 0..n {
                        val += gamma[rho][mu][lam] * gamma[lam][nu][sigma]
                            - gamma[rho][nu][lam] * gamma[lam][mu][sigma];
                    }
                    riemann[rho][sigma][mu][nu] = val;
                }
            }
        }
    }

    riemann
}

/// Ricci tensor R_μν = R^λ_μλν.
pub fn ricci_tensor(
    metric: &dyn Fn(&[f64]) -> Array2<f64>,
    coords: &[f64],
    h: f64,
) -> Array2<f64> {
    let n = coords.len();
    let riemann = riemann_tensor(metric, coords, h);
    let mut ricci = Array2::zeros((n, n));

    for mu in 0..n {
        for nu in 0..n {
            let mut sum = 0.0;
            for lam in 0..n {
                sum += riemann[lam][mu][lam][nu];
            }
            ricci[[mu, nu]] = sum;
        }
    }

    ricci
}

/// Ricci scalar R = g^{μν} R_μν.
pub fn ricci_scalar(
    metric: &dyn Fn(&[f64]) -> Array2<f64>,
    coords: &[f64],
    h: f64,
) -> f64 {
    let g = metric(coords);
    let g_inv = invert_metric(&g);
    let ricci = ricci_tensor(metric, coords, h);
    let n = coords.len();

    let mut r = 0.0;
    for mu in 0..n {
        for nu in 0..n {
            r += g_inv[[mu, nu]] * ricci[[mu, nu]];
        }
    }
    r
}

/// Einstein tensor G_μν = R_μν - ½ g_μν R.
pub fn einstein_tensor(
    metric: &dyn Fn(&[f64]) -> Array2<f64>,
    coords: &[f64],
    h: f64,
) -> Array2<f64> {
    let g = metric(coords);
    let n = coords.len();
    let ricci = ricci_tensor(metric, coords, h);
    let r = ricci_scalar(metric, coords, h);

    let mut einstein = Array2::zeros((n, n));
    for mu in 0..n {
        for nu in 0..n {
            einstein[[mu, nu]] = ricci[[mu, nu]] - 0.5 * g[[mu, nu]] * r;
        }
    }
    einstein
}

/// Simple 4x4 matrix inverse (Gauss-Jordan).
fn invert_metric(m: &Array2<f64>) -> Array2<f64> {
    let n = m.nrows();
    let mut aug = Array2::zeros((n, 2 * n));
    for i in 0..n {
        for j in 0..n {
            aug[[i, j]] = m[[i, j]];
        }
        aug[[i, n + i]] = 1.0;
    }

    for col in 0..n {
        // Partial pivoting
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
            continue;
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
    inv
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::general_relativity::metrics::{schwarzschild_metric_full, minkowski_metric};

    #[test]
    fn test_minkowski_christoffel_zero() {
        let mink = |_c: &[f64]| minkowski_metric();
        let gamma = christoffel_from_metric(&mink, &[0.0, 1.0, 1.0, 0.0], 1e-5);
        for lam in 0..4 {
            for mu in 0..4 {
                for nu in 0..4 {
                    assert!(gamma[lam][mu][nu].abs() < 1e-4,
                        "Non-zero Christoffel for Minkowski: gamma[{}][{}][{}] = {}", lam, mu, nu, gamma[lam][mu][nu]);
                }
            }
        }
    }

    /// Schwarzschild is a vacuum solution, so R ≡ 0 everywhere outside the
    /// horizon — at every radius *and* every polar angle.
    ///
    /// The metric closure must be a genuine function of the coordinates: a
    /// closure that ignores its argument has zero derivatives and passes
    /// trivially. It must also carry the `sin²θ` factor in `g_φφ`, otherwise
    /// the angular terms fail to cancel and R comes out as -2/r².
    #[test]
    fn test_schwarzschild_ricci_flat() {
        let metric = |c: &[f64]| schwarzschild_metric_full(c[1], c[2], 1.0);

        for &r in &[6.0_f64, 10.0, 25.0, 100.0] {
            for &theta in &[
                std::f64::consts::FRAC_PI_4,
                std::f64::consts::FRAC_PI_2,
                2.0,
            ] {
                let rs = ricci_scalar(&metric, &[0.0, r, theta, 0.0], 1e-3);
                // Curvature scale at this radius is ~2M/r³; demand the vacuum
                // residual sit far below it rather than below a fixed 0.1.
                let scale = 2.0 / (r * r * r);
                assert!(
                    rs.abs() < 1e-3 * scale.max(1e-8),
                    "R should vanish for Schwarzschild at r={r}, theta={theta}; got {rs:e}"
                );
            }
        }
    }

    /// Guards the specific regression: freezing θ = π/2 in the metric (so that
    /// ∂_θ g_φφ = 0) produces R = -2/r² on a vacuum solution.
    #[test]
    fn test_equatorial_only_metric_is_not_ricci_flat() {
        use crate::general_relativity::metrics::schwarzschild_metric;
        let frozen = |c: &[f64]| schwarzschild_metric(c[1], 1.0);
        let r0 = 10.0_f64;
        let rs = ricci_scalar(&frozen, &[0.0, r0, std::f64::consts::FRAC_PI_2, 0.0], 1e-3);
        assert!(
            (rs - (-2.0 / (r0 * r0))).abs() < 1e-4,
            "expected the documented -2/r² artefact, got {rs:e}"
        );
    }
}
