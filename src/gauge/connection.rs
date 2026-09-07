//! Gauge connection and covariant derivative.

/// Gauge connection A_μ (Lie algebra-valued, stored as a flat vector of n × n matrices).
///
/// For SU(2): A_μ = A_μ^a σ^a/2i, stored as (n_points × n_components × dim × dim).
pub fn gauge_connection(
    a_field: &[f64],         // A_μ^a(x) at each point: [n_points][n_components]
    n_components: usize,     // Number of gauge field components
    dim: usize,              // Matrix dimension (2 for SU(2), 3 for SU(3))
    n_points: usize,
) -> Vec<Vec<Vec<Vec<f64>>>> {
    // Return A_μ at each point as dim×dim matrices
    let mut result = vec![vec![vec![vec![0.0; dim]; dim]; n_components]; n_points];

    // For SU(2), A_μ = A_μ^a σ^a / (2i)
    if dim == 2 && n_components >= 3 {
        for p in 0..n_points {
            for mu in 0..n_components.min(4) {
                let a0 = a_field[p * n_components + 0.min(n_components - 1)];
                let a1 = a_field[p * n_components + 1.min(n_components - 1)];
                let a2 = a_field[p * n_components + 2.min(n_components - 1)];

                // A_μ = (a0 σ_x + a1 σ_y + a2 σ_z) / 2
                result[p][mu][0][0] = a2 / 2.0;
                result[p][mu][0][1] = (a0 - a1 * 0.0) / 2.0;  // simplified
                result[p][mu][1][0] = (a0 + a1 * 0.0) / 2.0;
                result[p][mu][1][1] = -a2 / 2.0;
            }
        }
    }

    result
}

/// Covariant derivative: D_μ φ = ∂_μ φ + i g A_μ φ.
///
/// For a scalar field φ (real-valued, simplified).
pub fn covariant_derivative(
    phi: &[f64],
    a_mu: &[f64],         // A_μ at the point
    dx: &[f64],           // Grid spacings
    coupling: f64,
    point_idx: usize,
    n_grid: usize,
) -> Vec<f64> {
    let n_comp = a_mu.len();
    let mut result = vec![0.0; n_comp];

    for mu in 0..n_comp {
        // ∂_μ φ (central difference)
        let idx_plus = (point_idx + n_grid).min(n_grid * n_grid - 1);
        let idx_minus = point_idx.saturating_sub(n_grid);
        let dphi = (phi[idx_plus] - phi[idx_minus]) / (2.0 * dx[mu.min(dx.len() - 1)]);

        // D_μ φ = ∂_μ φ + i g A_μ φ (imaginary part dropped for real fields)
        result[mu] = dphi + coupling * a_mu[mu] * phi[point_idx];
    }

    result
}

/// Parallel transport of a vector along a path using the gauge connection.
pub fn parallel_transport(
    vector: &[f64],
    connection_along_path: &[Vec<Vec<f64>>],  // A_μ at each point
    path: &[(f64, f64)],                       // (x, y) coordinates
    dt: f64,
    coupling: f64,
) -> Vec<f64> {
    let mut v = vector.to_vec();

    for i in 0..path.len() - 1 {
        let dx = path[i + 1].0 - path[i].0;
        let dy = path[i + 1].1 - path[i].1;

        // U = exp(-i g A_μ dx^μ) ≈ I - i g A_μ dx^μ
        let a = &connection_along_path[i];
        let phase = coupling * (a[0][0] * dx + a[1][0] * dy);

        // Apply rotation
        let cos_p = phase.cos();
        let sin_p = phase.sin();
        if v.len() >= 2 {
            let new_x = cos_p * v[0] - sin_p * v[1];
            let new_y = sin_p * v[0] + cos_p * v[1];
            v[0] = new_x;
            v[1] = new_y;
        }
    }

    v
}
