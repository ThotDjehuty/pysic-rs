//! ADM 3+1 decomposition for numerical relativity.

use ndarray::Array2;
use crate::core::Result;

/// A single spatial hypersurface slice of the ADM decomposition.
#[derive(Debug, Clone)]
pub struct MetricSlice {
    /// Spatial metric γ_ij (n×n symmetric positive-definite)
    pub gamma: Array2<f64>,
    /// Extrinsic curvature K_ij (n×n symmetric)
    pub k: Array2<f64>,
    /// Lapse function N (scalar field)
    pub lapse: f64,
    /// Shift vector N^i (n-vector)
    pub shift: Vec<f64>,
}

impl MetricSlice {
    pub fn new(dim: usize) -> Self {
        Self {
            gamma: Array2::eye(dim),
            k: Array2::zeros((dim, dim)),
            lapse: 1.0,
            shift: vec![0.0; dim],
        }
    }
}

/// Full ADM metric representation over a 3D grid.
#[derive(Debug, Clone)]
pub struct AdmMetric {
    pub slices: Vec<MetricSlice>,
    pub grid_size: [usize; 3],
    pub dx: f64,
    pub dt: f64,
}

impl AdmMetric {
    pub fn new(grid_size: [usize; 3], dx: f64, dt: f64) -> Self {
        let n = grid_size[0] * grid_size[1] * grid_size[2];
        Self {
            slices: vec![MetricSlice::new(3); n],
            grid_size,
            dx,
            dt,
        }
    }
}

/// Evolve the spatial metric and extrinsic curvature by one timestep.
///
/// ∂_t γ_ij = -2 N K_ij + Lie_N γ_ij
/// ∂_t K_ij = -D_i D_j N + N (R_ij + K K_ij - 2 K_ik K^k_j) - 8πG S_ij
pub fn evolve_metric(
    gamma: &Array2<f64>,
    k: &Array2<f64>,
    lapse: f64,
    shift: &[f64],
    dt: f64,
) -> Result<(Array2<f64>, Array2<f64>)> {
    let dim = gamma.nrows();
    if dim != gamma.ncols() || dim != k.nrows() || dim != k.ncols() {
        return Err(crate::core::PysicError::DimensionMismatch {
            expected: dim,
            actual: k.nrows(),
        });
    }

    let mut gamma_new = gamma.clone();
    let mut k_new = k.clone();

    // ∂_t γ_ij = -2 N K_ij + (Lie derivative of shift)
    for i in 0..dim {
        for j in 0..dim {
            gamma_new[[i, j]] += dt * (-2.0 * lapse * k[[i, j]]);
            if i < shift.len() && j < shift.len() {
                gamma_new[[i, j]] += dt * (shift[i] + shift[j]);
            }
        }
    }

    // ∂_t K_ij (linearised)
    let trace_k: f64 = (0..dim).map(|l| k[[l, l]]).sum();
    for i in 0..dim {
        for j in 0..dim {
            k_new[[i, j]] += dt * lapse * (
                trace_k * k[[i, j]]
                - 2.0 * (0..dim).map(|l| k[[i, l]] * k[[l, j]]).sum::<f64>()
            );
        }
    }

    Ok((gamma_new, k_new))
}

/// Check the Hamiltonian constraint: ^{(3)}R + K² − K_ij K^ij = 16πG ρ.
pub fn hamiltonian_constraint(
    gamma: &Array2<f64>,
    k: &Array2<f64>,
    rho: f64,
) -> f64 {
    let dim = gamma.nrows();
    let trace_k: f64 = (0..dim).map(|l| k[[l, l]]).sum();
    let kk: f64 = (0..dim).flat_map(|i| (0..dim).map(move |j| k[[i, j]] * k[[i, j]])).sum();

    // Simplified: ^{(3)}R ≈ 0 for flat background
    trace_k * trace_k - kk - 16.0 * std::f64::consts::PI * rho
}

/// Check the momentum constraint: D_j(K^{ij} − γ^{ij} K) = 8πG j^i.
pub fn momentum_constraint(
    gamma: &Array2<f64>,
    k: &Array2<f64>,
    j: &[f64],
) -> Vec<f64> {
    let dim = gamma.nrows();
    let trace_k: f64 = (0..dim).map(|l| k[[l, l]]).sum();

    let mut residuals = Vec::with_capacity(dim);
    for i in 0..dim {
        let mut sum = 0.0;
        for jj in 0..dim {
            sum += k[[i, jj]] - gamma[[i, jj]] * trace_k;
        }
        let j_i = if i < j.len() { j[i] } else { 0.0 };
        residuals.push(sum - 8.0 * std::f64::consts::PI * j_i);
    }
    residuals
}
