//! 3D Maxwell FDTD solver (Yee grid).
//!
//! Solves Maxwell's equations:
//!   ∂E/∂t = (1/ε) ∇×H - σ/ε E
//!   ∂H/∂t = -(1/μ) ∇×E
//!
//! Uses the Yee staggered grid: E on integer nodes, H on half-integer nodes.

/// 3D Maxwell FDTD state.
pub struct MaxwellState {
    /// Ex on (i+1/2, j, k) grid points
    pub ex: Vec<f64>,
    /// Ey on (i, j+1/2, k)
    pub ey: Vec<f64>,
    /// Ez on (i, j, k+1/2)
    pub ez: Vec<f64>,
    /// Hx on (i, j+1/2, k+1/2)
    pub hx: Vec<f64>,
    /// Hy on (i+1/2, j, k+1/2)
    pub hy: Vec<f64>,
    /// Hz on (i+1/2, j+1/2, k)
    pub hz: Vec<f64>,
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
}

impl MaxwellState {
    pub fn new(nx: usize, ny: usize, nz: usize) -> Self {
        Self {
            ex: vec![0.0; nx * ny * nz],
            ey: vec![0.0; nx * ny * nz],
            ez: vec![0.0; nx * ny * nz],
            hx: vec![0.0; nx * ny * nz],
            hy: vec![0.0; nx * ny * nz],
            hz: vec![0.0; nx * ny * nz],
            nx,
            ny,
            nz,
        }
    }
}

/// One step of the 3D Maxwell FDTD update (Yee algorithm).
///
/// # Arguments
/// * `state` - Current E and H fields
/// * `dt` - Time step (must satisfy CFL: c·dt ≤ dx/√3)
/// * `dx`, `dy`, `dz` - Spatial grid spacings
/// * `eps` - Permittivity (scalar, uniform)
/// * `mu` - Permeability (scalar, uniform)
/// * `sigma` - Conductivity (scalar, uniform)
pub fn maxwell_step_3d(
    state: &mut MaxwellState,
    dt: f64,
    dx: f64,
    dy: f64,
    dz: f64,
    eps: f64,
    mu: f64,
    sigma: f64,
) {
    let nx = state.nx;
    let ny = state.ny;
    let nz = state.nz;

    // Update H from E (half step)
    let dt_mu = dt / mu;
    for k in 0..nz - 1 {
        for j in 0..ny - 1 {
            for i in 0..nx - 1 {
                let idx = k * ny * nx + j * nx + i;
                // ∂Hx/∂t = -(1/μ)(∂Ey/∂z - ∂Ez/∂y)
                state.hx[idx] -= dt_mu * (
                    (state.ey[(k + 1) * ny * nx + j * nx + i] - state.ey[idx]) / dz
                    - (state.ez[k * ny * nx + (j + 1) * nx + i] - state.ez[idx]) / dy
                );
                // ∂Hy/∂t = -(1/μ)(∂Ez/∂x - ∂Ex/∂z)
                state.hy[idx] -= dt_mu * (
                    (state.ez[k * ny * nx + j * nx + i + 1] - state.ez[idx]) / dx
                    - (state.ex[(k + 1) * ny * nx + j * nx + i] - state.ex[idx]) / dz
                );
                // ∂Hz/∂t = -(1/μ)(∂Ex/∂y - ∂Ey/∂x)
                state.hz[idx] -= dt_mu * (
                    (state.ex[k * ny * nx + (j + 1) * nx + i] - state.ex[idx]) / dy
                    - (state.ey[k * ny * nx + j * nx + i + 1] - state.ey[idx]) / dx
                );
            }
        }
    }

    // Update E from H (full step)
    let dt_eps = dt / eps;
    let sigma_factor = 1.0 - sigma * dt / (2.0 * eps);

    for k in 1..nz {
        for j in 1..ny {
            for i in 1..nx {
                let idx = k * ny * nx + j * nx + i;
                let prev = (k - 1) * ny * nx + (j - 1) * nx + (i - 1);

                // ∂Ex/∂t = (1/ε)(∂Hz/∂y - ∂Hy/∂z) - σ/ε Ex
                state.ex[idx] = sigma_factor * state.ex[idx] + dt_eps * (
                    (state.hz[idx] - state.hz[prev]) / dy
                    - (state.hy[idx] - state.hy[prev]) / dz
                );
                // ∂Ey/∂t = (1/ε)(∂Hx/∂z - ∂Hz/∂x) - σ/ε Ey
                state.ey[idx] = sigma_factor * state.ey[idx] + dt_eps * (
                    (state.hx[idx] - state.hx[prev]) / dz
                    - (state.hz[idx] - state.hz[prev]) / dx
                );
                // ∂Ez/∂t = (1/ε)(∂Hy/∂x - ∂Hx/∂y) - σ/ε Ez
                state.ez[idx] = sigma_factor * state.ez[idx] + dt_eps * (
                    (state.hy[idx] - state.hy[prev]) / dx
                    - (state.hx[idx] - state.hx[prev]) / dy
                );
            }
        }
    }
}

/// Get the maximum stable time step for Maxwell FDTD (CFL condition).
pub fn maxwell_max_dt(dx: f64, dy: f64, dz: f64, c: f64) -> f64 {
    let inv = (1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)).sqrt();
    1.0 / (c * inv)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_maxwell_cfl() {
        let dt = maxwell_max_dt(0.01, 0.01, 0.01, 3e8);
        assert!(dt > 0.0);
        assert!(dt < 1.0);
    }
}
