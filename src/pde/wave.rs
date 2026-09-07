//! Wave equation solvers: ∂²u/∂t² = v² ∇²u.
//!
//! Finite-difference time-domain (FDTD) method.

/// 1D wave equation via FDTD (leapfrog in time).
///
/// ∂²u/∂t² = v² ∂²u/∂x²
///
/// # Arguments
/// * `u0` - Initial displacement
/// * `v0` - Initial velocity (usually zero)
/// * `dx` - Spatial grid spacing
/// * `dt` - Time step (must satisfy CFL: dt ≤ dx/v)
/// * `v` - Wave speed
/// * `n_steps` - Number of time steps
pub fn wave_fdtd_1d(
    u0: &[f64],
    v0: &[f64],
    dx: f64,
    dt: f64,
    v: f64,
    n_steps: usize,
) -> Vec<f64> {
    let n = u0.len();
    let c2 = (v * dt / dx).powi(2);

    // CFL check
    if c2 > 1.0 {
        eprintln!("Warning: CFL number > 1, solution may be unstable");
    }

    let mut u_prev = u0.to_vec();
    let mut u_curr = vec![0.0; n];

    // Initialize u_curr from u0 and v0 (first time step)
    for i in 1..n - 1 {
        u_curr[i] = u_prev[i] + dt * v0[i] + 0.5 * c2 * (u_prev[i - 1] - 2.0 * u_prev[i] + u_prev[i + 1]);
    }
    u_curr[0] = u_prev[0];
    u_curr[n - 1] = u_prev[n - 1];

    let mut u_next = vec![0.0; n];

    for _step in 0..n_steps {
        for i in 1..n - 1 {
            u_next[i] = 2.0 * u_curr[i] - u_prev[i] + c2 * (u_curr[i - 1] - 2.0 * u_curr[i] + u_curr[i + 1]);
        }
        // Dirichlet BC
        u_next[0] = 0.0;
        u_next[n - 1] = 0.0;

        u_prev = u_curr;
        u_curr = u_next.clone();
    }

    u_curr
}

/// 2D wave equation via FDTD (leapfrog in time).
///
/// ∂²u/∂t² = v² (∂²u/∂x² + ∂²u/∂y²)
pub fn wave_fdtd_2d(
    u0: &[f64],
    v0: &[f64],
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    dt: f64,
    v: f64,
    n_steps: usize,
) -> Vec<f64> {
    let cx2 = (v * dt / dx).powi(2);
    let cy2 = (v * dt / dy).powi(2);

    if cx2 + cy2 > 1.0 {
        eprintln!("Warning: CFL number > 1, solution may be unstable");
    }

    let mut u_prev = u0.to_vec();
    let mut u_curr = vec![0.0; nx * ny];

    // Initialize
    for j in 1..ny - 1 {
        for i in 1..nx - 1 {
            let idx = j * nx + i;
            u_curr[idx] = u_prev[idx] + dt * v0[idx]
                + 0.5 * (cx2 * (u_prev[idx - 1] - 2.0 * u_prev[idx] + u_prev[idx + 1])
                       + cy2 * (u_prev[(j - 1) * nx + i] - 2.0 * u_prev[idx] + u_prev[(j + 1) * nx + i]));
        }
    }

    let mut u_next = vec![0.0; nx * ny];

    for _step in 0..n_steps {
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = j * nx + i;
                u_next[idx] = 2.0 * u_curr[idx] - u_prev[idx]
                    + cx2 * (u_curr[idx - 1] - 2.0 * u_curr[idx] + u_curr[idx + 1])
                    + cy2 * (u_curr[(j - 1) * nx + i] - 2.0 * u_curr[idx] + u_curr[(j + 1) * nx + i]);
            }
        }
        u_prev = u_curr;
        u_curr = u_next.clone();
    }

    u_curr
}
