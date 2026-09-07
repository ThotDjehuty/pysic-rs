//! Geodesic equation and integration.

use crate::core::Result;

/// State vector for geodesic: [t, x, y, z, u^t, u^x, u^y, u^z].
pub type GeodesicState = Vec<f64>;

/// Right-hand side of the geodesic equation:
/// dx^μ/dτ = u^μ
/// du^μ/dτ = -Γ^μ_αβ u^α u^β
pub fn geodesic_equation_rhs(
    state: &[f64],
    christoffel: &[Vec<Vec<f64>>],
    n_dim: usize,
) -> Vec<f64> {
    let mut rhs = vec![0.0; 2 * n_dim];
    let x = &state[..n_dim];
    let u = &state[n_dim..];

    // dx^μ/dτ = u^μ
    for i in 0..n_dim {
        rhs[i] = u[i];
    }

    // du^μ/dτ = -Γ^μ_αβ u^α u^β
    for mu in 0..n_dim {
        let mut sum = 0.0;
        for alpha in 0..n_dim {
            for beta in 0..n_dim {
                sum += christoffel[mu][alpha][beta] * u[alpha] * u[beta];
            }
        }
        rhs[n_dim + mu] = -sum;
    }

    rhs
}

/// Integrate a geodesic using classical RK4.
///
/// # Arguments
/// * `initial_state` - Initial [x^μ, u^μ]
/// * `christoffel_fn` - Function returning Christoffel symbols at a point
/// * `n_dim` - Spacetime dimension (4 for GR)
/// * `dt` - Proper time step
/// * `n_steps` - Number of integration steps
pub fn integrate_geodesic<F>(
    initial_state: &GeodesicState,
    christoffel_fn: &F,
    n_dim: usize,
    dt: f64,
    n_steps: usize,
) -> Result<Vec<GeodesicState>>
where
    F: Fn(&[f64]) -> Vec<Vec<Vec<f64>>>,
{
    let mut trajectory = Vec::with_capacity(n_steps + 1);
    let mut state = initial_state.clone();
    trajectory.push(state.clone());

    for _ in 0..n_steps {
        let gamma = christoffel_fn(&state[..n_dim]);
        let k1 = geodesic_equation_rhs(&state, &gamma, n_dim);

        let state_k2: Vec<f64> = state.iter().zip(k1.iter()).map(|(s, k)| s + 0.5 * dt * k).collect();
        let gamma2 = christoffel_fn(&state_k2[..n_dim]);
        let k2 = geodesic_equation_rhs(&state_k2, &gamma2, n_dim);

        let state_k3: Vec<f64> = state.iter().zip(k2.iter()).map(|(s, k)| s + 0.5 * dt * k).collect();
        let gamma3 = christoffel_fn(&state_k3[..n_dim]);
        let k3 = geodesic_equation_rhs(&state_k3, &gamma3, n_dim);

        let state_k4: Vec<f64> = state.iter().zip(k3.iter()).map(|(s, k)| s + dt * k).collect();
        let gamma4 = christoffel_fn(&state_k4[..n_dim]);
        let k4 = geodesic_equation_rhs(&state_k4, &gamma4, n_dim);

        for i in 0..state.len() {
            state[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }

        trajectory.push(state.clone());
    }

    Ok(trajectory)
}

/// Compute the conserved quantities along a geodesic (for Schwarzschild):
/// E = -(1 - 2M/r) dt/dτ, L = r² dφ/dτ
pub fn geodesic_conserved_schw(
    state: &GeodesicState,
    mass: f64,
) -> (f64, f64) {
    let r = state[2]; // assuming x = (t, r, θ, φ) with θ=π/2
    let t_dot = state[4]; // dt/dτ
    let phi_dot = state[7]; // dφ/dτ
    let rs = 2.0 * mass;
    let e = -(1.0 - rs / r) * t_dot;
    let l = r * r * phi_dot;
    (e, l)
}
