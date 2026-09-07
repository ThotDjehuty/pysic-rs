//! Path integral methods: discretization, Wick rotation, effective potential.

use num_complex::Complex64;

/// Discretized path integral for a particle in a potential V(x).
///
/// Z = ∫ Dx exp(-S[x]/ℏ) where S = ∫ (½m ẋ² + V(x)) dt.
///
/// Approximates as a product of free-particle propagators and potential weights.
pub fn discretized_path_integral(
    x_initial: f64,
    x_final: f64,
    v: &dyn Fn(f64) -> f64,
    mass: f64,
    hbar: f64,
    n_slices: usize,
    t_total: f64,
) -> Complex64 {
    let dt = t_total / n_slices as f64;
    let dx = (x_final - x_initial) / n_slices as f64;
    let prefactor = (mass / (2.0 * std::f64::consts::PI * hbar * dt)).sqrt();
    let mut weight = Complex64::new(1.0, 0.0);

    let mut x = x_initial;
    for _ in 0..n_slices {
        let x_next = x + dx;
        // Free particle propagator
        let free = prefactor * Complex64::new(0.0, -mass * dx * dx / (2.0 * hbar * dt)).exp();
        // Potential weight
        let pot = Complex64::new(0.0, -v(x) * dt / hbar).exp();
        weight *= free * pot;
        x = x_next;
    }

    weight
}

/// Compute the effective potential via Legendre transform of the free energy.
/// V_eff(φ) = Γ(φ) / β where Γ is the 1PI effective action.
pub fn effective_potential(
    phi: &[f64],
    v_tree: &dyn Fn(f64) -> f64,
    loop_correction: &dyn Fn(f64) -> f64,
) -> Vec<f64> {
    phi.iter().map(|&p| v_tree(p) + loop_correction(p)).collect()
}

/// Wick rotation: t → -iτ (real time to Euclidean time).
/// Returns the Euclidean action S_E = ∫ (½m (dx/dτ)² + V(x)) dτ.
pub fn wick_rotation(
    x_trajectory: &[f64],
    v: &dyn Fn(f64) -> f64,
    mass: f64,
    dt: f64,
) -> f64 {
    let mut s_e = 0.0;
    let n = x_trajectory.len();

    for i in 0..n - 1 {
        let dx_dt = (x_trajectory[i + 1] - x_trajectory[i]) / dt;
        s_e += (0.5 * mass * dx_dt * dx_dt + v(x_trajectory[i])) * dt;
    }

    s_e
}

/// Instanton action: S_inst = ∫_{x_min}^{x_max} √(2m V(x)) dx.
/// For a double-well potential, this gives the tunneling rate.
pub fn instanton_action(
    v: &dyn Fn(f64) -> f64,
    mass: f64,
    x_min: f64,
    x_max: f64,
    n_slices: usize,
) -> f64 {
    let dx = (x_max - x_min) / n_slices as f64;
    let mut s = 0.0;

    for i in 0..n_slices {
        let x = x_min + (i as f64 + 0.5) * dx;
        let v_val = v(x).abs();
        s += (2.0 * mass * v_val).sqrt() * dx;
    }

    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_path_integral_harmonic() {
        // Harmonic oscillator V(x) = ½ω²x²
        let v = |x: f64| 0.5 * x * x;
        let z = discretized_path_integral(0.0, 0.0, &v, 1.0, 1.0, 100, 1.0);
        // Should be non-zero and real for this symmetric case
        assert!(z.norm() > 0.0);
    }

    #[test]
    fn test_instanton_double_well() {
        let v = |x: f64| (x * x - 1.0).powi(2);
        let s = instanton_action(&v, 1.0, -1.0, 1.0, 1000);
        assert!(s > 0.0);
    }
}
