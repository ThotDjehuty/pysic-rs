//! Berry phase and Berry curvature for quantum systems.

use num_complex::Complex64;

/// Berry phase for a closed path parameterized by t ∈ [0, 1]:
/// γ = Im ∮ ⟨n(t)|∇_R|n(t)⟩ · dR
///
/// Numerically: γ = -Im ln ∏_i ⟨n(R_i)|n(R_{i+1})⟩
pub fn berry_phase(
    states: &[Vec<Complex64>],  // |n(R_i)⟩ at each point
) -> f64 {
    let n = states.len();
    let mut overlap = Complex64::new(1.0, 0.0);

    for i in 0..n {
        let j = (i + 1) % n;
        // ⟨n(R_i)|n(R_{i+1})⟩
        let inner: Complex64 = states[i].iter().zip(states[j].iter()).map(|(a, b)| a.conj() * b).sum();
        overlap *= inner;
    }

    -overlap.arg()
}

/// Berry curvature F_{μν} = ∂_μ A_ν - ∂_ν A_μ
/// where A_μ = Im ⟨n|∂_μ n⟩ is the Berry connection.
pub fn berry_curvature(
    states_center: &[Complex64],
    states_mu_plus: &[Complex64],
    states_mu_minus: &[Complex64],
    states_nu_plus: &[Complex64],
    states_nu_minus: &[Complex64],
    d_mu: f64,
    d_nu: f64,
) -> f64 {
    // A_μ = Im ⟨n|∂_μ n⟩
    let inner_product = |a: &[Complex64], b: &[Complex64]| -> Complex64 {
        a.iter().zip(b.iter()).map(|(x, y)| x.conj() * y).sum()
    };

    let overlap = |a: &[Complex64], b: &[Complex64]| -> f64 {
        inner_product(a, b).arg()
    };

    let a_mu_plus = overlap(states_center, states_mu_plus);
    let a_mu_minus = overlap(states_center, states_mu_minus);
    let a_nu_plus = overlap(states_center, states_nu_plus);
    let a_nu_minus = overlap(states_center, states_nu_minus);

    let a_mu = (a_mu_plus - a_mu_minus) / (2.0 * d_mu);
    let a_nu = (a_nu_plus - a_nu_minus) / (2.0 * d_nu);

    // F_μν = ∂_μ A_ν - ∂_ν A_μ (approximated)
    let f_mu_nu = (a_nu_plus - a_nu_minus) / (2.0 * d_nu) - (a_mu_plus - a_mu_minus) / (2.0 * d_mu);
    f_mu_nu
}

/// Berry phase for a two-level system traversing a great circle on the Bloch sphere.
pub fn berry_phase_bloch(theta: f64) -> f64 {
    // For a pure state |ψ⟩ = cos(θ/2)|0⟩ + e^{iφ} sin(θ/2)|1⟩
    // traversing a loop that encloses a solid angle Ω = 2π(1 - cos θ)
    // the Berry phase is γ = -½ Ω = -π(1 - cos θ)
    -std::f64::consts::PI * (1.0 - theta.cos())
}

/// Spin-1/2 Berry phase: γ = -m Ω where m is the spin projection.
pub fn spin_berry_phase(m: f64, solid_angle: f64) -> f64 {
    -m * solid_angle
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    #[test]
    fn test_berry_phase_closed_loop() {
        // Uniform state: Berry phase should be zero
        let states: Vec<Vec<Complex64>> = (0..100)
            .map(|_| vec![Complex64::new(1.0, 0.0)])
            .collect();
        let phase = berry_phase(&states);
        assert!(phase.abs() < 1e-10);
    }

    #[test]
    fn test_berry_phase_bloch() {
        // Berry phase for θ = π/2 (great circle) should be -π
        let phase = berry_phase_bloch(std::f64::consts::FRAC_PI_2);
        assert!((phase - (-std::f64::consts::PI)).abs() < 1e-10);
    }
}
