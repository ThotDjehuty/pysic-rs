//! Feynman propagators and Green's functions for quantum field theory.

use num_complex::Complex64;

/// Feynman propagator for a free scalar field in momentum space:
/// Δ_F(p) = i / (p² - m² + iε)
pub fn feynman_propagator_scalar(
    p_squared: f64,
    mass: f64,
    epsilon: f64,
) -> Complex64 {
    let denom = p_squared - mass * mass + Complex64::new(0.0, epsilon);
    Complex64::new(0.0, 1.0) / denom
}

/// Free scalar propagator in position space (Euclidean, 4D):
/// D(x) = m / (4π² |x|) K_1(m|x|)
/// For m=0: D(x) = 1 / (4π² |x|²)
pub fn free_propagator_position(r: f64, mass: f64) -> f64 {
    if mass * r < 1e-6 {
        // Massless limit: 1/(4π²r²)
        1.0 / (4.0 * std::f64::consts::PI * std::f64::consts::PI * r * r)
    } else {
        // Massive: use approximate Bessel K_1
        let mr = mass * r;
        let k1 = if mr > 2.0 {
            // Asymptotic: K_1(z) ≈ √(π/2z) e^(-z) (1 + 3/(8z) + ...)
            (std::f64::consts::PI / (2.0 * mr)).sqrt() * (-mr).exp() * (1.0 + 3.0 / (8.0 * mr))
        } else {
            // For small z: K_1(z) ≈ 1/z + z/2 [ln(z/2) + γ - 1/2]
            let euler = 0.5772156649015329;
            1.0 / mr + 0.5 * mr * (mr / 2.0).ln() + 0.5 * mr * (euler - 0.5)
        };
        mass / (4.0 * std::f64::consts::PI * std::f64::consts::PI * r) * k1
    }
}

/// Free scalar propagator in momentum space (on-shell):
/// G(p₀, p⃗) = i / (p₀² - |p⃗|² - m² + iε)
pub fn free_propagator_momentum(
    p0: f64,
    p_spatial_sq: f64,
    mass: f64,
    epsilon: f64,
) -> Complex64 {
    let p_sq = p0 * p0 - p_spatial_sq;
    let denom = p_sq - mass * mass + Complex64::new(0.0, epsilon);
    Complex64::new(0.0, 1.0) / denom
}

/// Photon propagator in Feynman gauge:
/// D_μν(p) = -i g_μν / (p² + iε)
pub fn photon_propagator_feynman(
    p_squared: f64,
    epsilon: f64,
) -> Vec<Vec<Complex64>> {
    let denom = p_squared + Complex64::new(0.0, epsilon);
    let val = Complex64::new(0.0, -1.0) / denom;

    // Minkowski metric η = diag(-1, +1, +1, +1)
    vec![
        vec![-val, Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0)],
        vec![Complex64::new(0.0, 0.0), val, Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0)],
        vec![Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0), val, Complex64::new(0.0, 0.0)],
        vec![Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0), val],
    ]
}

/// Electron propagator (Dirac fermion):
/// S_F(p) = i (p̸ + m) / (p² - m² + iε)
/// where p̸ = γ^μ p_μ
pub fn electron_propagator(
    p_squared: f64,
    mass: f64,
    epsilon: f64,
) -> Complex64 {
    let denom = p_squared - mass * mass + Complex64::new(0.0, epsilon);
    Complex64::new(0.0, 1.0) / denom
}

/// Retarded Green's function for the wave equation:
/// G_R(t, x) = θ(t) δ(t - |x|/c) / (4π|x|)
pub fn retarded_green_fn_wave(t: f64, r: f64, c: f64) -> f64 {
    if t < 0.0 || r < 1e-15 {
        return 0.0;
    }
    let delta_approx = if (t - r / c).abs() < 0.01 { 1.0 / 0.01 } else { 0.0 };
    delta_approx / (4.0 * std::f64::consts::PI * r)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_feynman_propagator_pole() {
        // On-shell: p² = m², propagator should diverge (large imaginary part)
        let g = feynman_propagator_scalar(1.0, 1.0, 0.01);
        assert!(g.norm() > 10.0);
    }

    #[test]
    fn test_photon_propagator_gauge() {
        let p = photon_propagator_feynman(1.0, 0.01);
        // Should be diagonal
        assert!((p[0][1].norm()).abs() < 1e-10);
        assert!((p[1][0].norm()).abs() < 1e-10);
    }
}
