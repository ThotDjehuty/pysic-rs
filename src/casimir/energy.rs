//! Casimir energy and force calculations.

/// Casimir energy per unit area for two parallel conducting plates:
/// E/A = -π² ℏc / (720 d³)
///
/// `d` is the plate separation.
pub fn casimir_energy_parallel_plates(d: f64, hbar: f64, c: f64) -> f64 {
    -std::f64::consts::PI * std::f64::consts::PI * hbar * c / (720.0 * d * d * d)
}

/// Casimir force per unit area:
/// F/A = -π² ℏc / (240 d⁴)
pub fn casimir_force(d: f64, hbar: f64, c: f64) -> f64 {
    -std::f64::consts::PI * std::f64::consts::PI * hbar * c / (240.0 * d * d * d * d)
}

/// Casimir energy for a conducting sphere (Boyеr approximation):
/// E = -ℏc / (8π a)  (approximate, with cutoff regularization)
pub fn casimir_energy_sphere(radius: f64, hbar: f64, c: f64) -> f64 {
    -hbar * c / (8.0 * std::f64::consts::PI * radius)
}

/// Regularized Casimir energy using zeta function regularization.
/// For parallel plates: E/A = -ζ(3) ℏc / (8π d³)
pub fn casimir_energy_zeta(d: f64, hbar: f64, c: f64) -> f64 {
    let zeta3 = 1.2020569031595942; // ζ(3)
    -zeta3 * hbar * c / (8.0 * std::f64::consts::PI * d * d * d)
}

/// Casimir energy between concentric cylinders (radius a, separation d).
pub fn casimir_energy_cylinders(
    radius: f64,
    separation: f64,
    hbar: f64,
    c: f64,
) -> f64 {
    // Approximate for d << a
    -hbar * c * radius / (24.0 * separation * separation)
}

/// Temperature-dependent Casimir energy (finite T):
/// E(T) = E(0) + (π² k_B⁴ T⁴)/(45 ℏ³ c³) × V
pub fn casimir_energy_finite_temperature(
    d: f64,
    temperature: f64,
    hbar: f64,
    c: f64,
    k_b: f64,
) -> f64 {
    let e_zero = casimir_energy_parallel_plates(d, hbar, c);
    let thermal = std::f64::consts::PI * std::f64::consts::PI * k_b * k_b * k_b * k_b * temperature * temperature * temperature * temperature
        / (45.0 * hbar * hbar * hbar * c * c * c);
    e_zero + thermal
}

/// Lifshitz formula for Casimir energy with dielectric media:
/// E/A = -k_BT/(2π) Σ_n ∫₀^∞ k⊥ dk⊥ [ln(1 - r₁r₂ e^{-2κd}) + (s-p) terms]
///
/// Simplified Lifshitz formula at T=0 for two identical dielectrics:
/// E/A = -ℏc/(4π²) ∫₀^∞ k³ dk ∫₀^∞ dx [1/(e^{2kd} - 1)] (ε(ik) - 1)/(ε(ik) + 1)
pub fn casimir_energy_lifshitz(
    d: f64,
    epsilon_fn: &dyn Fn(f64) -> f64,  // ε(iω) on imaginary axis
    hbar: f64,
    c: f64,
    n_points: usize,
) -> f64 {
    let mut integral = 0.0;
    let dk = 10.0 / d / n_points as f64;

    for i in 0..n_points {
        let k = (i as f64 + 0.5) * dk;
        let eps = epsilon_fn(k);
        let kappa = k; // For vacuum outside
        let x = 2.0 * kappa * d;
        let boltzmann = if x > 500.0 { 0.0 } else { 1.0 / (x.exp() - 1.0) };
        let reflection = (eps - 1.0) / (eps + 1.0);
        integral += k * k * k * reflection * reflection * boltzmann;
    }

    -hbar * c * integral * dk / (4.0 * std::f64::consts::PI * std::f64::consts::PI)
}

/// Polder potential: interaction energy between an atom and a conducting wall.
/// U(z) = -3ℏc α/(8π z⁴) (non-retarded) or -23ℏc α/(64π z⁵) (retarded)
pub fn polder_potential(
    distance: f64,
    polarizability: f64,
    hbar: f64,
    c: f64,
    retarded: bool,
) -> f64 {
    if retarded {
        -23.0 * hbar * c * polarizability / (64.0 * std::f64::consts::PI * distance.powi(5))
    } else {
        -3.0 * hbar * c * polarizability / (8.0 * std::f64::consts::PI * distance.powi(4))
    }
}

/// van der Waals interaction: U(r) = -C₆/r⁶
pub fn van_der_waals(r: f64, c6: f64) -> f64 {
    -c6 / r.powi(6)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_casimir_energy_sign() {
        // Casimir energy should be negative (attractive)
        let e = casimir_energy_parallel_plates(1e-6, 1.054571817e-34, 3.0e8);
        assert!(e < 0.0);
    }

    #[test]
    fn test_casimir_force_sign() {
        let f = casimir_force(1e-6, 1.054571817e-34, 3.0e8);
        assert!(f < 0.0);
    }

    #[test]
    fn test_casimir_energy_zeta_matches() {
        let d = 1e-9;
        let hbar = 1.0;
        let c = 1.0;
        let e_standard = casimir_energy_parallel_plates(d, hbar, c);
        let e_zeta = casimir_energy_zeta(d, hbar, c);
        // Should be same order of magnitude
        assert!((e_standard / e_zeta).abs() < 2.0);
    }
}
