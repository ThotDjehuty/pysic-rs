//! Radiation formulas: Larmor, dipole, multipole expansion.

/// Larmor formula: P = (2/3) (q² a²)/(4π ε₀ c³)
/// Power radiated by an accelerating point charge.
pub fn larmor_formula(
    charge: f64,
    acceleration: f64,
    eps0: f64,
    c: f64,
) -> f64 {
    (2.0 / 3.0) * charge * charge * acceleration * acceleration
        / (4.0 * std::f64::consts::PI * eps0 * c * c * c)
}

/// Dipole radiation power: P = (μ₀ ω⁴ |p|²)/(12π c³)
/// where p is the electric dipole moment.
pub fn dipole_radiation(
    dipole_moment: f64,
    omega: f64,
    mu0: f64,
    c: f64,
) -> f64 {
    mu0 * omega.powi(4) * dipole_moment * dipole_moment
        / (12.0 * std::f64::consts::PI * c * c * c)
}

/// Radiation power for an oscillating dipole: p(t) = p₀ cos(ωt).
pub fn radiation_power_oscillating(
    p0: f64,
    omega: f64,
    mu0: f64,
    c: f64,
) -> f64 {
    dipole_radiation(p0, omega, mu0, c)
}

/// Angular distribution of dipole radiation:
/// dP/dΩ = (μ₀ ω⁴ |p|²)/(32π² c³) sin²θ
pub fn dipole_angular_distribution(
    theta: f64,
    dipole_moment: f64,
    omega: f64,
    mu0: f64,
    c: f64,
) -> f64 {
    let prefactor = mu0 * omega.powi(4) * dipole_moment * dipole_moment
        / (32.0 * std::f64::consts::PI * std::f64::consts::PI * c * c * c);
    prefactor * theta.sin().powi(2)
}

/// Total power radiated by a magnetic dipole:
/// P = (μ₀ m² ω⁴)/(12π c³)
pub fn magnetic_dipole_radiation(
    magnetic_moment: f64,
    omega: f64,
    mu0: f64,
    c: f64,
) -> f64 {
    mu0 * omega.powi(4) * magnetic_moment * magnetic_moment
        / (12.0 * std::f64::consts::PI * c * c * c)
}

/// Cyclotron radiation power: P = (2/3) r_e m c (γ⁴ ω_c² β²)
/// where r_e is the classical electron radius, ω_c is the cyclotron frequency.
pub fn cyclotron_radiation(
    gamma: f64,
    omega_c: f64,
    beta: f64,
    classical_radius: f64,
    mass: f64,
    c: f64,
) -> f64 {
    (2.0 / 3.0) * classical_radius * mass * c * gamma.powi(4) * omega_c * omega_c * beta * beta
}

/// Synchrotron power: P = (2/3) r_e c γ⁴ / ρ²
/// where ρ is the bending radius.
pub fn synchrotron_power(
    gamma: f64,
    rho: f64,
    classical_radius: f64,
    c: f64,
) -> f64 {
    (2.0 / 3.0) * classical_radius * c * gamma.powi(4) / (rho * rho)
}

/// Bremsstrahlung spectrum: dI/dω ∝ Z²/ω (for ω < ω_max).
pub fn bremsstrahlung_spectrum(
    omega: f64,
    omega_max: f64,
    z: f64,
) -> f64 {
    if omega > omega_max || omega < 1e-15 {
        0.0
    } else {
        z * z / omega
    }
}

/// Thomson scattering cross-section: σ_T = (8π/3) r_e².
pub fn thomson_cross_section(classical_radius: f64) -> f64 {
    (8.0 / 3.0) * std::f64::consts::PI * classical_radius * classical_radius
}

/// Compton scattering wavelength shift: Δλ = (h/(m_e c))(1 - cos θ).
pub fn compton_wavelength_shift(
    theta: f64,
    h: f64,
    m_e: f64,
    c: f64,
) -> f64 {
    (h / (m_e * c)) * (1.0 - theta.cos())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_thomson_cross_section() {
        let r_e = 2.8179403262e-15; // Classical electron radius
        let sigma = thomson_cross_section(r_e);
        assert!((sigma - 6.6524587321e-29).abs() < 1e-35);
    }
}
