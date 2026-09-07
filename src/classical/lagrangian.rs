//! Lagrangian mechanics and equations of motion.

/// Lagrangian: L = T - V.
pub fn lagrangian_l(
    velocity: &[f64],
    position: &[f64],
    mass: f64,
    potential: &dyn Fn(&[f64]) -> f64,
) -> f64 {
    kinetic_energy(velocity, mass) - potential(position)
}

/// Kinetic energy: T = ½ m Σ v_i².
pub fn kinetic_energy(velocity: &[f64], mass: f64) -> f64 {
    0.5 * mass * velocity.iter().map(|v| v * v).sum::<f64>()
}

/// Potential energy V(q).
pub fn potential_energy(position: &[f64], potential: &dyn Fn(&[f64]) -> f64) -> f64 {
    potential(position)
}

/// Euler-Lagrange equation: d/dt (∂L/∂q̇_i) - ∂L/∂q_i = 0.
///
/// Returns the acceleration q̈_i for each degree of freedom.
pub fn euler_lagrange(
    position: &[f64],
    velocity: &[f64],
    mass: f64,
    potential: &dyn Fn(&[f64]) -> f64,
    h: f64,
) -> Vec<f64> {
    let n = position.len();
    let mut acceleration = vec![0.0; n];

    for i in 0..n {
        // d/dt (∂L/∂v_i) = m a_i (for constant mass)
        // -∂L/∂q_i = -∂V/∂q_i = F_i
        let mut q_plus = position.to_vec();
        let mut q_minus = position.to_vec();
        q_plus[i] += h;
        q_minus[i] -= h;
        let force = -(potential(&q_plus) - potential(&q_minus)) / (2.0 * h);
        acceleration[i] = force / mass;
    }

    acceleration
}

/// Lagrangian for a charged particle: L = ½mv² + qv·A - qφ.
pub fn lagrangian_charged(
    velocity: &[f64],
    position: &[f64],
    mass: f64,
    charge: f64,
    a_field: &dyn Fn(&[f64]) -> Vec<f64>,
    scalar_potential: &dyn Fn(&[f64]) -> f64,
) -> f64 {
    let t = kinetic_energy(velocity, mass);
    let a = a_field(position);
    let v_a: f64 = velocity.iter().zip(a.iter()).map(|(v, a)| v * a).sum();
    let phi = scalar_potential(position);
    t + charge * v_a - charge * phi
}

/// Routhian: R(q, q̇, p) = Σ p_i q̇_i - L.
pub fn routhian(
    velocity: &[f64],
    momentum: &[f64],
    position: &[f64],
    mass: f64,
    potential: &dyn Fn(&[f64]) -> f64,
) -> f64 {
    let pv: f64 = momentum.iter().zip(velocity.iter()).map(|(p, v)| p * v).sum();
    let l = lagrangian_l(velocity, position, mass, potential);
    pv - l
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_harmonic_oscillator() {
        let v = |q: &[f64]| 0.5 * q[0] * q[0]; // V = ½x²
        let l = lagrangian_l(&[1.0], &[0.5], 1.0, &v);
        assert!((l - (0.5 - 0.125)).abs() < 1e-10);
    }
}
