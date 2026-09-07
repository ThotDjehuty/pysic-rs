//! Hilbert space operations: inner products, projectors, coherent states.

use num_complex::Complex64;

/// Inner product ⟨φ|ψ⟩.
pub fn inner_product(phi: &[Complex64], psi: &[Complex64]) -> Complex64 {
    phi.iter().zip(psi.iter()).map(|(p, q)| p.conj() * q).sum()
}

/// Outer product |ψ⟩⟨φ|.
pub fn outer_product(psi: &[Complex64], phi: &[Complex64]) -> Vec<Vec<Complex64>> {
    let n = psi.len();
    let m = phi.len();
    (0..n).map(|i| (0..m).map(|j| psi[i] * phi[j].conj()).collect()).collect()
}

/// Projection operator onto state |ψ⟩: P = |ψ⟩⟨ψ|.
pub fn projector(psi: &[Complex64]) -> Vec<Vec<Complex64>> {
    outer_product(psi, psi)
}

/// Coherent state |α⟩ in the Fock basis (first n_terms terms).
///
/// |α⟩ = e^{-|α|²/2} Σ_n αⁿ/√(n!) |n⟩
pub fn coherent_state(alpha: Complex64, n_terms: usize) -> Vec<Complex64> {
    let norm = (-alpha.norm_sqr() / 2.0).exp();
    let mut state = Vec::with_capacity(n_terms);
    let mut alpha_pow = Complex64::new(1.0, 0.0);
    let mut factorial: f64 = 1.0;

    for n in 0..n_terms {
        state.push(norm * alpha_pow / factorial.sqrt());
        alpha_pow *= alpha;
        factorial *= (n + 1) as f64;
    }

    state
}

/// Number state |n⟩ (Fock state).
pub fn number_state(n: usize, dim: usize) -> Vec<Complex64> {
    let mut state = vec![Complex64::new(0.0, 0.0); dim];
    if n < dim {
        state[n] = Complex64::new(1.0, 0.0);
    }
    state
}

/// Squeeze operator S(ζ) applied to vacuum |0⟩.
///
/// S(ζ) = exp(ζ* a†²/2 - ζ a²/2)
/// For vacuum input, produces squeezed vacuum.
pub fn squeeze_state(zeta: Complex64, n_terms: usize) -> Vec<Complex64> {
    let r = zeta.norm();
    let theta = zeta.arg();
    let mut state = Vec::with_capacity(n_terms);

    for n in 0..n_terms {
        if n % 2 == 0 {
            let m = n / 2;
            let sech_r = 1.0 / r.cosh();
            let coeff = sech_r * (-Complex64::new(0.0, 1.0) * theta * m as f64).exp()
                * (r.tanh()).powi(m as i32)
                * (factorial_approx(2 * m) as f64).sqrt()
                / (factorial_approx(m) as f64);
            state.push(coeff);
        } else {
            state.push(Complex64::new(0.0, 0.0));
        }
    }

    // Normalize
    let norm: f64 = state.iter().map(|c| c.norm_sqr()).sum::<f64>().sqrt();
    if norm > 1e-15 {
        state.iter_mut().for_each(|c| *c /= norm);
    }

    state
}

fn factorial_approx(n: usize) -> usize {
    match n {
        0 | 1 => 1,
        _ => n * factorial_approx(n - 1),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_inner_product() {
        let psi = vec![Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)];
        let ip = inner_product(&psi, &psi);
        assert!((ip.norm() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_coherent_state_normalization() {
        let alpha = Complex64::new(1.0, 0.5);
        let state = coherent_state(alpha, 20);
        let norm: f64 = state.iter().map(|c| c.norm_sqr()).sum();
        assert!((norm - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_number_state() {
        let state = number_state(2, 5);
        assert!((state[2].re - 1.0).abs() < 1e-10);
        assert!((state[0].re).abs() < 1e-10);
    }
}
