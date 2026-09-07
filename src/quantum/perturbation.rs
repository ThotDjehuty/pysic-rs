//! Rayleigh-Schrödinger perturbation theory.

/// First-order energy correction: E_n^(1) = ⟨n|V|n⟩.
pub fn rayleigh_schrodinger_1st(
    e0: &[f64],           // Unperturbed energies
    v_matrix: &[Vec<f64>], // Perturbation matrix
    n: usize,             // State index
) -> f64 {
    v_matrix[n][n]
}

/// Second-order energy correction:
/// E_n^(2) = Σ_{m≠n} |⟨m|V|n⟩|² / (E_n^(0) - E_m^(0))
pub fn rayleigh_schrodinger_2nd(
    e0: &[f64],
    v_matrix: &[Vec<f64>],
    n: usize,
) -> f64 {
    let mut sum = 0.0;
    for m in 0..e0.len() {
        if m == n {
            continue;
        }
        let denom = e0[n] - e0[m];
        if denom.abs() > 1e-15 {
            sum += v_matrix[m][n] * v_matrix[m][n] / denom;
        }
    }
    sum
}

/// First-order wavefunction correction:
/// |n^(1)⟩ = Σ_{m≠n} ⟨m|V|n⟩ / (E_n^(0) - E_m^(0)) |m^(0)⟩
pub fn perturbed_wavefunction_1st(
    e0: &[f64],
    v_matrix: &[Vec<f64>],
    states: &[Vec<f64>],
    n: usize,
) -> Vec<f64> {
    let dim = e0.len();
    let state_dim = states[0].len();
    let mut psi1 = vec![0.0; state_dim];

    for m in 0..dim {
        if m == n {
            continue;
        }
        let denom = e0[n] - e0[m];
        if denom.abs() > 1e-15 {
            let coeff = v_matrix[m][n] / denom;
            for i in 0..state_dim {
                psi1[i] += coeff * states[m][i];
            }
        }
    }

    psi1
}

/// Construct the perturbation matrix from a potential function.
/// V_{mn} = ∫ φ_m(x) V(x) φ_n(x) dx
pub fn perturbation_matrix(
    states: &[Vec<f64>],
    potential: &[f64],
    dx: f64,
) -> Vec<Vec<f64>> {
    let dim = states.len();
    let mut v_matrix = vec![vec![0.0; dim]; dim];

    for m in 0..dim {
        for n in 0..dim {
            let mut integral = 0.0;
            for i in 0..states[0].len() {
                integral += states[m][i] * potential[i] * states[n][i] * dx;
            }
            v_matrix[m][n] = integral;
        }
    }

    v_matrix
}
