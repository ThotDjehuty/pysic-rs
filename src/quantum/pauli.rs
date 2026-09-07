//! Pauli matrices and spin operators.

use num_complex::Complex64;

/// Pauli matrix σ_x.
pub fn pauli_x() -> [[Complex64; 2]; 2] {
    [[Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0)],
     [Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)]]
}

/// Pauli matrix σ_y.
pub fn pauli_y() -> [[Complex64; 2]; 2] {
    [[Complex64::new(0.0, 0.0), Complex64::new(0.0, -1.0)],
     [Complex64::new(0.0, 1.0), Complex64::new(0.0, 0.0)]]
}

/// Pauli matrix σ_z.
pub fn pauli_z() -> [[Complex64; 2]; 2] {
    [[Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)],
     [Complex64::new(0.0, 0.0), Complex64::new(-1.0, 0.0)]]
}

/// Return all three Pauli matrices as a tuple (σ_x, σ_y, σ_z).
pub fn pauli_matrices() -> ([[Complex64; 2]; 2], [[Complex64; 2]; 2], [[Complex64; 2]; 2]) {
    (pauli_x(), pauli_y(), pauli_z())
}

/// Spin-1/2 operators S_i = (ℏ/2) σ_i.
pub fn spin_operators(hbar: f64) -> ([[Complex64; 2]; 2], [[Complex64; 2]; 2], [[Complex64; 2]; 2]) {
    let half = Complex64::new(hbar / 2.0, 0.0);
    let (sx, sy, sz) = pauli_matrices();
    (
        sx.map(|row| row.map(|c| c * half)),
        sy.map(|row| row.map(|c| c * half)),
        sz.map(|row| row.map(|c| c * half)),
    )
}

/// Tensor product (Kronecker product) of two matrices.
pub fn tensor_product(a: &[[Complex64; 2]; 2], b: &[[Complex64; 2]; 2]) -> [[Complex64; 4]; 4] {
    let mut result = [[Complex64::new(0.0, 0.0); 4]; 4];
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                for l in 0..2 {
                    result[2 * i + k][2 * j + l] = a[i][j] * b[k][l];
                }
            }
        }
    }
    result
}

/// Compute the density matrix ρ = |ψ⟩⟨ψ| from a state vector.
pub fn density_matrix(psi: &[Complex64]) -> Vec<Vec<Complex64>> {
    let n = psi.len();
    let mut rho = vec![vec![Complex64::new(0.0, 0.0); n]; n];
    for i in 0..n {
        for j in 0..n {
            rho[i][j] = psi[i] * psi[j].conj();
        }
    }
    rho
}

/// Partial trace over subsystem B (last qubits) of a density matrix.
/// Assumes the density matrix is for n_qubits qubits (dimension 2^n_qubits).
/// Tr_B traces out the last qubits_keep qubits.
pub fn partial_trace(rho: &[Vec<Complex64>], n_total_qubits: usize, qubits_keep: usize) -> Vec<Vec<Complex64>> {
    let dim_total = 1 << n_total_qubits;
    let dim_keep = 1 << qubits_keep;
    let dim_trace = 1 << (n_total_qubits - qubits_keep);

    let mut result = vec![vec![Complex64::new(0.0, 0.0); dim_keep]; dim_keep];

    for i in 0..dim_keep {
        for j in 0..dim_keep {
            let mut sum = Complex64::new(0.0, 0.0);
            for k in 0..dim_trace {
                let full_i = i * dim_trace + k;
                let full_j = j * dim_trace + k;
                if full_i < dim_total && full_j < dim_total {
                    sum += rho[full_i][full_j];
                }
            }
            result[i][j] = sum;
        }
    }

    result
}

/// Purity of a density matrix: Tr(ρ²).
pub fn purity(rho: &[Vec<Complex64>]) -> f64 {
    let n = rho.len();
    let mut trace_sq = Complex64::new(0.0, 0.0);
    for i in 0..n {
        for j in 0..n {
            trace_sq += rho[i][j] * rho[j][i];
        }
    }
    trace_sq.re
}

/// Von Neumann entropy: S = -Tr(ρ ln ρ).
/// Requires eigenvalue decomposition of ρ.
pub fn von_neumann_entropy(rho: &[Vec<Complex64>]) -> f64 {
    let n = rho.len();
    // Simple case: pure state
    let pur = purity(rho);
    if (pur - 1.0).abs() < 1e-10 {
        return 0.0;
    }
    // For mixed states, use Tr(ρ ln ρ) approximation via eigenvalues
    // This is a simplified version — production code would use full diagonalization
    let diag_sum: f64 = (0..n).map(|i| {
        let p = rho[i][i].re.max(1e-15);
        -p * p.ln()
    }).sum();
    diag_sum
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_pauli_algebra() {
        let (sx, sy, sz) = pauli_matrices();
        // σ_x² = I
        let sx2 = matrix_mul_2x2(&sx, &sx);
        assert_relative_eq!(sx2[0][0].re, 1.0, epsilon = 1e-10);
        assert_relative_eq!(sx2[1][1].re, 1.0, epsilon = 1e-10);

        // [σ_x, σ_y] = 2i σ_z
        let commutator = matrix_commutator(&sx, &sy);
        let two_i_sz = sz.map(|row| row.map(|c| c * Complex64::new(0.0, 2.0)));
        for i in 0..2 {
            for j in 0..2 {
                assert_relative_eq!(commutator[i][j].re, two_i_sz[i][j].re, epsilon = 1e-10);
                assert_relative_eq!(commutator[i][j].im, two_i_sz[i][j].im, epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn test_density_matrix_purity() {
        let psi = vec![Complex64::new(1.0, 0.0), Complex64::new(0.0, 0.0)];
        let rho = density_matrix(&psi);
        assert_relative_eq!(purity(&rho), 1.0, epsilon = 1e-10);
    }
}

fn matrix_mul_2x2(a: &[[Complex64; 2]; 2], b: &[[Complex64; 2]; 2]) -> [[Complex64; 2]; 2] {
    let mut c = [[Complex64::new(0.0, 0.0); 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

fn matrix_commutator(a: &[[Complex64; 2]; 2], b: &[[Complex64; 2]; 2]) -> [[Complex64; 2]; 2] {
    let ab = matrix_mul_2x2(a, b);
    let ba = matrix_mul_2x2(b, a);
    [[ab[0][0] - ba[0][0], ab[0][1] - ba[0][1]],
     [ab[1][0] - ba[1][0], ab[1][1] - ba[1][1]]]
}
