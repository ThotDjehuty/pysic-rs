//! Lie algebras and gauge group structures.

/// Structure constants f^{abc} for SU(2): f^{abc} = ε^{abc}.
pub fn su2_structure_constants() -> [[[f64; 3]; 3]; 3] {
    let mut f = [[[0.0; 3]; 3]; 3];
    f[0][1][2] = 1.0; f[1][2][0] = 1.0; f[2][0][1] = 1.0;
    f[0][2][1] = -1.0; f[2][1][0] = -1.0; f[1][0][2] = -1.0;
    f
}

/// Structure constants f^{abc} for SU(3) (Gell-Mann basis).
pub fn su3_structure_constants() -> [[[f64; 8]; 8]; 8] {
    let mut f = [[[0.0; 8]; 8]; 8];

    // Non-zero structure constants for SU(3)
    // f^{123} = 1
    f[0][1][2] = 1.0; f[1][2][0] = 1.0; f[2][0][1] = 1.0;
    f[0][2][1] = -1.0; f[2][1][0] = -1.0; f[1][0][2] = -1.0;

    // f^{147} = f^{246} = f^{257} = f^{345} = 1/2
    f[0][3][6] = 0.5; f[0][6][3] = 0.5; f[6][3][0] = 0.5;
    f[0][6][3] = -0.5; f[6][0][3] = -0.5; f[3][0][6] = -0.5;
    f[3][6][0] = -0.5; f[6][0][3] = 0.5; f[0][3][6] = -0.5;

    // More non-zero components...
    f[1][3][5] = 0.5; f[1][5][3] = 0.5; f[5][3][1] = 0.5;
    f[1][4][6] = 0.5; f[1][6][4] = 0.5; f[6][4][1] = 0.5;
    f[2][3][4] = 0.5; f[2][4][3] = 0.5; f[4][3][2] = 0.5;

    // f^{458} = f^{678} = √3/2
    let sqrt3_2 = 3.0_f64.sqrt() / 2.0;
    f[3][4][7] = sqrt3_2; f[4][7][3] = sqrt3_2; f[7][3][4] = sqrt3_2;
    f[5][6][7] = sqrt3_2; f[6][7][5] = sqrt3_2; f[7][5][6] = sqrt3_2;

    f
}

/// Get structure constants for a given algebra.
pub fn structure_constants(algebra: &str) -> Vec<Vec<Vec<f64>>> {
    match algebra {
        "su2" => {
            let f = su2_structure_constants();
            f.iter().map(|a| a.iter().map(|b| b.to_vec()).collect()).collect()
        }
        "su3" => {
            let f = su3_structure_constants();
            f.iter().map(|a| a.iter().map(|b| b.to_vec()).collect()).collect()
        }
        _ => vec![],
    }
}

/// Gell-Mann matrices λ_a (a=1..8) for SU(3) fundamental representation.
pub fn gell_mann_matrices() -> Vec<Vec<num_complex::Complex64>> {
    use num_complex::Complex64;
    let zero = Complex64::new(0.0, 0.0);
    let one = Complex64::new(1.0, 0.0);

    vec![
        // λ1
        vec![zero, one, zero, one, zero, zero, zero, zero, zero],
        // λ2
        vec![zero, Complex64::new(0.0, -1.0), zero, Complex64::new(0.0, 1.0), zero, zero, zero, zero, zero],
        // λ3
        vec![one, zero, zero, zero, -one, zero, zero, zero, zero],
        // λ4
        vec![zero, zero, one, zero, zero, zero, one, zero, zero],
        // λ5
        vec![zero, zero, Complex64::new(0.0, -1.0), zero, zero, zero, Complex64::new(0.0, 1.0), zero, zero],
        // λ6
        vec![zero, zero, zero, zero, zero, one, zero, zero, one],
        // λ7
        vec![zero, zero, zero, zero, zero, Complex64::new(0.0, -1.0), zero, zero, Complex64::new(0.0, 1.0)],
        // λ8
        {
            let inv_sqrt3 = 1.0 / 3.0_f64.sqrt();
            vec![
                Complex64::new(inv_sqrt3, 0.0), zero, zero,
                zero, Complex64::new(inv_sqrt3, 0.0), zero,
                zero, zero, Complex64::new(-2.0 * inv_sqrt3, 0.0),
            ]
        },
    ]
}

/// SU(2) generators T_a = σ_a / 2 (fundamental representation).
pub fn su2_generators() -> [[[f64; 2]; 2]; 3] {
    [
        [[0.0, 0.5], [0.5, 0.0]],
        [[0.0, -0.5], [0.5, 0.0]],
        [[0.5, 0.0], [0.0, -0.5]],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_su2_jacobi() {
        let f = su2_structure_constants();
        // Jacobi identity: Σ_d (f^{ade} f^{bcd} + f^{bde} f^{cad} + f^{cde} f^{abd}) = 0
        for a in 0..3 {
            for b in 0..3 {
                for c in 0..3 {
                    let mut sum = 0.0;
                    for d in 0..3 {
                        for e in 0..3 {
                            sum += f[a][d][e] * f[b][c][d];
                            sum += f[b][d][e] * f[c][a][d];
                            sum += f[c][d][e] * f[a][b][d];
                        }
                    }
                    assert!(sum.abs() < 1e-10,
                        "Jacobi failed for a={}, b={}, c={}: sum={}", a, b, c, sum);
                }
            }
        }
    }
}
