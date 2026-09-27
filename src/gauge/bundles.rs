//! Lie algebras and gauge group structures.

/// Fill every permutation of (a,b,c) with the parity-correct sign.
///
/// The structure constants of a compact simple Lie algebra in an orthonormal
/// basis are **totally antisymmetric**, so one independent value determines all
/// six orderings: even permutations get `+v`, odd permutations get `-v`. Writing
/// them out by hand invites sign slips and duplicate assignments, so the table
/// below lists only the independent entries and this helper does the rest.
fn set_antisymmetric<const N: usize>(f: &mut [[[f64; N]; N]; N], a: usize, b: usize, c: usize, v: f64) {
    f[a][b][c] = v;
    f[b][c][a] = v;
    f[c][a][b] = v;
    f[a][c][b] = -v;
    f[c][b][a] = -v;
    f[b][a][c] = -v;
}

/// Structure constants f^{abc} for SU(2): f^{abc} = ε^{abc}.
pub fn su2_structure_constants() -> [[[f64; 3]; 3]; 3] {
    let mut f = [[[0.0; 3]; 3]; 3];
    set_antisymmetric(&mut f, 0, 1, 2, 1.0);
    f
}

/// Structure constants f^{abc} for SU(3) in the Gell-Mann basis.
///
/// Convention: T^a = λ^a/2 with [T^a, T^b] = i f^{abc} T^c. The independent
/// non-zero values (1-indexed, as they are usually tabulated) are
///
/// ```text
/// f^{123} = 1
/// f^{147} = 1/2    f^{156} = -1/2
/// f^{246} = 1/2    f^{257} =  1/2
/// f^{345} = 1/2    f^{367} = -1/2
/// f^{458} = √3/2   f^{678} =  √3/2
/// ```
///
/// with every other non-zero component fixed by total antisymmetry. That gives
/// 9 × 3! = 54 non-zero entries out of 8³ = 512.
pub fn su3_structure_constants() -> [[[f64; 8]; 8]; 8] {
    let mut f = [[[0.0; 8]; 8]; 8];
    let h = 0.5;
    let s = 3.0_f64.sqrt() / 2.0;

    // Indices below are 0-based: subtract one from the tabulated values above.
    set_antisymmetric(&mut f, 0, 1, 2, 1.0); // f^{123}
    set_antisymmetric(&mut f, 0, 3, 6, h); //  f^{147}
    set_antisymmetric(&mut f, 0, 4, 5, -h); // f^{156}
    set_antisymmetric(&mut f, 1, 3, 5, h); //  f^{246}
    set_antisymmetric(&mut f, 1, 4, 6, h); //  f^{257}
    set_antisymmetric(&mut f, 2, 3, 4, h); //  f^{345}
    set_antisymmetric(&mut f, 2, 5, 6, -h); // f^{367}
    set_antisymmetric(&mut f, 3, 4, 7, s); //  f^{458}
    set_antisymmetric(&mut f, 5, 6, 7, s); //  f^{678}

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
