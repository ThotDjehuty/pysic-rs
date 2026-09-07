//! Chern numbers and topological invariants.

use num_complex::Complex64;

/// First Chern number: c₁ = (1/2π) ∫ F_{12} d²k over the Brillouin zone.
pub fn chern_number(
    berry_curvature: &Vec<Vec<f64>>,  // F_{xy}(k_x, k_y) on a grid
    dk: f64,
) -> f64 {
    let mut integral = 0.0;
    for row in berry_curvature {
        for &val in row {
            integral += val * dk * dk;
        }
    }
    integral / (2.0 * std::f64::consts::PI)
}

/// Discrete Chern number using the lattice gauge theory formula.
///
/// c₁ = (1/2π) Σ_{plaquettes} arg(U_{12}(k))
/// where U_{12} = U_1(k) U_2(k+ê₁) U_1^{-1}(k+ê₂) U_2^{-1}(k)
/// is the Wilson loop around a plaquette.
pub fn chern_number_discrete(
    overlap_x: &[Vec<Complex64>],  // ⟨n(k)|n(k+ê₁)⟩
    overlap_y: &[Vec<Complex64>],  // ⟨n(k)|n(k+ê₂)⟩
) -> f64 {
    let nx = overlap_x.len();
    let ny = if nx > 0 { overlap_x[0].len() } else { 0 };
    let mut total = 0.0;

    for i in 0..nx {
        for j in 0..ny {
            let ip = (i + 1) % nx;
            let jp = (j + 1) % ny;

            // Wilson loop around plaquette (i,j)
            let u1 = overlap_x[i][j];
            let u2 = overlap_y[ip][j];
            let u3 = overlap_x[i][jp].conj();
            let u4 = overlap_y[i][j].conj();

            let wilson = u1 * u2 * u3 * u4;
            total += wilson.arg();
        }
    }

    total / (2.0 * std::f64::consts::PI)
}

/// Second Chern number (4D): c₂ = (1/32π²) ∫ ε^{μνρσ} F_{μν} F_{ρσ} d⁴k.
pub fn chern_number_4d(
    field_strength: &dyn Fn(&[f64; 4]) -> [[f64; 4]; 4],
    k_min: [f64; 4],
    k_max: [f64; 4],
    n_grid: usize,
) -> f64 {
    let dk: Vec<f64> = (0..4).map(|i| (k_max[i] - k_min[i]) / n_grid as f64).collect();
    let mut integral = 0.0;

    for i0 in 0..n_grid {
        for i1 in 0..n_grid {
            for i2 in 0..n_grid {
                for i3 in 0..n_grid {
                    let k = [
                        k_min[0] + (i0 as f64 + 0.5) * dk[0],
                        k_min[1] + (i1 as f64 + 0.5) * dk[1],
                        k_min[2] + (i2 as f64 + 0.5) * dk[2],
                        k_min[3] + (i3 as f64 + 0.5) * dk[3],
                    ];
                    let f = field_strength(&k);
                    let volume = dk.iter().product::<f64>();

                    // ε^{0123} F_{01} F_{23} + permutations
                    let mut integrand = 0.0;
                    integrand += f[0][1] * f[2][3] - f[0][2] * f[1][3] + f[0][3] * f[1][2];
                    integral += integrand * volume;
                }
            }
        }
    }

    integral / (32.0 * std::f64::consts::PI * std::f64::consts::PI)
}

/// TKNN formula: total Chern number for a multi-band system.
pub fn tknn_chern(
    bands: &[Vec<Vec<Complex64>>],  // Bloch states for each band
    k_points: &[[f64; 2]],
    dk: f64,
) -> Vec<f64> {
    bands.iter().map(|band| {
        // Simplified: compute Berry curvature numerically
        let mut total_curvature = 0.0;
        for &k in k_points {
            // Placeholder: in practice, compute F_{xy} at each k-point
            total_curvature += 0.0;
        }
        total_curvature * dk * dk / (2.0 * std::f64::consts::PI)
    }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chern_trivial() {
        // Zero Berry curvature → c₁ = 0
        let f = vec![vec![0.0; 10]; 10];
        let c = chern_number(&f, 0.1);
        assert!(c.abs() < 1e-10);
    }
}
