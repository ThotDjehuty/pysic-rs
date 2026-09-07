//! Field strength tensor F_μν = ∂_μ A_ν - ∂_ν A_μ + ig [A_μ, A_ν].

/// Compute the field strength tensor F_μν from gauge field components.
///
/// For an abelian theory (QED): F_μν = ∂_μ A_ν - ∂_ν A_μ
pub fn field_strength(
    a_field: &[Vec<f64>],  // A_μ at each grid point: [n_points][n_components]
    dx: f64,
    n_points: usize,
    n_components: usize,
) -> Vec<Vec<Vec<f64>>> {
    let mut f = vec![vec![vec![0.0; n_components]; n_components]; n_points];

    for p in 0..n_points {
        for mu in 0..n_components {
            for nu in 0..n_components {
                if mu == nu {
                    continue;
                }
                // ∂_μ A_ν - ∂_ν A_μ
                let idx_plus_mu = (p + n_components).min(n_points - 1);
                let idx_minus_mu = p.saturating_sub(n_components);
                let idx_plus_nu = (p + 1).min(n_points - 1);
                let idx_minus_nu = p.saturating_sub(1);

                let dmu_anu = (a_field[idx_plus_mu][nu] - a_field[idx_minus_mu][nu]) / (2.0 * dx);
                let dnu_amu = (a_field[idx_plus_nu][mu] - a_field[idx_minus_nu][mu]) / (2.0 * dx);

                f[p][mu][nu] = dmu_anu - dnu_amu;
            }
        }
    }

    f
}

/// Compute F_μν for a point charge at the origin (Coulomb field).
/// E_r = Q/(4π ε₀ r²), B = 0.
pub fn field_strength_point_charge(
    charge: f64,
    r: f64,
    eps0: f64,
) -> f64 {
    charge / (4.0 * std::f64::consts::PI * eps0 * r * r)
}

/// Electromagnetic duality: F_μν → *F_μν = ½ ε_μνρσ F^ρσ.
pub fn dual_field_strength(
    f_tensor: &Vec<Vec<f64>>,
) -> Vec<Vec<f64>> {
    let n = f_tensor.len();
    let mut f_dual = vec![vec![0.0; n]; n];

    // Levi-Civita ε_{0123} = +1 in Minkowski
    let epsilon = [[[0.0; 4]; 4]; 4];
    // Simplified: for 4D, dual of F_μν
    f_dual[0][1] = f_tensor[2][3];
    f_dual[0][2] = -f_tensor[1][3];
    f_dual[0][3] = f_tensor[1][2];
    f_dual[1][2] = f_tensor[0][3];
    f_dual[1][3] = -f_tensor[0][2];
    f_dual[2][3] = f_tensor[0][1];

    // Antisymmetric
    for i in 0..n {
        for j in i + 1..n {
            f_dual[j][i] = -f_dual[i][j];
        }
    }

    f_dual
}
