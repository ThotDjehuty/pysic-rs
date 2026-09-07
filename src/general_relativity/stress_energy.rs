//! Stress-energy tensors for matter sources.

use ndarray::Array2;

/// Perfect fluid stress-energy tensor:
/// T_μν = (ρ + p) u_μ u_ν + p g_μν
pub fn perfect_fluid_stress_energy(
    rho: f64,
    pressure: f64,
    u: &[f64],  // 4-velocity
    g: &Array2<f64>,
) -> Array2<f64> {
    let n = u.len();
    let mut t = Array2::zeros((n, n));

    for mu in 0..n {
        for nu in 0..n {
            t[[mu, nu]] = (rho + pressure) * u[mu] * u[nu] + pressure * g[[mu, nu]];
        }
    }

    t
}

/// Electromagnetic stress-energy tensor:
/// T_μν = F_μρ F_ν^ρ - ¼ g_μν F_ρσ F^ρσ
pub fn electromagnetic_stress_energy(
    f_tensor: &Array2<f64>,
    g: &Array2<f64>,
    g_inv: &Array2<f64>,
) -> Array2<f64> {
    let n = f_tensor.nrows();
    let mut t = Array2::zeros((n, n));

    // F_ρσ F^ρσ
    let mut f_sq = 0.0;
    for rho in 0..n {
        for sigma in 0..n {
            let mut f_up = 0.0;
            for mu in 0..n {
                f_up += g_inv[[mu, sigma]] * f_tensor[[rho, mu]];
            }
            f_sq += f_tensor[[rho, sigma]] * f_up;
        }
    }

    for mu in 0..n {
        for nu in 0..n {
            // F_μρ F_ν^ρ
            let mut ff = 0.0;
            for rho in 0..n {
                let f_nu_rho = 0.0;
                for sig in 0..n {
                    // F_ν^ρ = g^{ρσ} F_νσ
                    let _ = sig;
                }
                ff += f_tensor[[mu, rho]] * g_inv[[rho, rho]] * f_tensor[[nu, rho]];
            }
            t[[mu, nu]] = ff - 0.25 * g[[mu, nu]] * f_sq;
        }
    }

    t
}

/// Cosmological constant contribution: T^Λ_μν = -Λ/(8πG) g_μν.
pub fn cosmological_stress_energy(
    lambda: f64,
    g: &Array2<f64>,
) -> Array2<f64> {
    let n = g.nrows();
    let coeff = -lambda / (8.0 * std::f64::consts::PI * crate::constants::G());
    let mut t = Array2::zeros((n, n));
    for i in 0..n {
        for j in 0..n {
            t[[i, j]] = coeff * g[[i, j]];
        }
    }
    t
}
