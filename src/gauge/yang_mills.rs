//! Yang-Mills equations of motion and action.

/// Yang-Mills action: S = -¼ ∫ F^a_μν F^{a μν} d⁴x.
///
/// Simplified: computes the action density at a point.
pub fn yang_mills_action(
    f_tensor: &[Vec<f64>],  // F_μν at each point
    g_inv: &[Vec<f64>],     // Inverse metric
    n_components: usize,
    dx: f64,
) -> f64 {
    let mut action = 0.0;

    for mu in 0..n_components {
        for nu in 0..n_components {
            for rho in 0..n_components {
                for sigma in 0..n_components {
                    action += g_inv[mu][rho] * g_inv[nu][sigma]
                        * f_tensor[mu * n_components + nu][rho * n_components + sigma];
                }
            }
        }
    }

    -0.25 * action * dx.powi(4)
}

/// Yang-Mills equation of motion: D_μ F^μν = J^ν.
/// Returns the residual (should vanish on-shell).
pub fn yang_mills_eom(
    f_tensor: &[Vec<f64>],
    a_field: &[Vec<f64>],
    coupling: f64,
    dx: f64,
    n_components: usize,
) -> Vec<f64> {
    let mut eom = vec![0.0; n_components];

    for nu in 0..n_components {
        let mut sum = 0.0;
        for mu in 0..n_components {
            // ∂_μ F^μν
            let idx_plus = (mu + n_components).min(n_components * n_components - 1);
            let idx_minus = mu.saturating_sub(n_components);
            let df = (f_tensor[idx_plus * n_components + nu][mu * n_components + nu]
                - f_tensor[idx_minus * n_components + nu][mu * n_components + nu]) / (2.0 * dx);

            // ig [A_μ, F^μν] (commutator, simplified for abelian)
            let comm = coupling * (a_field[mu][0] * f_tensor[mu * n_components + nu][0 * n_components + nu]
                - f_tensor[mu * n_components + nu][0 * n_components + nu] * a_field[mu][0]);

            sum += df + comm;
        }
        eom[nu] = sum;
    }

    eom
}

/// Instanton action for SU(2): S_inst = 8π²/g².
pub fn instanton_action(g: f64) -> f64 {
    8.0 * std::f64::consts::PI * std::f64::consts::PI / (g * g)
}

/// 't Hooft symbol σ^a_μν for instanton construction.
pub fn thooft_symbol(a: usize, mu: usize, nu: usize) -> f64 {
    match (a, mu, nu) {
        (0, 1, 2) | (0, 2, 1) => if mu == 1 { 1.0 } else { -1.0 },
        (1, 0, 2) | (1, 2, 0) => if mu == 0 { -1.0 } else { 1.0 },
        (2, 0, 1) | (2, 1, 0) => if mu == 0 { 1.0 } else { -1.0 },
        // Self-dual: η̄ = 0 for anti-instanton
        _ => 0.0,
    }
}
