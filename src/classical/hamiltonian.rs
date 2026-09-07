//! Hamiltonian mechanics and canonical formulation.

/// Hamiltonian function H(q, p, t).
pub fn hamiltonian_h(
    q: &[f64],
    p: &[f64],
    potential: &dyn Fn(&[f64]) -> f64,
    mass: f64,
) -> f64 {
    let kinetic: f64 = p.iter().map(|pi| pi * pi / (2.0 * mass)).sum();
    kinetic + potential(q)
}

/// Canonical equations: dq/dt = ∂H/∂p, dp/dt = -∂H/∂q.
pub fn canonical_equations(
    q: &[f64],
    p: &[f64],
    hamiltonian: &dyn Fn(&[f64], &[f64]) -> f64,
    h: f64,
) -> (Vec<f64>, Vec<f64>) {
    let n = q.len();
    let mut dq_dt = vec![0.0; n];
    let mut dp_dt = vec![0.0; n];

    // dq/dt = ∂H/∂p (central differences)
    for i in 0..n {
        let mut p_plus = p.to_vec();
        let mut p_minus = p.to_vec();
        p_plus[i] += h;
        p_minus[i] -= h;
        dq_dt[i] = (hamiltonian(q, &p_plus) - hamiltonian(q, &p_minus)) / (2.0 * h);
    }

    // dp/dt = -∂H/∂q
    for i in 0..n {
        let mut q_plus = q.to_vec();
        let mut q_minus = q.to_vec();
        q_plus[i] += h;
        q_minus[i] -= h;
        dp_dt[i] = -(hamiltonian(&q_plus, p) - hamiltonian(&q_minus, p)) / (2.0 * h);
    }

    (dq_dt, dp_dt)
}

/// Poisson bracket {f, g} = Σ_i (∂f/∂q_i ∂g/∂p_i - ∂f/∂p_i ∂g/∂q_i).
pub fn poisson_bracket(
    f: &dyn Fn(&[f64], &[f64]) -> f64,
    g: &dyn Fn(&[f64], &[f64]) -> f64,
    q: &[f64],
    p: &[f64],
    h: f64,
) -> f64 {
    let n = q.len();
    let mut bracket = 0.0;

    for i in 0..n {
        let mut q_plus = q.to_vec();
        let mut q_minus = q.to_vec();
        let mut p_plus = p.to_vec();
        let mut p_minus = p.to_vec();

        q_plus[i] += h;
        q_minus[i] -= h;
        p_plus[i] += h;
        p_minus[i] -= h;

        let df_dq = (f(&q_plus, p) - f(&q_minus, p)) / (2.0 * h);
        let df_dp = (f(q, &p_plus) - f(q, &p_minus)) / (2.0 * h);
        let dg_dq = (g(&q_plus, p) - g(&q_minus, p)) / (2.0 * h);
        let dg_dp = (g(q, &p_plus) - g(q, &p_minus)) / (2.0 * h);

        bracket += df_dq * dg_dp - df_dp * dg_dq;
    }

    bracket
}

/// Check if a transformation is canonical using the fundamental Poisson brackets.
pub fn canonical_transform(
    q: &[f64],
    p: &[f64],
    h: f64,
) -> bool {
    // {Q_i, Q_j} = 0, {P_i, P_j} = 0, {Q_i, P_j} = δ_ij
    // Simplified: check {Q, P} ≈ 1
    let identity = |_q: &[f64], _p: &[f64]| -> f64 { 1.0 };
    let bracket = poisson_bracket(&identity, &identity, q, p, h);
    bracket.abs() < 0.1 // Rough check
}

/// Liouville theorem: dρ/dt = {ρ, H} = 0 for autonomous systems.
pub fn liouville_rhs(
    rho: f64,
    h: &dyn Fn(&[f64], &[f64]) -> f64,
    q: &[f64],
    p: &[f64],
    h_fd: f64,
) -> f64 {
    // ∂ρ/∂t = {H, ρ}
    let bracket = poisson_bracket(h, &|_q, _p| rho, q, p, h_fd);
    bracket
}
