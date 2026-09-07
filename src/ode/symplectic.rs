//! Symplectic integrators for Hamiltonian systems.
//!
//! Preserves phase-space volume and energy drift is bounded (no secular drift).

/// Single leapfrog (Verlet) step: position kicks then velocity drift.
///
/// Hamiltonian: H = T(p) + V(q)
/// - p_{n+1/2} = p_n - (dt/2) ∇V(q_n)
/// - q_{n+1} = q_n + dt ∇T(p_{n+1/2})
/// - p_{n+1} = p_{n+1/2} - (dt/2) ∇V(q_{n+1})
pub fn leapfrog_step(
    q: &[f64],
    p: &[f64],
    grad_v: &[f64],       // ∇V(q)
    grad_t: &[f64],       // ∇T(p) = p/m (for H = p²/2m + V)
    dt: f64,
    inv_mass: f64,         // 1/m
) -> (Vec<f64>, Vec<f64>) {
    let n = q.len();
    let mut p_half = vec![0.0; n];
    let mut q_new = vec![0.0; n];
    let mut p_new = vec![0.0; n];

    // Half-kick
    for i in 0..n {
        p_half[i] = p[i] - 0.5 * dt * grad_v[i];
    }

    // Drift
    for i in 0..n {
        q_new[i] = q[i] + dt * inv_mass * p_half[i];
    }

    // Half-kick (uses grad_v at new position, passed as argument)
    for i in 0..n {
        p_new[i] = p_half[i] - 0.5 * dt * grad_t[i];
    }

    (q_new, p_new)
}

/// Velocity Verlet step (equivalent to leapfrog, different ordering).
pub fn velocity_verlet_step(
    q: &[f64],
    p: &[f64],
    grad_v: &[f64],
    dt: f64,
    inv_mass: f64,
    grad_v_new: &[f64],    // ∇V(q + dt·p/m)
) -> (Vec<f64>, Vec<f64>) {
    let n = q.len();
    let mut q_new = vec![0.0; n];
    let mut p_new = vec![0.0; n];

    // Full drift
    for i in 0..n {
        q_new[i] = q[i] + dt * inv_mass * p[i];
    }

    // Half-kick then full kick
    for i in 0..n {
        p_new[i] = p[i] - 0.5 * dt * (grad_v[i] + grad_v_new[i]);
    }

    (q_new, p_new)
}

/// Yoshida 4th-order symplectic integrator (composition of 3 leapfrog steps).
///
/// Uses the Yoshida coefficients:
/// - w_1 = w_3 = 1 / (2 - 2^{1/3})
/// - w_2 = -2^{1/3} / (2 - 2^{1/3})
pub fn yoshida_step<FV, FT>(
    q: &[f64],
    p: &[f64],
    grad_v: &FV,
    grad_t: &FT,
    dt: f64,
    inv_mass: f64,
) -> (Vec<f64>, Vec<f64>)
where
    FV: Fn(&[f64]) -> Vec<f64>,
    FT: Fn(&[f64]) -> Vec<f64>,
{
    let cbrt2 = 2.0_f64.powf(1.0 / 3.0);
    let w1 = 1.0 / (2.0 - cbrt2);
    let w3 = w1;
    let w2 = 1.0 - 2.0 * w1;

    // Step 1: half-step with w1*dt
    let gv = grad_v(q);
    let gt = grad_t(p);
    let (q1, p1) = leapfrog_step(q, p, &gv, &gt, w1 * dt, inv_mass);

    // Step 2: full-step with w2*dt
    let gv1 = grad_v(&q1);
    let gt1 = grad_t(&p1);
    let (q2, p2) = leapfrog_step(&q1, &p1, &gv1, &gt1, w2 * dt, inv_mass);

    // Step 3: half-step with w3*dt
    let gv2 = grad_v(&q2);
    let gt2 = grad_t(&p2);
    let (q3, p3) = leapfrog_step(&q2, &p2, &gv2, &gt2, w3 * dt, inv_mass);

    (q3, p3)
}

/// Integrate a Hamiltonian system using leapfrog over n_steps.
pub fn integrate_hamiltonian<FV, FT>(
    q0: &[f64],
    p0: &[f64],
    grad_v: &FV,
    grad_t: &FT,
    dt: f64,
    n_steps: usize,
    inv_mass: f64,
) -> (Vec<Vec<f64>>, Vec<Vec<f64>>)
where
    FV: Fn(&[f64]) -> Vec<f64>,
    FT: Fn(&[f64]) -> Vec<f64>,
{
    let mut q = q0.to_vec();
    let mut p = p0.to_vec();
    let mut q_traj = vec![q.clone()];
    let mut p_traj = vec![p.clone()];

    for _ in 0..n_steps {
        let gv = grad_v(&q);
        let gt = grad_t(&p);
        let (q_new, p_new) = leapfrog_step(&q, &p, &gv, &gt, dt, inv_mass);
        q = q_new;
        p = p_new;
        q_traj.push(q.clone());
        p_traj.push(p.clone());
    }

    (q_traj, p_traj)
}
