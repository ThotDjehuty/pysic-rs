//! Classical 4th-order Runge-Kutta integrator.

/// Solve dy/dt = f(t, y) using RK4.
///
/// # Arguments
/// * `f` - Right-hand side: f(t, y) -> dy/dt
/// * `y0` - Initial state vector
/// * `t_span` - (t_start, t_end)
/// * `n_steps` - Number of integration steps
///
/// # Returns
/// * `times` - Time grid
/// * `trajectory` - Solution at each time point
pub fn rk4_solve<F>(
    f: &F,
    y0: &[f64],
    t_span: (f64, f64),
    n_steps: usize,
) -> (Vec<f64>, Vec<Vec<f64>>)
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    let (t0, tf) = t_span;
    let dt = (tf - t0) / n_steps as f64;
    let n = y0.len();

    let mut times = Vec::with_capacity(n_steps + 1);
    let mut trajectory = Vec::with_capacity(n_steps + 1);
    let mut y = y0.to_vec();

    times.push(t0);
    trajectory.push(y.clone());

    for step in 0..n_steps {
        let t = t0 + step as f64 * dt;
        let k1 = f(t, &y);
        let y_k2: Vec<f64> = y.iter().zip(k1.iter()).map(|(yi, ki)| yi + 0.5 * dt * ki).collect();
        let k2 = f(t + 0.5 * dt, &y_k2);
        let y_k3: Vec<f64> = y.iter().zip(k2.iter()).map(|(yi, ki)| yi + 0.5 * dt * ki).collect();
        let k3 = f(t + 0.5 * dt, &y_k3);
        let y_k4: Vec<f64> = y.iter().zip(k3.iter()).map(|(yi, ki)| yi + dt * ki).collect();
        let k4 = f(t + dt, &y_k4);

        for i in 0..n {
            y[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }

        times.push(t0 + (step + 1) as f64 * dt);
        trajectory.push(y.clone());
    }

    (times, trajectory)
}
