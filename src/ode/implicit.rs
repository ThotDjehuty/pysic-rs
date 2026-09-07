//! Implicit ODE solvers: backward Euler and Crank-Nicolson.

/// Backward Euler: y_{n+1} = y_n + h f(t_{n+1}, y_{n+1}).
/// Uses fixed-point iteration (Picard) for nonlinear systems.
pub fn backward_euler<F>(
    f: &F,
    y0: &[f64],
    t_span: (f64, f64),
    n_steps: usize,
    max_iter: usize,
    tol: f64,
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
        let t_new = t0 + (step + 1) as f64 * dt;
        let mut y_new = y.clone();

        // Fixed-point iteration
        for _iter in 0..max_iter {
            let k = f(t_new, &y_new);
            let y_pred: Vec<f64> = y.iter().zip(k.iter()).map(|(yi, ki)| yi + dt * ki).collect();

            let err: f64 = y_pred.iter().zip(y_new.iter()).map(|(a, b)| (a - b).powi(2)).sum::<f64>().sqrt();
            y_new = y_pred;

            if err < tol {
                break;
            }
        }

        y = y_new;
        times.push(t_new);
        trajectory.push(y.clone());
    }

    (times, trajectory)
}

/// Crank-Nicolson for ODEs: y_{n+1} = y_n + h/2 [f(t_n, y_n) + f(t_{n+1}, y_{n+1})].
pub fn crank_nicolson_ode<F>(
    f: &F,
    y0: &[f64],
    t_span: (f64, f64),
    n_steps: usize,
    max_iter: usize,
    tol: f64,
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
        let t_n = t0 + step as f64 * dt;
        let t_new = t0 + (step + 1) as f64 * dt;
        let k_n = f(t_n, &y);
        let mut y_new = y.clone();

        for _iter in 0..max_iter {
            let k_new = f(t_new, &y_new);
            let y_pred: Vec<f64> = y.iter().zip(k_n.iter()).zip(k_new.iter()).map(|((yi, kni), knwi)| yi + 0.5 * dt * (kni + knwi)).collect();

            let err: f64 = y_pred.iter().zip(y_new.iter()).map(|(a, b)| (a - b).powi(2)).sum::<f64>().sqrt();
            y_new = y_pred;

            if err < tol {
                break;
            }
        }

        y = y_new;
        times.push(t_new);
        trajectory.push(y.clone());
    }

    (times, trajectory)
}
