//! Dormand-Prince adaptive RK45 ODE solver.

/// Solution output from adaptive RK45.
pub struct OdeSolution {
    pub times: Vec<f64>,
    pub trajectory: Vec<Vec<f64>>,
    pub n_steps: usize,
    pub n_evals: usize,
}

/// Solve dy/dt = f(t, y) using adaptive Dormand-Prince RK45.
///
/// # Arguments
/// * `f` - Right-hand side: f(t, y) -> dy/dt
/// * `y0` - Initial state vector
/// * `t_span` - (t_start, t_end)
/// * `rtol` - Relative tolerance (default 1e-6)
/// * `atol` - Absolute tolerance (default 1e-9)
/// * `max_steps` - Maximum number of steps
pub fn rk45_solve<F>(
    f: &F,
    y0: &[f64],
    t_span: (f64, f64),
    rtol: f64,
    atol: f64,
    max_steps: usize,
) -> OdeSolution
where
    F: Fn(f64, &[f64]) -> Vec<f64>,
{
    let (t0, tf) = t_span;
    let n = y0.len();
    let mut t = t0;
    let mut y = y0.to_vec();
    let mut dt = (tf - t0) / 100.0; // Initial step size guess

    let mut times = vec![t];
    let mut trajectory = vec![y.clone()];
    let mut n_evals = 0;

    // Dormand-Prince coefficients
    let a21 = 1.0/5.0;
    let a31 = 3.0/40.0; let a32 = 9.0/40.0;
    let a41 = 44.0/45.0; let a42 = -56.0/15.0; let a43 = 32.0/9.0;
    let a51 = 19372.0/6561.0; let a52 = -25360.0/2187.0; let a53 = 64448.0/6561.0; let a54 = -212.0/729.0;
    let a61 = 9017.0/3168.0; let a62 = -355.0/33.0; let a63 = 46732.0/5247.0; let a64 = 49.0/176.0; let a65 = -5103.0/18656.0;

    let b1 = 35.0/384.0; let b3 = 500.0/1113.0; let b4 = 125.0/192.0; let b5 = -2187.0/6784.0; let b6 = 11.0/84.0;

    let c1 = 0.0; let c2 = 1.0/5.0; let c3 = 3.0/10.0; let c4 = 4.0/5.0; let c5 = 8.0/9.0;

    for _step in 0..max_steps {
        if t >= tf {
            break;
        }
        if t + dt > tf {
            dt = tf - t;
        }

        let k1 = f(t, &y);
        let y2: Vec<f64> = y.iter().zip(k1.iter()).map(|(yi, ki)| yi + a21 * dt * ki).collect();
        let k2 = f(t + c2 * dt, &y2);
        let y3: Vec<f64> = y.iter().zip(k1.iter()).zip(k2.iter()).map(|((yi, ki), k2i)| yi + dt * (a31 * ki + a32 * k2i)).collect();
        let k3 = f(t + c3 * dt, &y3);
        let y4: Vec<f64> = y.iter().zip(k1.iter()).zip(k2.iter()).zip(k3.iter()).map(|(((yi, ki), k2i), k3i)| yi + dt * (a41 * ki + a42 * k2i + a43 * k3i)).collect();
        let k4 = f(t + c4 * dt, &y4);
        let y5: Vec<f64> = y.iter().zip(k1.iter()).zip(k2.iter()).zip(k3.iter()).zip(k4.iter()).map(|((((yi, ki), k2i), k3i), k4i)| yi + dt * (a51 * ki + a52 * k2i + a53 * k3i + a54 * k4i)).collect();
        let k5 = f(t + c5 * dt, &y5);
        let y6: Vec<f64> = y.iter().zip(k1.iter()).zip(k2.iter()).zip(k3.iter()).zip(k4.iter()).zip(k5.iter()).map(|(((((yi, ki), k2i), k3i), k4i), k5i)| yi + dt * (a61 * ki + a62 * k2i + a63 * k3i + a64 * k4i + a65 * k5i)).collect();
        let k6 = f(t + dt, &y6);

        n_evals += 6;

        // 5th order solution (for error estimation)
        let y_new: Vec<f64> = y.iter().zip(k1.iter()).zip(k3.iter()).zip(k4.iter()).zip(k5.iter()).zip(k6.iter()).map(|(((((yi, ki), k3i), k4i), k5i), k6i)| yi + dt * (b1 * ki + b3 * k3i + b4 * k4i + b5 * k5i + b6 * k6i)).collect();

        // Error estimate (difference between 4th and 5th order)
        let err: f64 = y_new.iter().zip(y.iter()).zip(k1.iter()).zip(k3.iter()).zip(k4.iter()).zip(k5.iter()).zip(k6.iter()).map(|((((((yni, yi), ki), k3i), k4i), k5i), k6i)| {
            let e = (yni - (yi + dt * (5179.0/57600.0 * ki + 7571.0/16695.0 * k3i + 393.0/640.0 * k4i + -92097.0/339200.0 * k5i + 187.0/2100.0 * k6i + 41.0/840.0 * f(t + dt, &y6)[0])));
            e * e
        }).sum::<f64>().sqrt();

        // Step size control
        let scale = y_new.iter().map(|yi| atol + rtol * yi.abs()).sum::<f64>() / n as f64;
        let dt_new = if err > 0.0 {
            dt * (0.9 * scale / err).powf(0.2).min(5.0).max(0.2)
        } else {
            dt * 5.0
        };

        if err <= scale || dt <= 1e-15 {
            // Accept step
            t += dt;
            y = y_new;
            times.push(t);
            trajectory.push(y.clone());
        }

        dt = dt_new.min(tf - t);
        if dt < 1e-15 {
            break;
        }
    }

    let n_steps_val = times.len() - 1;
    OdeSolution {
        times,
        trajectory,
        n_steps: n_steps_val,
        n_evals,
    }
}
