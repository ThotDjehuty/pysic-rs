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
///
/// # Notes on the error estimator
///
/// Dormand--Prince 5(4) is FSAL: the seventh stage is `f(t+dt, y_new)`, and
/// the embedded fourth-order weights are
/// `[5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40]`.
/// The previous implementation used Fehlberg's `41/840` for the last weight,
/// evaluated `f` *inside* the per-component fold, and took component `[0]` of
/// that evaluation for every component; it then compared an unnormalised
/// L2 norm against a per-component mean of the tolerance, so the accept test
/// was not dimensionally meaningful and the integrator stalled. Those are the
/// `rk45_solve` defects reported in the paper.
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

    let span = tf - t0;
    let mut dt = span / 100.0;
    let dt_min = span.abs() * 1e-13;

    let mut times = vec![t];
    let mut trajectory = vec![y.clone()];
    let mut n_evals = 0usize;

    // Stage coefficients (a_ij)
    const A21: f64 = 1.0 / 5.0;
    const A31: f64 = 3.0 / 40.0;
    const A32: f64 = 9.0 / 40.0;
    const A41: f64 = 44.0 / 45.0;
    const A42: f64 = -56.0 / 15.0;
    const A43: f64 = 32.0 / 9.0;
    const A51: f64 = 19372.0 / 6561.0;
    const A52: f64 = -25360.0 / 2187.0;
    const A53: f64 = 64448.0 / 6561.0;
    const A54: f64 = -212.0 / 729.0;
    const A61: f64 = 9017.0 / 3168.0;
    const A62: f64 = -355.0 / 33.0;
    const A63: f64 = 46732.0 / 5247.0;
    const A64: f64 = 49.0 / 176.0;
    const A65: f64 = -5103.0 / 18656.0;

    // Fifth-order weights (also the seventh stage, FSAL)
    const B1: f64 = 35.0 / 384.0;
    const B3: f64 = 500.0 / 1113.0;
    const B4: f64 = 125.0 / 192.0;
    const B5: f64 = -2187.0 / 6784.0;
    const B6: f64 = 11.0 / 84.0;

    // Embedded fourth-order weights
    const E1: f64 = 5179.0 / 57600.0;
    const E3: f64 = 7571.0 / 16695.0;
    const E4: f64 = 393.0 / 640.0;
    const E5: f64 = -92097.0 / 339200.0;
    const E6: f64 = 187.0 / 2100.0;
    const E7: f64 = 1.0 / 40.0;

    // Node offsets
    const C2: f64 = 1.0 / 5.0;
    const C3: f64 = 3.0 / 10.0;
    const C4: f64 = 4.0 / 5.0;
    const C5: f64 = 8.0 / 9.0;

    let forward = tf >= t0;
    let mut steps = 0usize;

    while steps < max_steps && ((forward && t < tf) || (!forward && t > tf)) {
        if (forward && t + dt > tf) || (!forward && t + dt < tf) {
            dt = tf - t;
        }

        let k1 = f(t, &y);
        let y2: Vec<f64> = (0..n).map(|i| y[i] + dt * A21 * k1[i]).collect();
        let k2 = f(t + C2 * dt, &y2);
        let y3: Vec<f64> = (0..n)
            .map(|i| y[i] + dt * (A31 * k1[i] + A32 * k2[i]))
            .collect();
        let k3 = f(t + C3 * dt, &y3);
        let y4: Vec<f64> = (0..n)
            .map(|i| y[i] + dt * (A41 * k1[i] + A42 * k2[i] + A43 * k3[i]))
            .collect();
        let k4 = f(t + C4 * dt, &y4);
        let y5: Vec<f64> = (0..n)
            .map(|i| y[i] + dt * (A51 * k1[i] + A52 * k2[i] + A53 * k3[i] + A54 * k4[i]))
            .collect();
        let k5 = f(t + C5 * dt, &y5);
        let y6: Vec<f64> = (0..n)
            .map(|i| {
                y[i] + dt
                    * (A61 * k1[i] + A62 * k2[i] + A63 * k3[i] + A64 * k4[i] + A65 * k5[i])
            })
            .collect();
        let k6 = f(t + dt, &y6);

        // Fifth-order update.
        let y_new: Vec<f64> = (0..n)
            .map(|i| {
                y[i] + dt * (B1 * k1[i] + B3 * k3[i] + B4 * k4[i] + B5 * k5[i] + B6 * k6[i])
            })
            .collect();

        // Seventh stage, once, outside any per-component fold.
        let k7 = f(t + dt, &y_new);
        n_evals += 7;

        // Scaled RMS error over components.
        let mut acc = 0.0f64;
        for i in 0..n {
            let y4th = y[i]
                + dt * (E1 * k1[i]
                    + E3 * k3[i]
                    + E4 * k4[i]
                    + E5 * k5[i]
                    + E6 * k6[i]
                    + E7 * k7[i]);
            let sc = atol + rtol * y[i].abs().max(y_new[i].abs());
            let e = (y_new[i] - y4th) / sc;
            acc += e * e;
        }
        let err = (acc / n as f64).sqrt();

        if err <= 1.0 || dt.abs() <= dt_min {
            t += dt;
            y = y_new;
            times.push(t);
            trajectory.push(y.clone());
        }
        steps += 1;

        // PI-free step control with the standard 1/5 exponent and limiters.
        let factor = if err > 0.0 {
            (0.9 * err.powf(-0.2)).clamp(0.2, 5.0)
        } else {
            5.0
        };
        dt *= factor;
        if dt.abs() < dt_min {
            dt = if forward { dt_min } else { -dt_min };
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

#[cfg(test)]
mod tests {
    use super::*;

    /// dy/dt = -y, y(0) = 1 => y(t) = e^{-t}. Graded against the closed form.
    #[test]
    fn exponential_decay_matches_closed_form() {
        let f = |_t: f64, y: &[f64]| vec![-y[0]];
        let sol = rk45_solve(&f, &[1.0], (0.0, 1.0), 1e-10, 1e-12, 100_000);
        let y_end = *sol.trajectory.last().unwrap().first().unwrap();
        let exact = (-1.0f64).exp();
        assert!(
            (y_end - exact).abs() < 1e-9,
            "y(1) = {y_end}, exact {exact}"
        );
        assert!(sol.n_steps < 500, "took {} steps for a trivial ODE", sol.n_steps);
    }

    /// The harmonic oscillator conserves E = (y^2 + v^2)/2.
    #[test]
    fn harmonic_oscillator_conserves_energy() {
        let f = |_t: f64, y: &[f64]| vec![y[1], -y[0]];
        let sol = rk45_solve(&f, &[1.0, 0.0], (0.0, 20.0), 1e-11, 1e-13, 200_000);
        let last = sol.trajectory.last().unwrap();
        let e0 = 0.5;
        let e1 = 0.5 * (last[0] * last[0] + last[1] * last[1]);
        assert!((e1 - e0).abs() < 1e-7, "energy drifted from {e0} to {e1}");
    }

    /// Tightening the tolerance must reduce the error — the property a stalled
    /// integrator fails, because its answer does not depend on rtol at all.
    #[test]
    fn error_decreases_with_tolerance() {
        let f = |_t: f64, y: &[f64]| vec![-y[0]];
        let exact = (-1.0f64).exp();
        let loose = rk45_solve(&f, &[1.0], (0.0, 1.0), 1e-4, 1e-6, 100_000);
        let tight = rk45_solve(&f, &[1.0], (0.0, 1.0), 1e-11, 1e-13, 100_000);
        let e_loose = (loose.trajectory.last().unwrap()[0] - exact).abs();
        let e_tight = (tight.trajectory.last().unwrap()[0] - exact).abs();
        assert!(
            e_tight < e_loose,
            "tightening rtol did not help: {e_loose} -> {e_tight}"
        );
    }
}
