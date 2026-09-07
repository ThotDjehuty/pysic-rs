//! Poisson/Laplace equation solver via FFT.
//!
//! -∇²u = f on a 2D rectangular domain with zero Dirichlet BC.

/// Solve -∇²u = f on [0,Lx] × [0,Ly] with u=0 on boundary.
///
/// Uses discrete sine transform (DST-I) via FFT.
pub fn poisson_fft_2d(
    f: &[f64],
    nx: usize,
    ny: usize,
) -> Vec<f64> {
    let mut u = vec![0.0; nx * ny];

    // Compute DST-I of f
    let mut f_hat = vec![0.0; nx * ny];
    for k in 1..ny {
        for j in 1..nx {
            let mut sum = 0.0;
            for m in 1..ny {
                for i in 1..nx {
                    sum += f[m * nx + i]
                        * (std::f64::consts::PI * j as f64 * i as f64 / nx as f64).sin()
                        * (std::f64::consts::PI * k as f64 * m as f64 / ny as f64).sin();
                }
            }
            f_hat[k * nx + j] = 4.0 / (nx as f64 * ny as f64) * sum;
        }
    }

    // Solve in spectral space: u_hat = f_hat / (λ_j + λ_k)
    for k in 1..ny {
        for j in 1..nx {
            let lambda_j = 4.0 / (1.0 * 1.0) * (std::f64::consts::PI * j as f64 / (2.0 * nx as f64)).sin().powi(2);
            let lambda_k = 4.0 / (1.0 * 1.0) * (std::f64::consts::PI * k as f64 / (2.0 * ny as f64)).sin().powi(2);
            let denom = lambda_j + lambda_k;
            if denom > 1e-15 {
                u[k * nx + j] = f_hat[k * nx + j] / denom;
            }
        }
    }

    // Inverse DST-I
    let mut result = vec![0.0; nx * ny];
    for m in 1..ny {
        for i in 1..nx {
            let mut sum = 0.0;
            for k in 1..ny {
                for j in 1..nx {
                    sum += u[k * nx + j]
                        * (std::f64::consts::PI * j as f64 * i as f64 / nx as f64).sin()
                        * (std::f64::consts::PI * k as f64 * m as f64 / ny as f64).sin();
                }
            }
            result[m * nx + i] = sum;
        }
    }

    result
}

/// Solve the 2D Poisson equation -∇²u = f using SOR (Successive Over-Relaxation).
///
/// Useful for irregular domains or when FFT-based methods don't apply.
pub fn poisson_sor_2d(
    f: &[f64],
    nx: usize,
    ny: usize,
    omega: f64,
    max_iter: usize,
    tol: f64,
) -> Vec<f64> {
    let mut u = vec![0.0; nx * ny];
    let dx2 = 1.0;
    let dy2 = 1.0;

    for _iter in 0..max_iter {
        let mut max_change = 0.0f64;

        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = j * nx + i;
                let neighbors = u[idx - 1] + u[idx + 1] + u[(j - 1) * nx + i] + u[(j + 1) * nx + i];
                let u_new = 0.25 * (neighbors + f[idx] * dx2);
                let change = omega * (u_new - u[idx]);
                u[idx] += change;
                max_change = max_change.max(change.abs());
            }
        }

        if max_change < tol {
            break;
        }
    }

    u
}
