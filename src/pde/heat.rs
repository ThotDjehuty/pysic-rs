//! Heat/diffusion equation: ∂u/∂t = α ∇²u.
//!
//! Crank-Nicolson (1D) and ADI (2D) methods.

/// 1D heat equation via Crank-Nicolson (unconditionally stable, 2nd order).
///
/// ∂u/∂t = α ∂²u/∂x²
pub fn heat_crank_nicolson_1d(
    u0: &[f64],
    dx: f64,
    dt: f64,
    alpha: f64,
    n_steps: usize,
) -> Vec<f64> {
    let n = u0.len();
    let r = alpha * dt / (dx * dx);
    let mut u = u0.to_vec();

    // Set up tridiagonal system for Crank-Nicolson
    // (I - r/2 D²) u^{n+1} = (I + r/2 D²) u^n
    let n_inner = n - 2;

    for _step in 0..n_steps {
        // RHS: (I + r/2 D²) u^n
        let mut rhs = vec![0.0; n_inner];
        for i in 0..n_inner {
            let idx = i + 1;
            rhs[i] = (1.0 + r) * u[idx] - 0.5 * r * (u[idx - 1] + u[idx + 1]);
        }
        // Boundary conditions (Dirichlet: u[0] and u[n-1] fixed)
        rhs[0] += 0.5 * r * u[0];
        rhs[n_inner - 1] += 0.5 * r * u[n - 1];

        // Solve tridiagonal: (1+r) on diagonal, -r/2 on off-diagonals
        let mut diag = vec![1.0 + r; n_inner];
        let mut lower = vec![-0.5 * r; n_inner - 1];
        let mut upper = vec![-0.5 * r; n_inner - 1];

        // Thomas algorithm
        for i in 1..n_inner {
            let m = lower[i - 1] / diag[i - 1];
            diag[i] -= m * upper[i - 1];
            rhs[i] -= m * rhs[i - 1];
        }

        let mut x = vec![0.0; n_inner];
        x[n_inner - 1] = rhs[n_inner - 1] / diag[n_inner - 1];
        for i in (0..n_inner - 1).rev() {
            x[i] = (rhs[i] - upper[i] * x[i + 1]) / diag[i];
        }

        // Update
        u[0] = u[0]; // Dirichlet BC
        for i in 0..n_inner {
            u[i + 1] = x[i];
        }
        u[n - 1] = u[n - 1]; // Dirichlet BC
    }

    u
}

/// 2D heat equation via ADI (Alternating Direction Implicit).
///
/// ∂u/∂t = α (∂²u/∂x² + ∂²u/∂y²)
pub fn heat_crank_nicolson_2d(
    u0: &[f64],
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    dt: f64,
    alpha: f64,
    n_steps: usize,
) -> Vec<f64> {
    let mut u = u0.to_vec();
    let rx = alpha * dt / (2.0 * dx * dx);
    let ry = alpha * dt / (2.0 * dy * dy);

    for _step in 0..n_steps {
        let mut u_temp = u.clone();

        // X-sweep: implicit in x, explicit in y
        for j in 1..ny - 1 {
            let mut a = vec![-rx; nx - 2];
            let mut b = vec![1.0 + 2.0 * rx; nx - 2];
            let mut c = vec![-rx; nx - 2];
            let mut d = vec![0.0; nx - 2];

            for i in 0..nx - 2 {
                let idx = (j) * nx + (i + 1);
                d[i] = u[idx] + ry * (u[(j - 1) * nx + (i + 1)] - 2.0 * u[idx] + u[(j + 1) * nx + (i + 1)]);
            }
            d[0] += rx * u[j * nx];
            d[nx - 3] += rx * u[j * nx + nx - 1];

            // Thomas algorithm
            for i in 1..nx - 2 {
                let m = a[i - 1] / b[i - 1];
                b[i] -= m * c[i - 1];
                d[i] -= m * d[i - 1];
            }
            let mut x = vec![0.0; nx - 2];
            x[nx - 3] = d[nx - 3] / b[nx - 3];
            for i in (0..nx - 3).rev() {
                x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
            }
            for i in 0..nx - 2 {
                u_temp[j * nx + (i + 1)] = x[i];
            }
        }

        u = u_temp.clone();

        // Y-sweep: implicit in y, explicit in x
        for i in 1..nx - 1 {
            let mut a = vec![-ry; ny - 2];
            let mut b = vec![1.0 + 2.0 * ry; ny - 2];
            let mut c = vec![-ry; ny - 2];
            let mut d = vec![0.0; ny - 2];

            for j in 0..ny - 2 {
                let idx = (j + 1) * nx + i;
                d[j] = u[idx] + rx * (u[(j + 1) * nx + i - 1] - 2.0 * u[idx] + u[(j + 1) * nx + i + 1]);
            }
            d[0] += ry * u[i];
            d[ny - 3] += ry * u[(ny - 1) * nx + i];

            for j in 1..ny - 2 {
                let m = a[j - 1] / b[j - 1];
                b[j] -= m * c[j - 1];
                d[j] -= m * d[j - 1];
            }
            let mut y = vec![0.0; ny - 2];
            y[ny - 3] = d[ny - 3] / b[ny - 3];
            for j in (0..ny - 3).rev() {
                y[j] = (d[j] - c[j] * y[j + 1]) / b[j];
            }
            for j in 0..ny - 2 {
                u[(j + 1) * nx + i] = y[j];
            }
        }
    }

    u
}
