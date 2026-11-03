//! Time-dependent and time-independent Schrödinger equation solvers.

use num_complex::Complex64;

/// Split-step Fourier method for the 1D TDSE: iℏ ∂ψ/∂t = [-ℏ²/2m ∂²/∂x² + V(x)] ψ.
///
/// Uses operator splitting: kinetic half-step in k-space, potential step in x-space.
///
/// # Arguments
/// * `psi0` - Initial wavefunction (complex)
/// * `v` - Potential V(x)
/// * `dx` - Spatial grid spacing
/// * `dt` - Time step
/// * `n_steps` - Number of time steps
/// * `mass` - Particle mass (default 1.0 in natural units)
/// * `hbar` - Reduced Planck constant (default 1.0 in natural units)
pub fn schrodinger_split_step_1d(
    psi0: &[Complex64],
    v: &[f64],
    dx: f64,
    dt: f64,
    n_steps: usize,
    mass: f64,
    hbar: f64,
) -> Vec<Complex64> {
    let n = psi0.len();
    let mut psi = psi0.to_vec();

    // Precompute kinetic phase in k-space
    let dk = 2.0 * std::f64::consts::PI / (n as f64 * dx);
    let mut kinetic_phase = Vec::with_capacity(n);
    for i in 0..n {
        let k = if i <= n / 2 {
            i as f64 * dk
        } else {
            (i as f64 - n as f64) * dk
        };
        let ek = hbar * hbar * k * k / (2.0 * mass);
        kinetic_phase.push(Complex64::new(0.0, -ek * dt / hbar).exp());
    }

    // Precompute potential phase
    let mut potential_phase = Vec::with_capacity(n);
    for i in 0..n {
        potential_phase.push(Complex64::new(0.0, -v[i] * dt / (2.0 * hbar)).exp());
    }

    // FFT plan (in-place via rustfft)
    let mut fft = rustfft::FftPlanner::<f64>::new().plan_fft_forward(n);
    let mut ifft = rustfft::FftPlanner::<f64>::new().plan_fft_inverse(n);

    for _ in 0..n_steps {
        // Half potential step
        for i in 0..n {
            psi[i] *= potential_phase[i];
        }

        // FFT to k-space
        let mut buf: Vec<Complex64> = psi.clone();
        fft.process(&mut buf);

        // Kinetic step in k-space
        for i in 0..n {
            buf[i] *= kinetic_phase[i];
        }

        // IFFT back to x-space
        ifft.process(&mut buf);
        for i in 0..n {
            buf[i] /= n as f64;
        }
        psi = buf;

        // Half potential step
        for i in 0..n {
            psi[i] *= potential_phase[i];
        }
    }

    psi
}

/// Solve the 1D time-independent Schrödinger equation using matrix diagonalization.
///
/// Discretizes -ℏ²/2m d²ψ/dx² + V(x)ψ = Eψ using finite differences.
///
/// # Returns
/// * `energies` - Sorted eigenvalues (lowest first)
/// * `eigenstates` - Corresponding eigenvectors (normalized)
pub fn schrodinger_eigen_1d(
    v: &[f64],
    dx: f64,
    n_states: usize,
    mass: f64,
    hbar: f64,
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = v.len();
    let coeff = hbar * hbar / (2.0 * mass * dx * dx);

    // Symmetric tridiagonal Hamiltonian: diag 2*coeff + V_i, offdiag -coeff.
    let mut d: Vec<f64> = (0..n).map(|i| 2.0 * coeff + v[i]).collect();
    let mut e: Vec<f64> = vec![0.0; n];
    for i in 0..n.saturating_sub(1) {
        e[i] = -coeff;
    }

    // Accumulate eigenvectors in z (columns are eigenvectors).
    let mut z = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        z[i][i] = 1.0;
    }

    tridiagonal_ql_implicit(&mut d, &mut e, &mut z);

    // Sort ascending by eigenvalue and take the lowest n_states.
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| d[a].partial_cmp(&d[b]).unwrap());

    let take = n_states.min(n);
    let eigenvalues: Vec<f64> = idx.iter().take(take).map(|&i| d[i]).collect();
    let eigenvectors: Vec<Vec<f64>> = idx
        .iter()
        .take(take)
        .map(|&col| {
            let mut vec_col: Vec<f64> = (0..n).map(|row| z[row][col]).collect();
            let norm: f64 = vec_col.iter().map(|x| x * x).sum::<f64>().sqrt();
            if norm > 1e-300 {
                vec_col.iter_mut().for_each(|x| *x /= norm);
            }
            vec_col
        })
        .collect();

    (eigenvalues, eigenvectors)
}

/// QL with implicit shifts for a symmetric tridiagonal matrix (EISPACK
/// `imtql2` / Numerical Recipes `tqli`). `d` holds the diagonal and is
/// overwritten with the eigenvalues; `e` holds the subdiagonal in
/// `e[0..n-1]` and is destroyed; the columns of `z` are rotated into the
/// eigenvectors.
///
/// This replaces a "Jacobi sweep" that computed a malformed rotation angle and
/// updated the off-diagonals from themselves
/// (`h_off[p-1] = c*h_off[p-1] - s*h_off[p-1]`), so it never deflated and
/// returned values on the order of the matrix diagonal — i.e. `O(1/dx^2)`
/// rather than the physical spectrum. A plane-rotation sweep is in any case
/// the wrong shape for a tridiagonal, because each rotation creates fill-in
/// that has to be chased back onto the band.
fn tridiagonal_ql_implicit(d: &mut [f64], e: &mut [f64], z: &mut [Vec<f64>]) {
    let n = d.len();
    if n == 0 {
        return;
    }
    if n > 0 {
        e[n - 1] = 0.0;
    }

    let sign = |a: f64, b: f64| if b >= 0.0 { a.abs() } else { -a.abs() };

    for l in 0..n {
        let mut iter = 0;
        loop {
            // Look for a negligible subdiagonal element to split on.
            let mut m = l;
            while m + 1 < n {
                let dd = d[m].abs() + d[m + 1].abs();
                if e[m].abs() <= f64::EPSILON * dd {
                    break;
                }
                m += 1;
            }
            if m == l {
                break;
            }
            iter += 1;
            if iter > 50 {
                break; // refuse to spin forever
            }

            let mut g = (d[l + 1] - d[l]) / (2.0 * e[l]);
            let mut r = g.hypot(1.0);
            g = d[m] - d[l] + e[l] / (g + sign(r, g));
            let mut s = 1.0f64;
            let mut c = 1.0f64;
            let mut p = 0.0f64;

            let mut i = m;
            while i > l {
                i -= 1;
                let mut f = s * e[i];
                let b = c * e[i];
                r = f.hypot(g);
                e[i + 1] = r;
                if r == 0.0 {
                    d[i + 1] -= p;
                    e[m] = 0.0;
                    break;
                }
                s = f / r;
                c = g / r;
                g = d[i + 1] - p;
                r = (d[i] - g) * s + 2.0 * c * b;
                p = s * r;
                d[i + 1] = g + p;
                g = c * r - b;

                for k in 0..n {
                    f = z[k][i + 1];
                    z[k][i + 1] = s * z[k][i] + c * f;
                    z[k][i] = c * z[k][i] - s * f;
                }
            }
            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    #[test]
    fn test_split_step_normalization() {
        let n = 256;
        let dx = 0.1;
        let x: Vec<f64> = (0..n).map(|i| (i as f64 - n as f64 / 2.0) * dx).collect();

        // Gaussian wavepacket
        let sigma = 2.0;
        let k0 = 1.0;
        let mut psi0: Vec<Complex64> = x.iter().map(|&xi| {
            let norm = 1.0 / (sigma * std::f64::consts::PI.sqrt()).sqrt();
            Complex64::new(norm * (-xi * xi / (2.0 * sigma * sigma)).exp(), 0.0)
                * Complex64::new(0.0, k0 * xi).exp()
        }).collect();

        // Free potential
        let v = vec![0.0; n];
        let psi = schrodinger_split_step_1d(&psi0, &v, dx, 0.01, 100, 1.0, 1.0);

        // Check normalization
        let norm2: f64 = psi.iter().map(|p| p.norm_sqr() * dx).sum();
        assert!((norm2 - 1.0).abs() < 0.01, "norm = {}", norm2);
    }
}
