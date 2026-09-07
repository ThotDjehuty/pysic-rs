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

    // Build Hamiltonian matrix (tridiagonal)
    let mut h_diag = vec![0.0; n];
    let mut h_off = vec![0.0; n - 1];

    for i in 0..n {
        h_diag[i] = 2.0 * coeff + v[i];
    }
    for i in 0..n - 1 {
        h_off[i] = -coeff;
    }

    // Jacobi eigenvalue algorithm (for symmetric tridiagonal)
    let mut eigenvalues = h_diag.clone();
    let mut eigenvectors = vec![vec![0.0; n]; n];
    for i in 0..n {
        eigenvectors[i][i] = 1.0;
    }

    for _iter in 0..200 {
        // Find largest off-diagonal
        let mut max_val = 0.0;
        let mut max_idx = 0;
        for i in 0..n - 1 {
            if h_off[i].abs() > max_val {
                max_val = h_off[i].abs();
                max_idx = i;
            }
        }
        if max_val < 1e-12 {
            break;
        }

        // Givens rotation to zero out h_off[max_idx]
        let p = max_idx;
        let q = p + 1;
        let theta = 0.5 * (2.0 * eigenvalues[p] * h_off[p] - (eigenvalues[p] - eigenvalues[q]).hypot(2.0 * h_off[p])).atan2(eigenvalues[p] - eigenvalues[q] + 2.0 * eigenvalues[p] * h_off[p] / (eigenvalues[p] - eigenvalues[q]).max(1e-15));
        let c = theta.cos();
        let s = theta.sin();

        // Apply rotation to eigenvalues and off-diagonals
        let dp = eigenvalues[p];
        let dq = eigenvalues[q];
        eigenvalues[p] = c * c * dp + s * s * dq - 2.0 * s * c * h_off[p];
        eigenvalues[q] = s * s * dp + c * c * dq + 2.0 * s * c * h_off[p];

        if p > 0 {
            h_off[p - 1] = c * h_off[p - 1] - s * if p > 0 { h_off[p - 1] } else { 0.0 };
        }
        if q < n - 1 {
            h_off[q] = c * h_off[q] - s * if q < n - 1 { h_off[q] } else { 0.0 };
        }
        h_off[p] = 0.0;

        // Apply rotation to eigenvectors
        for i in 0..n {
            let vp = eigenvectors[i][p];
            let vq = eigenvectors[i][q];
            eigenvectors[i][p] = c * vp - s * vq;
            eigenvectors[i][q] = s * vp + c * vq;
        }
    }

    // Sort by eigenvalue
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&a, &b| eigenvalues[a].partial_cmp(&eigenvalues[b]).unwrap());

    let sorted_eigenvalues: Vec<f64> = indices.iter().map(|&i| eigenvalues[i]).take(n_states).collect();
    let sorted_eigenvectors: Vec<Vec<f64>> = indices.iter().map(|&i| {
        let mut v = eigenvectors[i].clone();
        let norm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm > 1e-15 {
            v.iter_mut().for_each(|x| *x /= norm);
        }
        v
    }).take(n_states).collect();

    (sorted_eigenvalues, sorted_eigenvectors)
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
