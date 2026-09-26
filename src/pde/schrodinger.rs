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

/// Implicit QL eigen-decomposition of a real symmetric tridiagonal matrix
/// (EISPACK `tql2` / Numerical Recipes `tqli`, with Wilkinson shifts).
///
/// `d` holds the diagonal and is overwritten with the eigenvalues; `e` holds
/// the sub-diagonal in `e[0..n-1]` and is destroyed. `z` starts as the identity
/// and accumulates the rotations, so on return **column `j` of `z`** is the
/// unit eigenvector belonging to `d[j]`.
///
/// Returns `false` if any eigenvalue failed to converge within 30 sweeps.
fn tql2(d: &mut [f64], e: &mut [f64], z: &mut [Vec<f64>]) -> bool {
    let n = d.len();
    if n == 0 {
        return true;
    }
    if n == 1 {
        return true;
    }

    // Shift the sub-diagonal so that e[i] couples d[i] and d[i-1].
    for i in 1..n {
        e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;

    for l in 0..n {
        let mut iter = 0;
        loop {
            // Look for a single small sub-diagonal element to split the matrix.
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
            if iter > 30 {
                return false;
            }

            // Wilkinson shift.
            let mut g = (d[l + 1] - d[l]) / (2.0 * e[l]);
            let mut r = g.hypot(1.0);
            g = d[m] - d[l] + e[l] / (g + if g >= 0.0 { r.abs() } else { -r.abs() });

            let mut s = 1.0;
            let mut c = 1.0;
            let mut p = 0.0;
            let mut underflow = false;

            for i in (l..m).rev() {
                let mut f = s * e[i];
                let b = c * e[i];
                r = f.hypot(g);
                e[i + 1] = r;
                if r == 0.0 {
                    // Recover from underflow.
                    d[i + 1] -= p;
                    e[m] = 0.0;
                    underflow = true;
                    break;
                }
                s = f / r;
                c = g / r;
                g = d[i + 1] - p;
                r = (d[i] - g) * s + 2.0 * c * b;
                p = s * r;
                d[i + 1] = g + p;
                g = c * r - b;

                // Accumulate the rotation into the eigenvector matrix.
                for k in 0..n {
                    f = z[k][i + 1];
                    z[k][i + 1] = s * z[k][i] + c * f;
                    z[k][i] = c * z[k][i] - s * f;
                }
            }

            if underflow {
                continue;
            }
            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
        }
    }
    true
}

/// Solve the 1D time-independent Schrödinger equation by direct
/// diagonalization of the finite-difference Hamiltonian.
///
/// Discretizes `-ℏ²/2m ψ'' + V(x)ψ = Eψ` on a uniform grid with **Dirichlet
/// boundary conditions** (ψ vanishes one node outside each end), giving the
/// symmetric tridiagonal matrix
///
/// ```text
/// H_ii = ℏ²/(m dx²) + V_i ,      H_i,i±1 = -ℏ²/(2 m dx²).
/// ```
///
/// The spectrum is obtained with the implicit QL algorithm, so for a flat
/// potential the result reproduces the exact discrete eigenvalues
/// `E_k = (2ℏ²/m dx²) sin²(kπ / 2(n+1))` to machine precision, and converges
/// to the continuum particle-in-a-box energies as `dx → 0`.
///
/// # Returns
/// * `energies` — the `n_states` lowest eigenvalues, ascending.
/// * `eigenstates` — the matching eigenvectors, each with unit Euclidean norm
///   (`Σ ψ_i² = 1`). Divide by `sqrt(dx)` for the continuum normalization
///   `∫|ψ|² dx = 1`.
pub fn schrodinger_eigen_1d(
    v: &[f64],
    dx: f64,
    n_states: usize,
    mass: f64,
    hbar: f64,
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = v.len();
    if n == 0 || n_states == 0 {
        return (Vec::new(), Vec::new());
    }
    let coeff = hbar * hbar / (2.0 * mass * dx * dx);

    // Symmetric tridiagonal Hamiltonian.
    let mut d: Vec<f64> = (0..n).map(|i| 2.0 * coeff + v[i]).collect();
    let mut e: Vec<f64> = vec![0.0; n];
    for i in 1..n {
        e[i] = -coeff;
    }

    // z starts as the identity and ends holding the eigenvectors as columns.
    let mut z: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut row = vec![0.0; n];
            row[i] = 1.0;
            row
        })
        .collect();

    tql2(&mut d, &mut e, &mut z);

    // Sort ascending and take the lowest n_states.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| d[a].partial_cmp(&d[b]).unwrap_or(std::cmp::Ordering::Equal));
    let keep = n_states.min(n);

    let energies: Vec<f64> = order.iter().take(keep).map(|&j| d[j]).collect();
    let eigenstates: Vec<Vec<f64>> = order
        .iter()
        .take(keep)
        .map(|&j| {
            // Column j of z is the eigenvector for d[j].
            let mut psi: Vec<f64> = (0..n).map(|k| z[k][j]).collect();
            let norm: f64 = psi.iter().map(|x| x * x).sum::<f64>().sqrt();
            if norm > 1e-300 {
                psi.iter_mut().for_each(|x| *x /= norm);
            }
            // Fix the global sign so the result is deterministic: make the
            // first significant component positive.
            if let Some(&first) = psi.iter().find(|x| x.abs() > 1e-12) {
                if first < 0.0 {
                    psi.iter_mut().for_each(|x| *x = -*x);
                }
            }
            psi
        })
        .collect();

    (energies, eigenstates)
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

    /// The finite-difference Hamiltonian with a flat potential has the exact
    /// analytic spectrum E_k = (2ℏ²/m dx²) sin²(kπ / 2(n+1)). Reproducing it to
    /// machine precision pins the eigensolver itself, with no discretization
    /// error in the way.
    #[test]
    fn test_eigen_flat_potential_matches_exact_discrete_spectrum() {
        let n = 80;
        let dx = 0.05;
        let (hbar, mass) = (1.0, 1.0);
        let v = vec![0.0; n];
        let n_states = 6;

        let (energies, states) = schrodinger_eigen_1d(&v, dx, n_states, mass, hbar);
        assert_eq!(energies.len(), n_states);

        let pref = 2.0 * hbar * hbar / (mass * dx * dx);
        for k in 1..=n_states {
            let theta = (k as f64) * std::f64::consts::PI / (2.0 * (n as f64 + 1.0));
            let exact = pref * theta.sin() * theta.sin();
            let got = energies[k - 1];
            assert!(
                (got - exact).abs() < 1e-9 * exact.max(1.0),
                "state {k}: got {got}, exact {exact}"
            );
        }

        // Ascending order, and eigenvectors are unit-norm.
        for w in energies.windows(2) {
            assert!(w[0] <= w[1], "energies not sorted: {:?}", energies);
        }
        for psi in &states {
            let norm: f64 = psi.iter().map(|x| x * x).sum();
            assert!((norm - 1.0).abs() < 1e-10, "eigenvector norm = {norm}");
        }
    }

    /// The ground state of a box has no interior node; the k-th excited state
    /// has exactly k. This catches eigenvectors that are silently mismatched
    /// to their eigenvalues.
    #[test]
    fn test_eigen_eigenvectors_have_correct_node_count() {
        let n = 120;
        let dx = 1.0 / (n as f64 + 1.0);
        let v = vec![0.0; n];
        let (_, states) = schrodinger_eigen_1d(&v, dx, 4, 1.0, 1.0);

        for (k, psi) in states.iter().enumerate() {
            let nodes = psi
                .windows(2)
                .filter(|w| w[0] * w[1] < 0.0)
                .count();
            assert_eq!(nodes, k, "state {k} has {nodes} interior nodes");
        }
    }

    /// The canonical physics check: the harmonic oscillator V = ½mω²x² has
    /// E_n = (n + ½)ℏω exactly. Node counting is done above a relative
    /// amplitude threshold, because in the classically forbidden tails |ψ| is
    /// ~1e-300 and its sign is pure floating-point noise.
    #[test]
    fn test_eigen_harmonic_oscillator_spectrum() {
        // Kept modest: tql2 accumulates eigenvectors in O(n³), and `cargo test`
        // runs unoptimized.
        let n = 400;
        let dx = 20.0 / n as f64;
        let x: Vec<f64> = (0..n).map(|i| (i as f64 - n as f64 / 2.0) * dx).collect();
        let v: Vec<f64> = x.iter().map(|&xi| 0.5 * xi * xi).collect();

        let (energies, states) = schrodinger_eigen_1d(&v, dx, 5, 1.0, 1.0);

        for (k, &e) in energies.iter().enumerate() {
            let exact = k as f64 + 0.5;
            assert!(
                (e - exact).abs() < 5e-3,
                "level {k}: got {e}, exact {exact}"
            );
        }

        for (k, psi) in states.iter().enumerate() {
            let max = psi.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
            let sig: Vec<f64> = psi
                .iter()
                .copied()
                .filter(|x| x.abs() > 1e-8 * max)
                .collect();
            let nodes = sig.windows(2).filter(|w| w[0] * w[1] < 0.0).count();
            assert_eq!(nodes, k, "level {k} has {nodes} nodes, expected {k}");
        }
    }

    /// Continuum check: as dx → 0 the box energies approach n²π²ℏ²/(2mL²).
    #[test]
    fn test_eigen_converges_to_particle_in_a_box() {
        let l = 1.0_f64;
        let mut last_err = f64::INFINITY;
        for &n in &[100_usize, 200, 400] {
            let dx = l / (n as f64 + 1.0);
            let v = vec![0.0; n];
            let (energies, _) = schrodinger_eigen_1d(&v, dx, 3, 1.0, 1.0);
            let exact_1 = std::f64::consts::PI.powi(2) / (2.0 * l * l);
            let err = (energies[0] - exact_1).abs() / exact_1;
            assert!(err < last_err, "error did not decrease: {err} vs {last_err}");
            last_err = err;
        }
        assert!(last_err < 1e-4, "final relative error {last_err} too large");
    }

    /// The Dirac split-step propagator is a product of unitary factors, so the
    /// norm must be conserved. Before the mass-rotation fix this test blew up
    /// to ~1e57.
    #[test]
    fn test_dirac_split_step_conserves_norm() {
        use crate::pde::dirac::dirac_split_step_1d;
        let n_grid = 256;
        let dx = 0.1;
        let x: Vec<f64> = (0..n_grid).map(|i| (i as f64 - n_grid as f64 / 2.0) * dx).collect();

        // Interleaved 2-spinor [ψ_L, ψ_R] with a Gaussian upper component.
        let mut psi0 = Vec::with_capacity(2 * n_grid);
        for &xi in &x {
            let g = (-xi * xi / 2.0).exp();
            psi0.push(Complex64::new(g, 0.0));
            psi0.push(Complex64::new(0.0, 0.0));
        }
        let v = vec![0.0; n_grid];
        let norm0: f64 = psi0.iter().map(|c| c.norm_sqr()).sum::<f64>() * dx;

        for &n_steps in &[10_usize, 100, 500] {
            let psi = dirac_split_step_1d(&psi0, &v, dx, 1e-3, n_steps, 1.0, 1.0, 1.0);
            let norm: f64 = psi.iter().map(|c| c.norm_sqr()).sum::<f64>() * dx;
            let drift = (norm - norm0).abs() / norm0;
            assert!(drift < 1e-10, "after {n_steps} steps norm drift = {drift:e}");
        }
    }

    /// A zero-mass, zero-potential spinor component is a pure chiral mover: the
    /// upper component must translate rigidly at speed c without changing shape.
    #[test]
    fn test_dirac_massless_component_translates_at_c() {
        use crate::pde::dirac::dirac_split_step_1d;
        let n_grid = 512;
        let dx = 0.05;
        let c = 1.0;
        let x: Vec<f64> = (0..n_grid).map(|i| (i as f64 - n_grid as f64 / 2.0) * dx).collect();

        let mut psi0 = Vec::with_capacity(2 * n_grid);
        for &xi in &x {
            let g = (-(xi + 4.0) * (xi + 4.0) / (2.0 * 1.0_f64.powi(2))).exp();
            psi0.push(Complex64::new(g, 0.0)); // ψ_L
            psi0.push(Complex64::new(0.0, 0.0)); // ψ_R
        }
        let v = vec![0.0; n_grid];

        let dt = 1e-3;
        let n_steps = 2000;
        let t = dt * n_steps as f64;
        let psi = dirac_split_step_1d(&psi0, &v, dx, dt, n_steps, c, 0.0, 1.0);

        // With m = 0 the L component obeys ∂_t ψ_L = -c ∂_x ψ_L, so the packet
        // centred at -4 moves to -4 + c t.
        let centre_expected = -4.0 + c * t;
        let mut num = 0.0;
        let mut den = 0.0;
        for (i, &xi) in x.iter().enumerate() {
            let w = psi[2 * i].norm_sqr();
            num += xi * w;
            den += w;
        }
        let centre_got = num / den;
        assert!(
            (centre_got - centre_expected).abs() < 5.0 * dx,
            "massless packet centre {centre_got}, expected {centre_expected}"
        );
    }
}
