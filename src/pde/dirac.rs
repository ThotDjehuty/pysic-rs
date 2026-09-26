//! 1D Dirac equation solver using split-step Fourier method.
//!
//! iℏ ∂ψ/∂t = [c α p + β m c² + V(x)] ψ
//!
//! where α = σ_x ⊗ σ_z, β = σ_z ⊗ I (Dirac representation).

use num_complex::Complex64;

/// Split-step Fourier solver for the 1D Dirac equation with a scalar potential.
///
/// # Arguments
/// * `psi0` - Initial 2-component spinor (length 2*n: [ψ_L, ψ_R] interleaved)
/// * `v` - Scalar potential V(x)
/// * `dx` - Spatial grid spacing
/// * `dt` - Time step
/// * `n_steps` - Number of time steps
/// * `c` - Speed of light (default 1.0 in natural units)
/// * `mass` - Particle mass (default 1.0 in natural units)
/// * `hbar` - Reduced Planck constant (default 1.0)
pub fn dirac_split_step_1d(
    psi0: &[Complex64],
    v: &[f64],
    dx: f64,
    dt: f64,
    n_steps: usize,
    c_speed: f64,
    mass: f64,
    hbar: f64,
) -> Vec<Complex64> {
    let n_grid = v.len();
    let n = psi0.len(); // 2 * n_grid for 2-component spinor
    assert_eq!(n, 2 * n_grid, "psi0 must have length 2 * v.len()");

    let mut psi = psi0.to_vec();

    // Extract left and right components
    let mut psi_l: Vec<Complex64> = (0..n_grid).map(|i| psi[2 * i]).collect();
    let mut psi_r: Vec<Complex64> = (0..n_grid).map(|i| psi[2 * i + 1]).collect();

    // Precompute kinetic phase in k-space
    let dk = 2.0 * std::f64::consts::PI / (n_grid as f64 * dx);
    let mut kinetic_phase = Vec::with_capacity(n_grid);
    for i in 0..n_grid {
        let k = if i <= n_grid / 2 {
            i as f64 * dk
        } else {
            (i as f64 - n_grid as f64) * dk
        };
        let phase = Complex64::new(0.0, -c_speed * k * dt / hbar).exp();
        kinetic_phase.push(phase);
    }

    // Precompute mass phase
    let mass_phase = Complex64::new(0.0, -mass * c_speed * c_speed * dt / hbar).exp();

    let mut fft = rustfft::FftPlanner::<f64>::new().plan_fft_forward(n_grid);
    let mut ifft = rustfft::FftPlanner::<f64>::new().plan_fft_inverse(n_grid);

    for _step in 0..n_steps {
        // Half potential step
        for i in 0..n_grid {
            let pot_phase = Complex64::new(0.0, -v[i] * dt / (2.0 * hbar)).exp();
            psi_l[i] *= pot_phase;
            psi_r[i] *= pot_phase;
        }

        // Half mass rotation: exp(-i σ_x m c² (dt/2) / ħ) = cos(θ/2) I - i sin(θ/2) σ_x.
        // This is unitary by construction, so it conserves ∫|ψ|² exactly.
        let theta = mass * c_speed * c_speed * dt / hbar;
        let cos_half = (theta / 2.0).cos();
        let i_sin_half = Complex64::new(0.0, (theta / 2.0).sin());
        for i in 0..n_grid {
            let new_l = cos_half * psi_l[i] - i_sin_half * psi_r[i];
            let new_r = -i_sin_half * psi_l[i] + cos_half * psi_r[i];
            psi_l[i] = new_l;
            psi_r[i] = new_r;
        }

        // FFT to k-space
        let mut buf_l: Vec<Complex64> = psi_l.clone();
        let mut buf_r: Vec<Complex64> = psi_r.clone();
        fft.process(&mut buf_l);
        fft.process(&mut buf_r);

        // Kinetic step in k-space (helicity basis)
        for i in 0..n_grid {
            let p = hbar * (if i <= n_grid / 2 { i as f64 * dk } else { (i as f64 - n_grid as f64) * dk });
            let phase_pos = kinetic_phase[i];
            let phase_neg = Complex64::new(0.0, c_speed * p * dt / hbar).exp();

            let new_l = phase_pos * buf_l[i];
            let new_r = phase_neg * buf_r[i];
            buf_l[i] = new_l;
            buf_r[i] = new_r;
        }

        // IFFT back to x-space
        ifft.process(&mut buf_l);
        ifft.process(&mut buf_r);
        for i in 0..n_grid {
            psi_l[i] = buf_l[i] / n_grid as f64;
            psi_r[i] = buf_r[i] / n_grid as f64;
        }

        // Second half mass rotation, completing the Strang splitting
        // V/2 · M/2 · K · M/2 · V/2, which is second-order accurate in dt.
        for i in 0..n_grid {
            let new_l = cos_half * psi_l[i] - i_sin_half * psi_r[i];
            let new_r = -i_sin_half * psi_l[i] + cos_half * psi_r[i];
            psi_l[i] = new_l;
            psi_r[i] = new_r;
        }

        // Half potential step
        for i in 0..n_grid {
            let pot_phase = Complex64::new(0.0, -v[i] * dt / (2.0 * hbar)).exp();
            psi_l[i] *= pot_phase;
            psi_r[i] *= pot_phase;
        }
    }

    // Reinterleave
    let mut result = vec![Complex64::new(0.0, 0.0); n];
    for i in 0..n_grid {
        result[2 * i] = psi_l[i];
        result[2 * i + 1] = psi_r[i];
    }
    result
}
