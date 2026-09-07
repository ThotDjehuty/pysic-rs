//! Spectral methods: derivatives, power spectral density estimation.

/// Compute the nth-order derivative of f using spectral (FFT) method.
///
/// df^n/dx^n = IFFT[(ik)^n · FFT[f]]
pub fn spectral_derivative(f: &[f64], dx: f64, order: usize) -> Vec<f64> {
    let n = f.len();
    let mut spectrum = super::fft::fft_real(f);
    let dk = 2.0 * std::f64::consts::PI / (n as f64 * dx);

    for i in 0..n {
        let k = if i <= n / 2 {
            i as f64 * dk
        } else {
            (i as f64 - n as f64) * dk
        };
        // (ik)^n
        let phase = k * (order as f64 * std::f64::consts::FRAC_PI_2); // i^n = e^(i n π/2)
        let magnitude = k.powi(order as i32);
        let factor = num_complex::Complex64::from_polar(magnitude, phase);
        spectrum[i] *= factor;
    }

    let result = super::fft::ifft_complex(&spectrum);
    result.iter().map(|c| c.re).collect()
}

/// Welch's method for estimating the power spectral density.
///
/// Splits signal into overlapping segments, windows each, computes FFT, averages.
pub fn welch_psd(signal: &[f64], nperseg: usize, noverlap: usize) -> (Vec<f64>, Vec<f64>) {
    let n = signal.len();
    if nperseg > n {
        return (vec![], vec![]);
    }

    let step = nperseg - noverlap;
    let n_segments = (n - nperseg) / step + 1;

    // Hann window
    let window: Vec<f64> = (0..nperseg).map(|i| {
        0.5 * (1.0 - (2.0 * std::f64::consts::PI * i as f64 / nperseg as f64).cos())
    }).collect();
    let window_norm: f64 = window.iter().map(|w| w * w).sum::<f64>();

    let mut psd = vec![0.0; nperseg / 2 + 1];

    for seg in 0..n_segments {
        let start = seg * step;
        let segment: Vec<f64> = signal[start..start + nperseg].iter()
            .zip(window.iter())
            .map(|(s, w)| s * w)
            .collect();

        let spectrum = super::fft::fft_real(&segment);
        for i in 0..=nperseg / 2 {
            psd[i] += spectrum[i].norm_sqr();
        }
    }

    // Average and normalize
    for val in psd.iter_mut() {
        *val /= n_segments as f64 * window_norm;
    }

    let freqs = super::fft::fft_frequencies(nperseg, 1.0);
    (freqs, psd)
}

/// Compute the power spectrum via periodogram.
pub fn periodogram(signal: &[f64], dx: f64) -> (Vec<f64>, Vec<f64>) {
    let n = signal.len();
    let spectrum = super::fft::fft_real(signal);
    let freqs = super::fft::fft_frequencies(n, dx);
    let psd: Vec<f64> = spectrum.iter().take(n / 2 + 1).map(|c| c.norm_sqr() / (n as f64 * dx)).collect();
    (freqs, psd)
}
