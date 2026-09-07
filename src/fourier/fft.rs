//! FFT operations using the rustfft crate.

use num_complex::Complex64;

/// Compute the FFT of a complex signal.
pub fn fft_complex(signal: &[Complex64]) -> Vec<Complex64> {
    let n = signal.len();
    let mut buffer = signal.to_vec();
    let fft = rustfft::FftPlanner::<f64>::new().plan_fft_forward(n);
    fft.process(&mut buffer);
    buffer
}

/// Compute the inverse FFT of a complex spectrum.
pub fn ifft_complex(spectrum: &[Complex64]) -> Vec<Complex64> {
    let n = spectrum.len();
    let mut buffer = spectrum.to_vec();
    let ifft = rustfft::FftPlanner::<f64>::new().plan_fft_inverse(n);
    ifft.process(&mut buffer);
    for x in buffer.iter_mut() {
        *x /= n as f64;
    }
    buffer
}

/// Compute the FFT of a real signal (returns complex).
pub fn fft_real(signal: &[f64]) -> Vec<Complex64> {
    let complex_signal: Vec<Complex64> = signal.iter().map(|&x| Complex64::new(x, 0.0)).collect();
    fft_complex(&complex_signal)
}

/// Compute the inverse FFT of a real spectrum, returning the real part.
pub fn ifft_real(spectrum: &[Complex64]) -> Vec<f64> {
    let result = ifft_complex(spectrum);
    result.iter().map(|c| c.re).collect()
}

/// Compute the power spectral density (PSD) of a real signal.
pub fn power_spectral_density(signal: &[f64]) -> Vec<f64> {
    let spectrum = fft_real(signal);
    let n = spectrum.len();
    spectrum.iter().take(n / 2 + 1).map(|c| c.norm_sqr() / n as f64).collect()
}

/// Compute frequency bins for FFT of a real signal with given sampling rate.
pub fn fft_frequencies(n: usize, dt: f64) -> Vec<f64> {
    let df = 1.0 / (n as f64 * dt);
    (0..=n / 2).map(|i| i as f64 * df).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_fft_ifft_roundtrip() {
        let signal: Vec<Complex64> = (0..64).map(|i| Complex64::new(i as f64, 0.0)).collect();
        let spectrum = fft_complex(&signal);
        let reconstructed = ifft_complex(&spectrum);
        for (a, b) in signal.iter().zip(reconstructed.iter()) {
            assert_relative_eq!(a.re, b.re, epsilon = 1e-10);
            assert_relative_eq!(a.im, b.im, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_fft_sine_wave() {
        let n = 256;
        let dt = 0.01;
        let freq = 5.0;
        let signal: Vec<f64> = (0..n).map(|i| (2.0 * std::f64::consts::PI * freq * i as f64 * dt).sin()).collect();
        let psd = power_spectral_density(&signal);
        let freqs = fft_frequencies(n, dt);

        // Find peak
        let peak_idx = psd.iter().enumerate().max_by(|a, b| a.1.partial_cmp(b.1).unwrap()).unwrap().0;
        assert!((freqs[peak_idx] - freq).abs() < 1.0, "Peak at {} Hz, expected ~{} Hz", freqs[peak_idx], freq);
    }
}
