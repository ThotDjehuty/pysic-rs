//! Python bindings for Fourier analysis.

use pyo3::prelude::*;
use num_complex::Complex64;

pub fn register_python_functions(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(fft_py, m)?)?;
    m.add_function(wrap_pyfunction!(ifft_py, m)?)?;
    m.add_function(wrap_pyfunction!(spectral_derivative_py, m)?)?;
    m.add_function(wrap_pyfunction!(welch_psd_py, m)?)?;
    Ok(())
}

#[pyfunction]
fn fft_py(real: Vec<f64>, imag: Vec<f64>) -> PyResult<(Vec<f64>, Vec<f64>)> {
    let signal: Vec<Complex64> = real.iter().zip(imag.iter()).map(|(r, i)| Complex64::new(*r, *i)).collect();
    let result = super::fft::fft_complex(&signal);
    Ok((result.iter().map(|c| c.re).collect(), result.iter().map(|c| c.im).collect()))
}

#[pyfunction]
fn ifft_py(real: Vec<f64>, imag: Vec<f64>) -> PyResult<(Vec<f64>, Vec<f64>)> {
    let spectrum: Vec<Complex64> = real.iter().zip(imag.iter()).map(|(r, i)| Complex64::new(*r, *i)).collect();
    let result = super::fft::ifft_complex(&spectrum);
    Ok((result.iter().map(|c| c.re).collect(), result.iter().map(|c| c.im).collect()))
}

#[pyfunction]
fn spectral_derivative_py(f: Vec<f64>, dx: f64, order: usize) -> PyResult<Vec<f64>> {
    Ok(super::spectral::spectral_derivative(&f, dx, order))
}

#[pyfunction]
fn welch_psd_py(signal: Vec<f64>, nperseg: usize, noverlap: usize) -> PyResult<(Vec<f64>, Vec<f64>)> {
    Ok(super::spectral::welch_psd(&signal, nperseg, noverlap))
}
