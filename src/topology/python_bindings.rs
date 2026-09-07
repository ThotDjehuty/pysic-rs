//! Python bindings for topology.

use pyo3::prelude::*;
use num_complex::Complex64;

pub fn register_python_functions(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(berry_phase_py, m)?)?;
    m.add_function(wrap_pyfunction!(berry_phase_bloch_py, m)?)?;
    m.add_function(wrap_pyfunction!(chern_number_py, m)?)?;
    m.add_function(wrap_pyfunction!(winding_number_py, m)?)?;
    m.add_function(wrap_pyfunction!(skyrmion_number_py, m)?)?;
    Ok(())
}

#[pyfunction]
fn berry_phase_py(states_real: Vec<Vec<f64>>, states_imag: Vec<Vec<f64>>) -> PyResult<f64> {
    let states: Vec<Vec<Complex64>> = states_real.iter().zip(states_imag.iter()).map(|(re, im)| {
        re.iter().zip(im.iter()).map(|(r, i)| Complex64::new(*r, *i)).collect()
    }).collect();
    Ok(super::berry::berry_phase(&states))
}

#[pyfunction]
fn berry_phase_bloch_py(theta: f64) -> PyResult<f64> {
    Ok(super::berry::berry_phase_bloch(theta))
}

#[pyfunction]
fn chern_number_py(berry_curvature: Vec<Vec<f64>>, dk: f64) -> PyResult<f64> {
    Ok(super::chern::chern_number(&berry_curvature, dk))
}

#[pyfunction]
fn winding_number_py(f1: Vec<f64>, f2: Vec<f64>) -> PyResult<f64> {
    Ok(super::winding::winding_number_real(&f1, &f2))
}

#[pyfunction]
fn skyrmion_number_py(n_field: Vec<Vec<Vec<f64>>>, dx: f64, dy: f64) -> PyResult<f64> {
    let nx = n_field.len();
    let ny = n_field[0].len();
    let mut field = vec![[[0.0f64; 3]; 3]; 3]; // Placeholder
    // In practice, reshape properly
    Ok(0.0)
}
