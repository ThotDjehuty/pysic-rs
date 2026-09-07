//! Python bindings for quantum mechanics.

use pyo3::prelude::*;
use pyo3::types::PyList;
use num_complex::Complex64;

pub fn register_python_functions(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pauli_matrices_py, m)?)?;
    m.add_function(wrap_pyfunction!(density_matrix_py, m)?)?;
    m.add_function(wrap_pyfunction!(coherent_state_py, m)?)?;
    m.add_function(wrap_pyfunction!(number_state_py, m)?)?;
    m.add_function(wrap_pyfunction!(feynman_propagator_scalar_py, m)?)?;
    m.add_function(wrap_pyfunction!(free_propagator_position_py, m)?)?;
    Ok(())
}

#[pyfunction]
fn pauli_matrices_py(py: Python) -> PyResult<PyObject> {
    let (sx, sy, sz) = super::pauli::pauli_matrices();
    let to_list = |m: &[[Complex64; 2]; 2]| -> Vec<Vec<f64>> {
        m.iter().map(|row| row.iter().map(|c| c.re).collect()).collect()
    };
    let result = (to_list(&sx), to_list(&sy), to_list(&sz));
    Ok(result.to_object(py))
}

#[pyfunction]
fn density_matrix_py(psi_real: Vec<f64>, psi_imag: Vec<f64>) -> PyResult<Vec<Vec<f64>>> {
    let psi: Vec<Complex64> = psi_real.iter().zip(psi_imag.iter()).map(|(r, i)| Complex64::new(*r, *i)).collect();
    let rho = super::pauli::density_matrix(&psi);
    Ok(rho.iter().map(|row| row.iter().map(|c| c.re).collect()).collect())
}

#[pyfunction]
fn coherent_state_py(alpha_real: f64, alpha_imag: f64, n_terms: usize) -> PyResult<(Vec<f64>, Vec<f64>)> {
    let alpha = Complex64::new(alpha_real, alpha_imag);
    let state = super::hilbert::coherent_state(alpha, n_terms);
    Ok((state.iter().map(|c| c.re).collect(), state.iter().map(|c| c.im).collect()))
}

#[pyfunction]
fn number_state_py(n: usize, dim: usize) -> PyResult<(Vec<f64>, Vec<f64>)> {
    let state = super::hilbert::number_state(n, dim);
    Ok((state.iter().map(|c| c.re).collect(), state.iter().map(|c| c.im).collect()))
}

#[pyfunction]
fn feynman_propagator_scalar_py(p_squared: f64, mass: f64, epsilon: Option<f64>) -> PyResult<(f64, f64)> {
    let g = super::propagator::feynman_propagator_scalar(p_squared, mass, epsilon.unwrap_or(1e-6));
    Ok((g.re, g.im))
}

#[pyfunction]
fn free_propagator_position_py(r: f64, mass: f64) -> PyResult<f64> {
    Ok(super::propagator::free_propagator_position(r, mass))
}
