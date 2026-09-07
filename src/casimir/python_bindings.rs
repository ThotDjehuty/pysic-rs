//! Python bindings for Casimir effect.

use pyo3::prelude::*;

pub fn register_python_functions(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(casimir_energy_parallel_plates_py, m)?)?;
    m.add_function(wrap_pyfunction!(casimir_force_py, m)?)?;
    m.add_function(wrap_pyfunction!(polder_potential_py, m)?)?;
    Ok(())
}

#[pyfunction]
fn casimir_energy_parallel_plates_py(d: f64, hbar: f64, c: f64) -> PyResult<f64> {
    Ok(super::energy::casimir_energy_parallel_plates(d, hbar, c))
}

#[pyfunction]
fn casimir_force_py(d: f64, hbar: f64, c: f64) -> PyResult<f64> {
    Ok(super::energy::casimir_force(d, hbar, c))
}

#[pyfunction]
fn polder_potential_py(distance: f64, polarizability: f64, hbar: f64, c: f64, retarded: bool) -> PyResult<f64> {
    Ok(super::energy::polder_potential(distance, polarizability, hbar, c, retarded))
}
