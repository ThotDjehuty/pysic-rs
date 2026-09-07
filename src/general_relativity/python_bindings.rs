//! Python bindings for general relativity.

use pyo3::prelude::*;

pub fn register_python_functions(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(schwarzschild_metric_py, m)?)?;
    m.add_function(wrap_pyfunction!(christoffel_from_metric_py, m)?)?;
    m.add_function(wrap_pyfunction!(ricci_scalar_py, m)?)?;
    m.add_function(wrap_pyfunction!(einstein_tensor_py, m)?)?;
    m.add_function(wrap_pyfunction!(hamiltonian_constraint_py, m)?)?;
    m.add_function(wrap_pyfunction!(momentum_constraint_py, m)?)?;
    Ok(())
}

#[pyfunction]
fn schwarzschild_metric_py(r: f64, mass: f64) -> PyResult<Vec<Vec<f64>>> {
    let g = super::metrics::schwarzschild_metric(r, mass);
    Ok(g.rows().into_iter().map(|row| row.to_vec()).collect())
}

#[pyfunction]
fn christoffel_from_metric_py(
    coords: Vec<f64>,
    h: Option<f64>,
) -> PyResult<Vec<Vec<Vec<f64>>>> {
    // Use Schwarzschild as default metric
    let metric = |c: &[f64]| super::metrics::schwarzschild_metric(c[1], 1.0);
    Ok(super::christoffel::christoffel_from_metric(&metric, &coords, h.unwrap_or(1e-5)))
}

#[pyfunction]
fn ricci_scalar_py(coords: Vec<f64>, h: Option<f64>) -> PyResult<f64> {
    let metric = |c: &[f64]| super::metrics::schwarzschild_metric(c[1], 1.0);
    Ok(super::christoffel::ricci_scalar(&metric, &coords, h.unwrap_or(1e-4)))
}

#[pyfunction]
fn einstein_tensor_py(coords: Vec<f64>, h: Option<f64>) -> PyResult<Vec<Vec<f64>>> {
    let metric = |c: &[f64]| super::metrics::schwarzschild_metric(c[1], 1.0);
    let g = super::christoffel::einstein_tensor(&metric, &coords, h.unwrap_or(1e-4));
    Ok(g.rows().into_iter().map(|row| row.to_vec()).collect())
}

#[pyfunction]
fn hamiltonian_constraint_py(
    gamma: Vec<Vec<f64>>,
    k: Vec<Vec<f64>>,
    rho: f64,
) -> PyResult<f64> {
    use ndarray::Array2;
    let n = gamma.len();
    let gamma_arr = Array2::from_shape_fn((n, n), |(i, j)| gamma[i][j]);
    let k_arr = Array2::from_shape_fn((n, n), |(i, j)| k[i][j]);
    Ok(super::adm::hamiltonian_constraint(&gamma_arr, &k_arr, rho))
}

#[pyfunction]
fn momentum_constraint_py(
    gamma: Vec<Vec<f64>>,
    k: Vec<Vec<f64>>,
    j: Vec<f64>,
) -> PyResult<Vec<f64>> {
    use ndarray::Array2;
    let n = gamma.len();
    let gamma_arr = Array2::from_shape_fn((n, n), |(i, j)| gamma[i][j]);
    let k_arr = Array2::from_shape_fn((n, n), |(i, j)| k[i][j]);
    Ok(super::adm::momentum_constraint(&gamma_arr, &k_arr, &j))
}
