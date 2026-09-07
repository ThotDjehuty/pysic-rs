//! Python bindings for gauge theory.

use pyo3::prelude::*;

pub fn register_python_functions(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(su2_structure_constants_py, m)?)?;
    m.add_function(wrap_pyfunction!(su3_structure_constants_py, m)?)?;
    m.add_function(wrap_pyfunction!(field_strength_point_charge_py, m)?)?;
    m.add_function(wrap_pyfunction!(instanton_action_py, m)?)?;
    Ok(())
}

#[pyfunction]
fn su2_structure_constants_py() -> PyResult<Vec<Vec<Vec<f64>>>> {
    let arr = super::bundles::su2_structure_constants();
    Ok(arr.iter().map(|a| a.iter().map(|b| b.to_vec()).collect()).collect())
}

#[pyfunction]
fn su3_structure_constants_py() -> PyResult<Vec<Vec<Vec<f64>>>> {
    let arr = super::bundles::su3_structure_constants();
    Ok(arr.iter().map(|a| a.iter().map(|b| b.to_vec()).collect()).collect())
}

#[pyfunction]
fn field_strength_point_charge_py(charge: f64, r: f64, eps0: f64) -> PyResult<f64> {
    Ok(super::curvature::field_strength_point_charge(charge, r, eps0))
}

#[pyfunction]
fn instanton_action_py(g: f64) -> PyResult<f64> {
    Ok(super::yang_mills::instanton_action(g))
}
