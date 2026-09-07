//! Python bindings for classical mechanics.

use pyo3::prelude::*;

pub fn register_python_functions(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(euler_equations_py, m)?)?;
    m.add_function(wrap_pyfunction!(inertia_tensor_py, m)?)?;
    m.add_function(wrap_pyfunction!(euler_to_rotation_py, m)?)?;
    Ok(())
}

#[pyfunction]
fn euler_equations_py(
    omega: Vec<f64>,
    moments: Vec<f64>,
    torque: Vec<f64>,
) -> PyResult<Vec<f64>> {
    Ok(super::rigid_body::euler_equations(&omega, &moments, &torque))
}

#[pyfunction]
fn inertia_tensor_py(
    masses: Vec<f64>,
    positions: Vec<Vec<f64>>,
) -> PyResult<Vec<Vec<f64>>> {
    let pos: Vec<[f64; 3]> = positions.iter().map(|p| [p[0], p[1], p[2]]).collect();
    let i = super::rigid_body::inertia_tensor(&masses, &pos);
    Ok(i.iter().map(|row| row.to_vec()).collect())
}

#[pyfunction]
fn euler_to_rotation_py(phi: f64, theta: f64, psi: f64) -> PyResult<Vec<Vec<f64>>> {
    let r = super::rigid_body::euler_to_rotation(phi, theta, psi);
    Ok(r.iter().map(|row| row.to_vec()).collect())
}
