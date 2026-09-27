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

/// Skyrmion number Q = (1/4π) ∫ n̂ · (∂_x n̂ × ∂_y n̂) d²x.
///
/// `n_field` is indexed `[ix][iy][component]` with exactly 3 components per
/// grid point. Central differences are used, so the grid must be at least
/// 3 × 3.
#[pyfunction]
fn skyrmion_number_py(n_field: Vec<Vec<Vec<f64>>>, dx: f64, dy: f64) -> PyResult<f64> {
    use pyo3::exceptions::PyValueError;

    let nx = n_field.len();
    if nx < 3 {
        return Err(PyValueError::new_err(format!(
            "n_field must have at least 3 rows for central differences, got {nx}"
        )));
    }
    let ny = n_field[0].len();
    if ny < 3 {
        return Err(PyValueError::new_err(format!(
            "n_field must have at least 3 columns for central differences, got {ny}"
        )));
    }

    let mut field: Vec<Vec<[f64; 3]>> = Vec::with_capacity(nx);
    for (ix, row) in n_field.iter().enumerate() {
        if row.len() != ny {
            return Err(PyValueError::new_err(format!(
                "n_field is ragged: row {ix} has {} entries, expected {ny}",
                row.len()
            )));
        }
        let mut out_row: Vec<[f64; 3]> = Vec::with_capacity(ny);
        for (iy, v) in row.iter().enumerate() {
            if v.len() != 3 {
                return Err(PyValueError::new_err(format!(
                    "n_field[{ix}][{iy}] has {} components, expected 3",
                    v.len()
                )));
            }
            out_row.push([v[0], v[1], v[2]]);
        }
        field.push(out_row);
    }

    Ok(super::winding::skyrmion_number(&field, dx, dy))
}
