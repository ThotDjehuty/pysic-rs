//! Python bindings for electromagnetism.

use pyo3::prelude::*;

pub fn register_python_functions(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(green_fn_static_py, m)?)?;
    m.add_function(wrap_pyfunction!(larmor_formula_py, m)?)?;
    m.add_function(wrap_pyfunction!(dipole_radiation_py, m)?)?;
    m.add_function(wrap_pyfunction!(dipole_angular_distribution_py, m)?)?;
    m.add_function(wrap_pyfunction!(compton_wavelength_shift_py, m)?)?;
    Ok(())
}

#[pyfunction]
fn green_fn_static_py(r: f64) -> PyResult<f64> {
    Ok(super::green_fn::green_fn_static(r))
}

#[pyfunction]
fn larmor_formula_py(charge: f64, acceleration: f64, eps0: f64, c: f64) -> PyResult<f64> {
    Ok(super::radiation::larmor_formula(charge, acceleration, eps0, c))
}

#[pyfunction]
fn dipole_radiation_py(dipole_moment: f64, omega: f64, mu0: f64, c: f64) -> PyResult<f64> {
    Ok(super::radiation::dipole_radiation(dipole_moment, omega, mu0, c))
}

/// Angular distribution dP/dΩ = μ₀ω⁴|p|² sin²θ /(32π² c).
///
/// Integrating over the sphere reproduces `dipole_radiation`, since
/// ∫sin²θ dΩ = 8π/3.
#[pyfunction]
fn dipole_angular_distribution_py(
    theta: f64,
    dipole_moment: f64,
    omega: f64,
    mu0: f64,
    c: f64,
) -> PyResult<f64> {
    Ok(super::radiation::dipole_angular_distribution(
        theta, dipole_moment, omega, mu0, c,
    ))
}

#[pyfunction]
fn compton_wavelength_shift_py(theta: f64, h: f64, m_e: f64, c: f64) -> PyResult<f64> {
    Ok(super::radiation::compton_wavelength_shift(theta, h, m_e, c))
}
