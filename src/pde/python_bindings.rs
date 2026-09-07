//! Python bindings for PDE solvers.

use pyo3::prelude::*;
use num_complex::Complex64;

pub fn register_python_functions(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(schrodinger_split_step_1d_py, m)?)?;
    m.add_function(wrap_pyfunction!(schrodinger_eigen_1d_py, m)?)?;
    m.add_function(wrap_pyfunction!(dirac_split_step_1d_py, m)?)?;
    m.add_function(wrap_pyfunction!(heat_crank_nicolson_1d_py, m)?)?;
    m.add_function(wrap_pyfunction!(heat_crank_nicolson_2d_py, m)?)?;
    m.add_function(wrap_pyfunction!(wave_fdtd_1d_py, m)?)?;
    m.add_function(wrap_pyfunction!(wave_fdtd_2d_py, m)?)?;
    m.add_function(wrap_pyfunction!(poisson_fft_2d_py, m)?)?;
    Ok(())
}

#[pyfunction]
fn schrodinger_split_step_1d_py(
    py: Python,
    psi0_real: Vec<f64>,
    psi0_imag: Vec<f64>,
    v: Vec<f64>,
    dx: f64,
    dt: f64,
    n_steps: usize,
    mass: Option<f64>,
    hbar: Option<f64>,
) -> PyResult<(Vec<f64>, Vec<f64>)> {
    let psi0: Vec<Complex64> = psi0_real.iter().zip(psi0_imag.iter()).map(|(r, i)| Complex64::new(*r, *i)).collect();
    let result = super::schrodinger::schrodinger_split_step_1d(
        &psi0, &v, dx, dt, n_steps, mass.unwrap_or(1.0), hbar.unwrap_or(1.0),
    );
    let real: Vec<f64> = result.iter().map(|c| c.re).collect();
    let imag: Vec<f64> = result.iter().map(|c| c.im).collect();
    Ok((real, imag))
}

#[pyfunction]
fn schrodinger_eigen_1d_py(
    v: Vec<f64>,
    dx: f64,
    n_states: usize,
    mass: Option<f64>,
    hbar: Option<f64>,
) -> PyResult<(Vec<f64>, Vec<Vec<f64>>)> {
    Ok(super::schrodinger::schrodinger_eigen_1d(
        &v, dx, n_states, mass.unwrap_or(1.0), hbar.unwrap_or(1.0),
    ))
}

#[pyfunction]
fn dirac_split_step_1d_py(
    py: Python,
    psi0_real: Vec<f64>,
    psi0_imag: Vec<f64>,
    v: Vec<f64>,
    dx: f64,
    dt: f64,
    n_steps: usize,
    c_speed: Option<f64>,
    mass: Option<f64>,
    hbar: Option<f64>,
) -> PyResult<(Vec<f64>, Vec<f64>)> {
    let psi0: Vec<Complex64> = psi0_real.iter().zip(psi0_imag.iter()).map(|(r, i)| Complex64::new(*r, *i)).collect();
    let result = super::dirac::dirac_split_step_1d(
        &psi0, &v, dx, dt, n_steps, c_speed.unwrap_or(1.0), mass.unwrap_or(1.0), hbar.unwrap_or(1.0),
    );
    let real: Vec<f64> = result.iter().map(|c| c.re).collect();
    let imag: Vec<f64> = result.iter().map(|c| c.im).collect();
    Ok((real, imag))
}

#[pyfunction]
fn heat_crank_nicolson_1d_py(
    u0: Vec<f64>,
    dx: f64,
    dt: f64,
    alpha: f64,
    n_steps: usize,
) -> PyResult<Vec<f64>> {
    Ok(super::heat::heat_crank_nicolson_1d(&u0, dx, dt, alpha, n_steps))
}

#[pyfunction]
fn heat_crank_nicolson_2d_py(
    u0: Vec<f64>,
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    dt: f64,
    alpha: f64,
    n_steps: usize,
) -> PyResult<Vec<f64>> {
    Ok(super::heat::heat_crank_nicolson_2d(&u0, nx, ny, dx, dy, dt, alpha, n_steps))
}

#[pyfunction]
fn wave_fdtd_1d_py(
    u0: Vec<f64>,
    v0: Vec<f64>,
    dx: f64,
    dt: f64,
    v: f64,
    n_steps: usize,
) -> PyResult<Vec<f64>> {
    Ok(super::wave::wave_fdtd_1d(&u0, &v0, dx, dt, v, n_steps))
}

#[pyfunction]
fn wave_fdtd_2d_py(
    u0: Vec<f64>,
    v0: Vec<f64>,
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    dt: f64,
    v: f64,
    n_steps: usize,
) -> PyResult<Vec<f64>> {
    Ok(super::wave::wave_fdtd_2d(&u0, &v0, nx, ny, dx, dy, dt, v, n_steps))
}

#[pyfunction]
fn poisson_fft_2d_py(
    f: Vec<f64>,
    nx: usize,
    ny: usize,
) -> PyResult<Vec<f64>> {
    Ok(super::poisson::poisson_fft_2d(&f, nx, ny))
}
