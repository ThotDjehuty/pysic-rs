//! Python bindings for ODE solvers.

use pyo3::prelude::*;
use pyo3::types::PyList;

pub fn register_python_functions(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(rk4_solve_py, m)?)?;
    m.add_function(wrap_pyfunction!(rk45_solve_py, m)?)?;
    m.add_function(wrap_pyfunction!(backward_euler_py, m)?)?;
    m.add_function(wrap_pyfunction!(leapfrog_integrate_py, m)?)?;
    Ok(())
}

#[pyfunction]
fn rk4_solve_py(
    py: Python,
    f: PyObject,
    y0: Vec<f64>,
    t_start: f64,
    t_end: f64,
    n_steps: usize,
) -> PyResult<(Vec<f64>, Vec<Vec<f64>>)> {
    let f_ref = |t: f64, y: &[f64]| -> Vec<f64> {
        Python::with_gil(|py| {
            let args = (t, y.to_vec());
            let result = f.call1(py, args).unwrap();
            result.extract::<Vec<f64>>(py).unwrap()
        })
    };
    Ok(super::rk4_solve(&f_ref, &y0, (t_start, t_end), n_steps))
}

#[pyfunction]
fn rk45_solve_py(
    py: Python,
    f: PyObject,
    y0: Vec<f64>,
    t_start: f64,
    t_end: f64,
    rtol: Option<f64>,
    atol: Option<f64>,
    max_steps: Option<usize>,
) -> PyResult<PyObject> {
    let f_ref = |t: f64, y: &[f64]| -> Vec<f64> {
        Python::with_gil(|py| {
            let args = (t, y.to_vec());
            let result = f.call1(py, args).unwrap();
            result.extract::<Vec<f64>>(py).unwrap()
        })
    };
    let sol = super::rk45_solve(
        &f_ref,
        &y0,
        (t_start, t_end),
        rtol.unwrap_or(1e-6),
        atol.unwrap_or(1e-9),
        max_steps.unwrap_or(100_000),
    );
    Ok(py.None()) // TODO: return a proper struct
}

#[pyfunction]
fn backward_euler_py(
    py: Python,
    f: PyObject,
    y0: Vec<f64>,
    t_start: f64,
    t_end: f64,
    n_steps: usize,
) -> PyResult<(Vec<f64>, Vec<Vec<f64>>)> {
    let f_ref = |t: f64, y: &[f64]| -> Vec<f64> {
        Python::with_gil(|py| {
            let args = (t, y.to_vec());
            let result = f.call1(py, args).unwrap();
            result.extract::<Vec<f64>>(py).unwrap()
        })
    };
    Ok(super::backward_euler(&f_ref, &y0, (t_start, t_end), n_steps, 20, 1e-10))
}

#[pyfunction]
fn leapfrog_integrate_py(
    py: Python,
    grad_v: PyObject,
    grad_t: PyObject,
    q0: Vec<f64>,
    p0: Vec<f64>,
    dt: f64,
    n_steps: usize,
    inv_mass: f64,
) -> PyResult<(Vec<Vec<f64>>, Vec<Vec<f64>>)> {
    let gv_ref = |q: &[f64]| -> Vec<f64> {
        Python::with_gil(|py| {
            let args = (q.to_vec(),);
            let result = grad_v.call1(py, args).unwrap();
            result.extract::<Vec<f64>>(py).unwrap()
        })
    };
    let gt_ref = |p: &[f64]| -> Vec<f64> {
        Python::with_gil(|py| {
            let args = (p.to_vec(),);
            let result = grad_t.call1(py, args).unwrap();
            result.extract::<Vec<f64>>(py).unwrap()
        })
    };
    Ok(super::integrate_hamiltonian(&q0, &p0, &gv_ref, &gt_ref, dt, n_steps, inv_mass))
}
