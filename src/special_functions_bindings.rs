//! Python bindings for the special-function suite.
//!
//! These were measured in Rust but never reachable from Python, although the
//! documented binding surface advertises `ps.bessel_j0(x)` and `ps.beta(x, 2.0)`.
//! The semantics implemented here are the documented ones:
//!
//!   * a float returns a float, any sequence returns a list (shape-preserving);
//!   * two-argument functions broadcast scalar against sequence;
//!   * a length mismatch between two sequences **raises** rather than truncating.
//!
//! The last point matters: silently truncating to the shorter operand is the
//! kind of defect this library's validation discipline exists to catch.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::special_functions as sf;

/// Accepts either a scalar or any Python sequence of floats.
#[derive(FromPyObject)]
enum Arg {
    Scalar(f64),
    Seq(Vec<f64>),
}

/// Elementwise lift of a one-argument function.
fn map1(py: Python<'_>, x: Arg, f: impl Fn(f64) -> f64) -> PyResult<PyObject> {
    match x {
        Arg::Scalar(v) => Ok(f(v).into_py(py)),
        Arg::Seq(v) => Ok(v.into_iter().map(f).collect::<Vec<f64>>().into_py(py)),
    }
}

/// Elementwise lift of a two-argument function, with broadcasting.
fn map2(py: Python<'_>, a: Arg, b: Arg, f: impl Fn(f64, f64) -> f64) -> PyResult<PyObject> {
    match (a, b) {
        (Arg::Scalar(x), Arg::Scalar(y)) => Ok(f(x, y).into_py(py)),
        (Arg::Seq(xs), Arg::Scalar(y)) => {
            Ok(xs.into_iter().map(|x| f(x, y)).collect::<Vec<f64>>().into_py(py))
        }
        (Arg::Scalar(x), Arg::Seq(ys)) => {
            Ok(ys.into_iter().map(|y| f(x, y)).collect::<Vec<f64>>().into_py(py))
        }
        (Arg::Seq(xs), Arg::Seq(ys)) => {
            if xs.len() != ys.len() {
                return Err(PyValueError::new_err(format!(
                    "length mismatch: {} vs {} (pysic-rs broadcasts scalars, but will not \
                     truncate two sequences of different length)",
                    xs.len(),
                    ys.len()
                )));
            }
            Ok(xs
                .iter()
                .zip(ys.iter())
                .map(|(&x, &y)| f(x, y))
                .collect::<Vec<f64>>()
                .into_py(py))
        }
    }
}

pub fn register_python_functions(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(gamma_py, m)?)?;
    m.add_function(wrap_pyfunction!(ln_gamma_py, m)?)?;
    m.add_function(wrap_pyfunction!(beta_py, m)?)?;
    m.add_function(wrap_pyfunction!(erf_py, m)?)?;
    m.add_function(wrap_pyfunction!(erfc_py, m)?)?;
    m.add_function(wrap_pyfunction!(bessel_j0_py, m)?)?;
    m.add_function(wrap_pyfunction!(bessel_j1_py, m)?)?;
    m.add_function(wrap_pyfunction!(bessel_y0_py, m)?)?;
    m.add_function(wrap_pyfunction!(bessel_y1_py, m)?)?;
    m.add_function(wrap_pyfunction!(legendre_p_py, m)?)?;
    m.add_function(wrap_pyfunction!(legendre_plm_py, m)?)?;
    m.add_function(wrap_pyfunction!(spherical_harmonic_py, m)?)?;
    m.add_function(wrap_pyfunction!(chebyshev_t_py, m)?)?;
    m.add_function(wrap_pyfunction!(chebyshev_u_py, m)?)?;
    m.add_function(wrap_pyfunction!(airy_ai_py, m)?)?;
    m.add_function(wrap_pyfunction!(expint_e1_py, m)?)?;
    m.add_function(wrap_pyfunction!(expint_en_py, m)?)?;
    m.add_function(wrap_pyfunction!(zeta_py, m)?)?;
    Ok(())
}

// ── one-argument, elementwise ───────────────────────────────────────────

#[pyfunction]
fn gamma_py(py: Python<'_>, x: Arg) -> PyResult<PyObject> {
    map1(py, x, sf::gamma)
}

#[pyfunction]
fn ln_gamma_py(py: Python<'_>, x: Arg) -> PyResult<PyObject> {
    map1(py, x, sf::ln_gamma)
}

#[pyfunction]
fn erf_py(py: Python<'_>, x: Arg) -> PyResult<PyObject> {
    map1(py, x, sf::erf)
}

#[pyfunction]
fn erfc_py(py: Python<'_>, x: Arg) -> PyResult<PyObject> {
    map1(py, x, sf::erfc)
}

#[pyfunction]
fn bessel_j0_py(py: Python<'_>, x: Arg) -> PyResult<PyObject> {
    map1(py, x, sf::bessel_j0)
}

#[pyfunction]
fn bessel_j1_py(py: Python<'_>, x: Arg) -> PyResult<PyObject> {
    map1(py, x, sf::bessel_j1)
}

#[pyfunction]
fn bessel_y0_py(py: Python<'_>, x: Arg) -> PyResult<PyObject> {
    map1(py, x, sf::bessel_y0)
}

#[pyfunction]
fn bessel_y1_py(py: Python<'_>, x: Arg) -> PyResult<PyObject> {
    map1(py, x, sf::bessel_y1)
}

#[pyfunction]
fn airy_ai_py(py: Python<'_>, x: Arg) -> PyResult<PyObject> {
    map1(py, x, sf::airy_ai)
}

#[pyfunction]
fn expint_e1_py(py: Python<'_>, x: Arg) -> PyResult<PyObject> {
    map1(py, x, sf::expint_e1)
}

#[pyfunction]
fn zeta_py(py: Python<'_>, s: Arg) -> PyResult<PyObject> {
    map1(py, s, sf::zeta)
}

// ── two-argument, broadcasting ──────────────────────────────────────────

#[pyfunction]
fn beta_py(py: Python<'_>, a: Arg, b: Arg) -> PyResult<PyObject> {
    map2(py, a, b, sf::beta)
}

// ── integer-order families: order is a scalar, argument is elementwise ──

#[pyfunction]
fn legendre_p_py(py: Python<'_>, l: usize, x: Arg) -> PyResult<PyObject> {
    map1(py, x, |v| sf::legendre_p(l, v))
}

#[pyfunction]
fn legendre_plm_py(py: Python<'_>, l: usize, m: i32, x: Arg) -> PyResult<PyObject> {
    map1(py, x, |v| sf::legendre_plm(l, m, v))
}

#[pyfunction]
fn chebyshev_t_py(py: Python<'_>, n: usize, x: Arg) -> PyResult<PyObject> {
    map1(py, x, |v| sf::chebyshev_t(n, v))
}

#[pyfunction]
fn chebyshev_u_py(py: Python<'_>, n: usize, x: Arg) -> PyResult<PyObject> {
    map1(py, x, |v| sf::chebyshev_u(n, v))
}

#[pyfunction]
fn expint_en_py(py: Python<'_>, n: usize, x: Arg) -> PyResult<PyObject> {
    map1(py, x, |v| sf::expint_en(n, v))
}

/// Spherical harmonic `Y_l^m(theta, phi)`, returned as a `(re, im)` tuple.
#[pyfunction]
fn spherical_harmonic_py(l: usize, m: i32, theta: f64, phi: f64) -> PyResult<(f64, f64)> {
    let z = sf::spherical_harmonic(l, m, theta, phi);
    Ok((z.re, z.im))
}
