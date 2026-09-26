//! Python bindings for the special-functions suite.
//!
//! The Rust implementations in [`crate::special_functions`] were previously
//! unreachable from Python, which is why the documentation carried a "not yet
//! bound" note. Every entry point below is a thin wrapper; the numerics live in
//! the core module.
//!
//! Scalar functions accept either a float or a sequence, so they can be called
//! directly on a grid for plotting:
//!
//! ```python
//! import numpy as np, pysicrs as ps
//! x = np.linspace(0.1, 5, 400)
//! y = ps.gamma(x.tolist())        # -> list[float]
//! ```

use pyo3::prelude::*;

use crate::special_functions as sf;

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

/// Accept a float or any sequence of floats, returning the same shape.
fn map1(py: Python, arg: &Bound<'_, PyAny>, f: impl Fn(f64) -> f64) -> PyResult<PyObject> {
    if let Ok(x) = arg.extract::<f64>() {
        return Ok(f(x).into_py(py));
    }
    let xs: Vec<f64> = arg.extract()?;
    let ys: Vec<f64> = xs.into_iter().map(&f).collect();
    Ok(ys.into_py(py))
}

/// Gamma function Γ(x), via Lanczos approximation.
#[pyfunction]
#[pyo3(name = "gamma_py", signature = (x))]
fn gamma_py(py: Python, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, sf::gamma)
}

/// Natural logarithm of |Γ(x)|. Use for large arguments where Γ overflows.
#[pyfunction]
#[pyo3(name = "ln_gamma_py", signature = (x))]
fn ln_gamma_py(py: Python, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, sf::ln_gamma)
}

/// Beta function B(a, b) = Γ(a)Γ(b)/Γ(a+b).
#[pyfunction]
#[pyo3(name = "beta_py", signature = (a, b))]
fn beta_py(a: f64, b: f64) -> PyResult<f64> {
    Ok(sf::beta(a, b))
}

/// Error function erf(x) = (2/√π)∫₀ˣ e^{-t²} dt.
#[pyfunction]
#[pyo3(name = "erf_py", signature = (x))]
fn erf_py(py: Python, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, sf::erf)
}

/// Complementary error function erfc(x) = 1 - erf(x).
#[pyfunction]
#[pyo3(name = "erfc_py", signature = (x))]
fn erfc_py(py: Python, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, sf::erfc)
}

/// Bessel function of the first kind, order 0.
#[pyfunction]
#[pyo3(name = "bessel_j0_py", signature = (x))]
fn bessel_j0_py(py: Python, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, sf::bessel_j0)
}

/// Bessel function of the first kind, order 1.
#[pyfunction]
#[pyo3(name = "bessel_j1_py", signature = (x))]
fn bessel_j1_py(py: Python, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, sf::bessel_j1)
}

/// Bessel function of the second kind, order 0. Singular at x = 0.
#[pyfunction]
#[pyo3(name = "bessel_y0_py", signature = (x))]
fn bessel_y0_py(py: Python, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, sf::bessel_y0)
}

/// Bessel function of the second kind, order 1. Singular at x = 0.
#[pyfunction]
#[pyo3(name = "bessel_y1_py", signature = (x))]
fn bessel_y1_py(py: Python, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, sf::bessel_y1)
}

/// Legendre polynomial P_l(x), by the stable upward recurrence.
#[pyfunction]
#[pyo3(name = "legendre_p_py", signature = (l, x))]
fn legendre_p_py(py: Python, l: usize, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, |v| sf::legendre_p(l, v))
}

/// Associated Legendre function P_l^m(x).
#[pyfunction]
#[pyo3(name = "legendre_plm_py", signature = (l, m, x))]
fn legendre_plm_py(py: Python, l: usize, m: i32, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, |v| sf::legendre_plm(l, m, v))
}

/// Spherical harmonic Y_l^m(θ, φ), returned as the pair (re, im).
#[pyfunction]
#[pyo3(name = "spherical_harmonic_py", signature = (l, m, theta, phi))]
fn spherical_harmonic_py(l: usize, m: i32, theta: f64, phi: f64) -> PyResult<(f64, f64)> {
    let z = sf::spherical_harmonic(l, m, theta, phi);
    Ok((z.re, z.im))
}

/// Chebyshev polynomial of the first kind, T_n(x) = cos(n arccos x).
#[pyfunction]
#[pyo3(name = "chebyshev_t_py", signature = (n, x))]
fn chebyshev_t_py(py: Python, n: usize, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, |v| sf::chebyshev_t(n, v))
}

/// Chebyshev polynomial of the second kind, U_n(x).
#[pyfunction]
#[pyo3(name = "chebyshev_u_py", signature = (n, x))]
fn chebyshev_u_py(py: Python, n: usize, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, |v| sf::chebyshev_u(n, v))
}

/// Airy function Ai(x), the decaying solution of y'' = x y.
#[pyfunction]
#[pyo3(name = "airy_ai_py", signature = (x))]
fn airy_ai_py(py: Python, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, sf::airy_ai)
}

/// Exponential integral E₁(x) = ∫₁^∞ e^{-xt}/t dt, for x > 0.
#[pyfunction]
#[pyo3(name = "expint_e1_py", signature = (x))]
fn expint_e1_py(py: Python, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, sf::expint_e1)
}

/// Generalised exponential integral Eₙ(x).
#[pyfunction]
#[pyo3(name = "expint_en_py", signature = (n, x))]
fn expint_en_py(py: Python, n: usize, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, |v| sf::expint_en(n, v))
}

/// Riemann zeta function ζ(s) on the real line.
#[pyfunction]
#[pyo3(name = "zeta_py", signature = (s))]
fn zeta_py(py: Python, s: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, s, sf::zeta)
}
