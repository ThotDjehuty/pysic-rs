//! Python bindings for the special-functions suite.
//!
//! The Rust implementations in [`crate::special_functions`] were previously
//! unreachable from Python, which is why the documentation carried a "not yet
//! bound" note. Every entry point below is a thin wrapper; the numerics live in
//! the core module.
//!
//! All of them are elementwise: pass a float and get a float back, pass any
//! sequence (list, tuple, NumPy array) and get a list back, so a function can
//! be evaluated straight onto a plotting grid.
//!
//! ```python
//! import numpy as np, pysicrs as ps
//! x = np.linspace(0.1, 25, 800)
//! ps.bessel_j0(x)      # -> list[float]
//! ps.bessel_j0(2.5)    # -> float
//! ```
//!
//! Two-argument functions broadcast NumPy-style: either argument may be scalar,
//! and two sequences must be the same length.

use pyo3::exceptions::PyValueError;
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

/// A float argument that may also arrive as a sequence.
enum Arg {
    Scalar(f64),
    Vector(Vec<f64>),
}

impl Arg {
    fn from(obj: &Bound<'_, PyAny>) -> PyResult<Self> {
        if let Ok(x) = obj.extract::<f64>() {
            return Ok(Arg::Scalar(x));
        }
        Ok(Arg::Vector(obj.extract()?))
    }

    fn len(&self) -> Option<usize> {
        match self {
            Arg::Scalar(_) => None,
            Arg::Vector(v) => Some(v.len()),
        }
    }

    fn at(&self, i: usize) -> f64 {
        match self {
            Arg::Scalar(x) => *x,
            Arg::Vector(v) => v[i],
        }
    }
}

/// Apply `f` elementwise, preserving scalar-in / scalar-out.
fn map1(py: Python, arg: &Bound<'_, PyAny>, f: impl Fn(f64) -> f64) -> PyResult<PyObject> {
    match Arg::from(arg)? {
        Arg::Scalar(x) => Ok(f(x).into_py(py)),
        Arg::Vector(v) => {
            let out: Vec<f64> = v.into_iter().map(&f).collect();
            Ok(out.into_py(py))
        }
    }
}

/// Apply `f` elementwise over two arguments, broadcasting a scalar against a
/// sequence. Two sequences must agree in length.
fn map2(
    py: Python,
    a: &Bound<'_, PyAny>,
    b: &Bound<'_, PyAny>,
    f: impl Fn(f64, f64) -> f64,
) -> PyResult<PyObject> {
    let (a, b) = (Arg::from(a)?, Arg::from(b)?);
    let n = match (a.len(), b.len()) {
        (None, None) => return Ok(f(a.at(0), b.at(0)).into_py(py)),
        (Some(n), None) | (None, Some(n)) => n,
        (Some(n), Some(m)) if n == m => n,
        (Some(n), Some(m)) => {
            return Err(PyValueError::new_err(format!(
                "length mismatch: {n} vs {m}"
            )))
        }
    };
    let out: Vec<f64> = (0..n).map(|i| f(a.at(i), b.at(i))).collect();
    Ok(out.into_py(py))
}

/// Gamma function Γ(x), via the Lanczos approximation.
///
/// Poles at the non-positive integers; use `ln_gamma` for large arguments
/// where Γ overflows.
#[pyfunction]
fn gamma_py(py: Python, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, sf::gamma)
}

/// Natural logarithm of |Γ(x)|, finite where Γ itself overflows.
#[pyfunction]
fn ln_gamma_py(py: Python, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, sf::ln_gamma)
}

/// Beta function B(a, b) = Γ(a)Γ(b)/Γ(a+b).
///
/// Either argument may be a sequence; a scalar broadcasts against it.
#[pyfunction]
fn beta_py(py: Python, a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map2(py, a, b, sf::beta)
}

/// Error function erf(x) = (2/√π)∫₀ˣ e^{-t²} dt.
#[pyfunction]
fn erf_py(py: Python, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, sf::erf)
}

/// Complementary error function erfc(x) = 1 − erf(x).
#[pyfunction]
fn erfc_py(py: Python, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, sf::erfc)
}

/// Bessel function of the first kind, order 0.
#[pyfunction]
fn bessel_j0_py(py: Python, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, sf::bessel_j0)
}

/// Bessel function of the first kind, order 1.
#[pyfunction]
fn bessel_j1_py(py: Python, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, sf::bessel_j1)
}

/// Bessel function of the second kind, order 0. NaN for x ≤ 0.
#[pyfunction]
fn bessel_y0_py(py: Python, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, sf::bessel_y0)
}

/// Bessel function of the second kind, order 1. NaN for x ≤ 0.
#[pyfunction]
fn bessel_y1_py(py: Python, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, sf::bessel_y1)
}

/// Legendre polynomial P_l(x), by the stable upward recurrence.
#[pyfunction]
fn legendre_p_py(py: Python, l: usize, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, |v| sf::legendre_p(l, v))
}

/// Associated Legendre function P_l^m(x).
#[pyfunction]
fn legendre_plm_py(py: Python, l: usize, m: i32, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, |v| sf::legendre_plm(l, m, v))
}

/// Spherical harmonic Y_l^m(θ, φ), returned as the pair `(re, im)`.
///
/// `theta` may be a sequence, in which case both halves of the pair are lists
/// of the same length.
#[pyfunction]
fn spherical_harmonic_py(
    py: Python,
    l: usize,
    m: i32,
    theta: &Bound<'_, PyAny>,
    phi: f64,
) -> PyResult<PyObject> {
    match Arg::from(theta)? {
        Arg::Scalar(t) => {
            let z = sf::spherical_harmonic(l, m, t, phi);
            Ok((z.re, z.im).into_py(py))
        }
        Arg::Vector(ts) => {
            let (re, im): (Vec<f64>, Vec<f64>) = ts
                .into_iter()
                .map(|t| {
                    let z = sf::spherical_harmonic(l, m, t, phi);
                    (z.re, z.im)
                })
                .unzip();
            Ok((re, im).into_py(py))
        }
    }
}

/// Chebyshev polynomial of the first kind, T_n(x) = cos(n arccos x).
#[pyfunction]
fn chebyshev_t_py(py: Python, n: usize, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, |v| sf::chebyshev_t(n, v))
}

/// Chebyshev polynomial of the second kind, U_n(x).
#[pyfunction]
fn chebyshev_u_py(py: Python, n: usize, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, |v| sf::chebyshev_u(n, v))
}

/// Airy function Ai(x), the decaying solution of y″ = x y.
#[pyfunction]
fn airy_ai_py(py: Python, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, sf::airy_ai)
}

/// Exponential integral E₁(x) = ∫₁^∞ e^{−xt}/t dt, for x > 0.
#[pyfunction]
fn expint_e1_py(py: Python, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, sf::expint_e1)
}

/// Generalised exponential integral Eₙ(x).
#[pyfunction]
fn expint_en_py(py: Python, n: usize, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, x, |v| sf::expint_en(n, v))
}

/// Riemann zeta function ζ(s) on the real line. Pole at s = 1.
#[pyfunction]
fn zeta_py(py: Python, s: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    map1(py, s, sf::zeta)
}
