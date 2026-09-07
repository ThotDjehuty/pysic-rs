//! Pysic-rs — High-Performance Mathematical Physics Engine
//! =======================================================
//!
//! Provides CPU-only, generic implementations of standard mathematical physics:
//!
//! - **Special functions**: Bessel, Legendre, Airy, spherical harmonics, etc.
//! - **Linear algebra**: Cholesky, SVD, eigenvalues, tensor operations
//! - **Numerical calculus**: differentiation, integration, interpolation
//! - **ODE solvers**: RK4, RK45 (Dormand-Prince), implicit, symplectic (Leapfrog, Verlet)
//! - **PDE solvers**: Schrödinger, Dirac, heat, wave, Poisson, Maxwell FDTD
//! - **Fourier analysis**: FFT, spectral derivatives, power spectral density
//! - **Quantum mechanics**: Pauli matrices, density matrices, path integrals, propagators
//! - **General relativity**: metrics, Christoffel, Riemann, geodesics, ADM decomposition
//! - **Gauge theory**: connections, curvature, Yang-Mills, Lie algebras
//! - **Classical mechanics**: Hamiltonian/Lagrangian, rigid body, symplectic integration
//! - **Topology**: Berry phase, Chern numbers, winding numbers
//! - **Electromagnetism**: Green's functions, multipole expansion, radiation
//! - **Casimir effect**: energy and pressure between conductive plates
//!
//! All implementations are validated textbook physics. No experimental or
//! project-specific models are included.

#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3::types::PyModule;

pub mod core;
pub mod constants;

// Foundational modules
pub mod special_functions;
pub mod linalg;
pub mod calculus;

// Physics modules
pub mod ode;
pub mod pde;
pub mod fourier;
pub mod quantum;
pub mod general_relativity;
pub mod gauge;
pub mod classical;
pub mod topology;
pub mod em;
pub mod casimir;

/// Pysic-rs Python module
#[cfg(feature = "python")]
#[pymodule]
fn _core(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Constants
    m.add_function(wrap_pyfunction!(constants_py, m)?)?;

    // ODE
    ode::python_bindings::register_python_functions(m)?;

    // PDE
    pde::python_bindings::register_python_functions(m)?;

    // Fourier
    fourier::python_bindings::register_python_functions(m)?;

    // Quantum
    quantum::python_bindings::register_python_functions(m)?;

    // General Relativity
    general_relativity::python_bindings::register_python_functions(m)?;

    // Gauge
    gauge::python_bindings::register_python_functions(m)?;

    // Classical
    classical::python_bindings::register_python_functions(m)?;

    // Topology
    topology::python_bindings::register_python_functions(m)?;

    // EM
    em::python_bindings::register_python_functions(m)?;

    // Casimir
    casimir::python_bindings::register_python_functions(m)?;

    Ok(())
}

/// Return a dict of all physical constants.
#[cfg(feature = "python")]
#[pyfunction]
fn constants_py(_py: Python) -> std::collections::HashMap<&'static str, f64> {
    let mut m = std::collections::HashMap::new();
    m.insert("c", constants::c());
    m.insert("hbar", constants::hbar());
    m.insert("h", constants::h());
    m.insert("G", constants::G());
    m.insert("k_B", constants::k_B());
    m.insert("e", constants::e());
    m.insert("m_e", constants::m_e());
    m.insert("m_p", constants::m_p());
    m.insert("eps_0", constants::eps_0());
    m.insert("mu_0", constants::mu_0());
    m.insert("alpha", constants::alpha());
    m.insert("sigma_SB", constants::sigma_SB());
    m.insert("a_0", constants::a_0());
    m.insert("lambda_C", constants::lambda_C());
    m.insert("R_inf", constants::R_inf());
    m
}
