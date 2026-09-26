//! General relativity: metrics, Christoffel, Riemann, geodesics, ADM decomposition.

pub mod metrics;
pub mod christoffel;
pub mod geodesics;
pub mod adm;
pub mod stress_energy;

pub use metrics::{schwarzschild_metric, schwarzschild_metric_full, kerr_metric, flrw_metric, minkowski_metric};
pub use christoffel::{christoffel_from_metric, riemann_tensor, ricci_tensor, ricci_scalar, einstein_tensor};
pub use geodesics::{geodesic_equation_rhs, integrate_geodesic};
pub use adm::{MetricSlice, AdmMetric, evolve_metric, hamiltonian_constraint, momentum_constraint};
pub use stress_energy::{perfect_fluid_stress_energy, electromagnetic_stress_energy};

#[cfg(feature = "python")]
pub mod python_bindings;
