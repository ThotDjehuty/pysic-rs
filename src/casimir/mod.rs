//! Casimir effect: vacuum energy between conducting plates.

pub mod energy;

pub use energy::{casimir_energy_parallel_plates, casimir_energy_sphere, casimir_force};

#[cfg(feature = "python")]
pub mod python_bindings;
