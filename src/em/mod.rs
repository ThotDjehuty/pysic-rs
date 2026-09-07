//! Electromagnetism: Green's functions, radiation, classical E&M.

pub mod green_fn;
pub mod radiation;

pub use green_fn::{green_fn_static, green_fn_retarded, green_fn_advanced};
pub use radiation::{dipole_radiation, larmor_formula};

#[cfg(feature = "python")]
pub mod python_bindings;
