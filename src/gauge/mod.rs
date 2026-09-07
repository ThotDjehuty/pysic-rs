//! Gauge theory: connections, field strength, Yang-Mills, Lie algebras.

pub mod connection;
pub mod curvature;
pub mod yang_mills;
pub mod bundles;

pub use connection::{gauge_connection, covariant_derivative};
pub use curvature::field_strength;
pub use yang_mills::{yang_mills_action, yang_mills_eom};
pub use bundles::{su2_structure_constants, su3_structure_constants, structure_constants};

#[cfg(feature = "python")]
pub mod python_bindings;
