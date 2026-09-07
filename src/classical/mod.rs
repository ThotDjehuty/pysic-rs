//! Classical mechanics: Hamiltonian, Lagrangian, rigid body dynamics.

pub mod hamiltonian;
pub mod lagrangian;
pub mod rigid_body;

pub use hamiltonian::{hamiltonian_h, canonical_equations, poisson_bracket, canonical_transform};
pub use lagrangian::{lagrangian_l, euler_lagrange, kinetic_energy, potential_energy};
pub use rigid_body::{euler_equations, torque_free, inertia_tensor, angular_momentum};

#[cfg(feature = "python")]
pub mod python_bindings;
