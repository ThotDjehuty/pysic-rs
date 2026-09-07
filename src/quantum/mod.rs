//! Quantum mechanics: Pauli matrices, density matrices, path integrals, propagators.

pub mod pauli;
pub mod hilbert;
pub mod perturbation;
pub mod path_integral;
pub mod propagator;

pub use pauli::{pauli_matrices, spin_operators, density_matrix, partial_trace, purity};
pub use hilbert::{inner_product, projector, coherent_state, number_state};
pub use perturbation::{rayleigh_schrodinger_1st, rayleigh_schrodinger_2nd};
pub use path_integral::{discretized_path_integral, effective_potential, wick_rotation};
pub use propagator::{feynman_propagator_scalar, free_propagator_momentum, free_propagator_position};

#[cfg(feature = "python")]
pub mod python_bindings;
