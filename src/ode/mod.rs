//! ODE solvers: RK4, RK45 (Dormand-Prince), implicit, symplectic.

pub mod rk4;
pub mod rk45;
pub mod implicit;
pub mod symplectic;

pub use rk4::rk4_solve;
pub use rk45::{rk45_solve, OdeSolution};
pub use implicit::{backward_euler, crank_nicolson_ode};
pub use symplectic::{leapfrog_step, velocity_verlet_step, yoshida_step, integrate_hamiltonian};

#[cfg(feature = "python")]
pub mod python_bindings;
