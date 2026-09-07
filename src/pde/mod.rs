//! PDE solvers for fundamental physics equations.

pub mod schrodinger;
pub mod dirac;
pub mod heat;
pub mod wave;
pub mod poisson;
pub mod maxwell_fdtd;

pub use schrodinger::{schrodinger_split_step_1d, schrodinger_eigen_1d};
pub use dirac::dirac_split_step_1d;
pub use heat::{heat_crank_nicolson_1d, heat_crank_nicolson_2d};
pub use wave::{wave_fdtd_1d, wave_fdtd_2d};
pub use poisson::poisson_fft_2d;

#[cfg(feature = "python")]
pub mod python_bindings;
