//! Core error types and result definitions for pysic-rs.

use thiserror::Error;

/// Custom error type for pysic operations.
#[derive(Error, Debug, Clone)]
pub enum PysicError {
    #[error("Invalid parameter: {0}")]
    InvalidParameter(String),

    #[error("Invalid input: {0}")]
    InvalidInput(String),

    #[error("Dimension mismatch: expected {expected}, got {actual}")]
    DimensionMismatch { expected: usize, actual: usize },

    #[error("Empty data provided")]
    EmptyData,

    #[error("Convergence failed after {iterations} iterations: {message}")]
    ConvergenceFailed {
        iterations: usize,
        message: String,
    },

    #[error("Numerical error: {0}")]
    NumericalError(String),

    #[error("Computation error: {0}")]
    ComputationError(String),

    #[error("Singular matrix")]
    SingularMatrix,

    #[error("Matrix must be square, got {rows}x{cols}")]
    NonSquareMatrix { rows: usize, cols: usize },

    #[error("Boundary condition error: {0}")]
    BoundaryError(String),

    #[error("CFL condition violated: dt={dt} > max_dt={max_dt}")]
    CflViolation { dt: f64, max_dt: f64 },

    #[error("Tensor index out of bounds: index {index} >= dimension {dim}")]
    IndexOutOfBounds { index: usize, dim: usize },

    #[error("Signature error: expected {expected}-dimensional Lorentzian, got signature {got}")]
    SignatureError { expected: String, got: String },
}

/// Result type for pysic operations.
pub type Result<T> = std::result::Result<T, PysicError>;
