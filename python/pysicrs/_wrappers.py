"""High-level Python wrappers for pysic-rs with NumPy array support."""

import numpy as np
from pysicrs._core import (
    constants as _constants,
    gamma as _gamma,
    erf as _erf,
    erfc as _erfc,
    cholesky as _cholesky,
    inverse as _inverse,
    determinant as _determinant,
    fft as _fft,
    ifft as _ifft,
)


def constants() -> dict:
    """Return all physical constants as a dictionary."""
    return _constants()


def gamma(x):
    """Gamma function Γ(x)."""
    return _gamma(x)


def erf(x):
    """Error function erf(x)."""
    return _erf(x)


def erfc(x):
    """Complementary error function erfc(x)."""
    return _erfc(x)


def cholesky(a):
    """Cholesky decomposition A = L L^T."""
    return np.array(_cholesky(a.tolist()))


def inverse(a):
    """Matrix inverse via Gauss-Jordan elimination."""
    return np.array(_inverse(a.tolist()))


def determinant(a):
    """Matrix determinant via LU decomposition."""
    return _determinant(a.tolist())


def fft(x):
    """Complex FFT."""
    real = x.real.tolist() if np.iscomplexobj(x) else x.tolist()
    imag = x.imag.tolist() if np.iscomplexobj(x) else [0.0] * len(x)
    re, im = _fft(real, imag)
    return np.array(re) + 1j * np.array(im)


def ifft(x):
    """Inverse complex FFT."""
    real = x.real.tolist()
    imag = x.imag.tolist()
    re, im = _ifft(real, imag)
    return np.array(re) + 1j * np.array(im)
