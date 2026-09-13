"""Smoke tests for the pysicrs Python bindings.

These tests exercise the public Python API surface documented in the README
and the Sphinx documentation. They run against the built ``pysicrs`` module
(installed via ``maturin develop`` / ``maturin build``).
"""

import math

import pysicrs


def test_module_exports_core_functions():
    for name in (
        "schwarzschild_metric",
        "fft",
        "ifft",
        "casimir_energy_parallel_plates",
        "su3_structure_constants",
        "constants",
        "ricci_scalar",
    ):
        assert hasattr(pysicrs, name), f"missing pysicrs.{name}"


def test_version():
    assert pysicrs.__version__ == "0.1.0"


def test_schwarzschild_metric():
    g = pysicrs.schwarzschild_metric(r=10.0, mass=1.0)
    assert len(g) == 4
    assert all(len(row) == 4 for row in g)
    assert g[0][0] < 0.0 and g[1][1] > 0.0


def test_fft_sine_produces_flat_spectrum():
    x = [math.sin(2 * math.pi * k / 256) for k in range(256)]
    real, imag = pysicrs.fft(x, [0.0] * 256)
    assert len(real) == 256 and len(imag) == 256
    for value in real:
        assert abs(value) < 1e-8


def test_fft_ifft_roundtrip():
    x = [math.sin(2 * math.pi * k / 16) for k in range(16)]
    real, imag = pysicrs.fft(x, [0.0] * 16)
    back_re, back_im = pysicrs.ifft(real, imag)
    for original, recovered in zip(x, back_re):
        assert abs(original - recovered) < 1e-8


def test_casimir_energy_sign_and_units():
    energy = pysicrs.casimir_energy_parallel_plates(
        d=1e-6, hbar=1.054571817e-34, c=2.99792458e8
    )
    assert energy < 0.0
    assert abs(energy) < 1e-5


def test_constants_contains_physical_values():
    constants = pysicrs.constants()
    assert isinstance(constants, dict)
    assert len(constants) > 0