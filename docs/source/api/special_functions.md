# API Reference: Special Functions

The Python surface exposes the special functions at top level (`from pysicrs import gamma`).
Some functions exist in the Rust core but are not yet bound from Python; bindings are added
module by module. See also [Algorithms: Special Functions](../algorithms/special_functions.md).

```python
import pysicrs
```

## Constants

```python
constants() -> dict[str, float]
```

Returns measured physical constants as a dictionary:

| key | value | symbol |
|-----|-------|--------|
| `c` | 2.997925e+08 m/s | speed of light |
| `hbar` | 1.054572e-34 J·s | reduced Planck |
| `h` | 6.626070e-34 J·s | Planck |
| `G` | 6.674e-11 m³ kg⁻¹ s⁻² | gravitational constant |
| `k_B` | 1.380649e-23 J/K | Boltzmann |
| `e` | 1.602177e-19 C | elementary charge |
| `m_e` | 9.109534e-31 kg | electron mass |
| `m_p` | 1.672622e-27 kg | proton mass |
| `eps_0` | 8.854188e-12 F/m | vacuum permittivity |
| `mu_0` | 1.256637e-06 N/A² | vacuum permeability |
| `alpha` | 7.297353e-03 | fine-structure constant |
| `sigma_SB` | 5.670374e-08 W m⁻² K⁻⁴ | Stefan–Boltzmann |
| `a_0` | 5.291772e-11 m | Bohr radius |
| `lambda_C` | 2.426310e-12 m | Compton wavelength |
| `R_inf` | 1.097373e+07 m⁻¹ | Rydberg constant |

## Bound exports (linalg / special-functions family)

The currently Python-bound entry points in this family are convenience helpers; the full
Rust surface (Gamma, Bessel, Legendre, erf, ζ, …) is under active binding. Every bound
routine is validated against its closed-form reference value in the test suite.

```python
compton_wavelength_shift(theta, h, m_e, c) -> float
```

Compton shift $\Delta\lambda = \frac{h}{m_e c}(1-\cos\theta)$, with `theta` in radians.

---

## Rust (not yet Python-bound) — see the algorithms page

The following Rust-exported functions are described on the
[Special Functions algorithms page](../algorithms/special_functions.md): `gamma`,
`ln_gamma`, `beta`, `bessel_j0/j1`, `bessel_y0/y1`, `legendre_p`, `legendre_plm`,
`spherical_harmonic`, `chebyshev_t/u`, `airy_ai`, `erf`, `erfc`, `expint_e1`, `expint_en`,
`zeta`.