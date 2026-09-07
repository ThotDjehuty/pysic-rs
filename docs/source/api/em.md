# API Reference: Electromagnetism

The Python surface exposes electromagnetism routines at top level
(`from pysicrs import larmor_formula`). See also
[Algorithms: Electromagnetism](../algorithms/em.md).

```python
import pysicrs
```

## Green's functions

### `green_fn_static`

```python
green_fn_static(r: float) -> float
```

Static 3D Green function $\frac{1}{4\pi r}$.

## Radiation formulas

### `larmor_formula`

```python
larmor_formula(charge: float, acceleration: float, eps0: float, c: float) -> float
```

$P = \frac{2}{3}\frac{q^2 a^2}{4\pi\varepsilon_0 c^3}$.

```python
e, eps0, c = 1.602176634e-19, 8.8541878128e-12, 2.99792458e8
print(larmor_formula(e, 1.0, eps0, c))   # ~3.44e-50 W for a = 1 m/s²
```

### `dipole_radiation`

```python
dipole_radiation(dipole_moment: float, omega: float, mu0: float, c: float) -> float
```

$P = \frac{\mu_0\,\omega^4\,p_0^2}{12\pi c^3}$.

### `field_strength_point_charge`

```python
field_strength_point_charge(charge: float, r: float, eps0: float) -> float
```

Coulomb field $E = \frac{q}{4\pi\varepsilon_0 r^2}$.

### `compton_wavelength_shift`

```python
compton_wavelength_shift(theta: float, h: float, m_e: float, c: float) -> float
```

Compton shift $\Delta\lambda = \frac{h}{m_e c}(1-\cos\theta)$, `theta` in radians.

## Rust-only (not yet bound)

| function | purpose |
|----------|---------|
| `green_fn_retarded` / `green_fn_advanced` | causal propagators |
| `dipole_angular_distribution` | $\frac{dP}{d\Omega}\propto\sin^2\theta$ |
| `magnetic_dipole_radiation` | $\frac{\mu_0 m^2\omega^4}{12\pi c^3}$ |
| `radiation_power_oscillating` | cosine-time dependence |