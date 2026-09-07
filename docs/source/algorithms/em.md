# Electromagnetism

Pysic-rs implements the classical results of electrodynamics — Green's functions and the
canonical radiation-formula family — in natural SI units with explicit physical constants.

> **Python binding status:** `green_fn_static`, `larmor_formula`, `dipole_radiation`,
> `field_strength_point_charge`, and `compton_wavelength_shift` are bound. The retarded/
> advanced Green's functions and the angular-distribution helpers are Rust-only. See the
> [API page](../api/em.md).

---

## Mathematical Foundations

### Green's Functions of the Helmholtz/Poisson Operator

The Green function $G(\mathbf r, \mathbf r')$ solves
$\nabla^2 G = -\delta(\mathbf r - \mathbf r')$ (or the time-dependent analog) and builds
solutions by convolution. **Proofs:* standard textbook derivations (Jackson, Griffiths).

**Static (3D Poisson):**

$$
G_{\text{static}}(\mathbf r) = \frac{1}{4\pi |\mathbf r - \mathbf r'|}
$$

**Retarded (causal):**

$$
G_{\text{ret}}(\mathbf r, t; \mathbf r', t') = \frac{\delta\!\left(t' - t + \frac{|\mathbf r - \mathbf r'|}{c}\right)}{4\pi |\mathbf r - \mathbf r'|}
$$

**Advanced:**

$$
G_{\text{adv}}(\mathbf r, t; \mathbf r', t') = \frac{\delta\!\left(t' - t - \frac{|\mathbf r - \mathbf r'|}{c}\right)}{4\pi |\mathbf r - \mathbf r'|}
$$

**Proven properties:** $G_\text{ret}$ vanishes for $t' < t$ (causality); the retarded minus
advanced combination is the commutator (Feynman) propagator; the static solution is the
$c\to\infty$ limit.

### Radiation

**Larmor formula** (power radiated by an accelerating point charge; **proven** from the
Poynting flux):

$$
P = \frac{2}{3}\frac{q^2 a^2}{4\pi\varepsilon_0 c^3}
$$

For an oscillating dipole $p(t) = p_0\cos(\omega t)$, the **dipole radiation power**:

$$
P = \frac{\mu_0\,\omega^4\,p_0^2}{12\pi c^3}
$$

**Proven** ($\omega^4$ dipole law, Rayleigh scattering tail). The **angular distribution**:

$$
\frac{dP}{d\Omega} = \frac{\mu_0\,\omega^4\,p_0^2}{32\pi^2 c^3}\,\sin^2\theta
$$

The $\sin^2\theta$ pattern is the classic dipole lobe.

**Magnetic dipole** analog:

$$
P = \frac{\mu_0\,m^2\,\omega^4}{12\pi c^3}
$$

---

## Routines

| Routine | Formula |
|---------|---------|
| `green_fn_static(r)` | $\frac{1}{4\pi r}$ |
| `green_fn_retarded(r, t)` | causal $\delta$-propagator |
| `green_fn_advanced(r, t)` | anti-causal |
| `larmor_formula(q, a, eps0, c)` | $\frac{2}{3}\frac{q^2a^2}{4\pi\varepsilon_0c^3}$ |
| `dipole_radiation(p0, omega, mu0, c)` | $\frac{\mu_0\omega^4p_0^2}{12\pi c^3}$ |
| `radiation_power_oscillating(p0, omega, ...)` | dipole with cosine time dependence |
| `dipole_angular_distribution(theta, p0, ...)` | $\frac{dP}{d\Omega}$ |
| `magnetic_dipole_radiation(m, omega, ...)` | $\frac{\mu_0m^2\omega^4}{12\pi c^3}$ |

---

## Usage Examples

### Larmor radiation of an accelerating electron

```python
from pysicrs import larmor_formula

e = 1.602176634e-19     # C
eps0 = 8.8541878128e-12 # F/m
c = 2.99792458e8        # m/s
m_e = 9.1093837015e-31  # kg

g_accel = 9.81
a = g_accel / 2         # some acceleration scale
print(f"P = {larmor_formula(e, a, eps0, c):.3e} W")
```

### Dipole radiation power scales as ω⁴

```python
from pysicrs import dipole_radiation

p0, mu0, c = 1e-29, 1.25663706212e-6, 2.99792458e8
for f in (1e8, 1e9):
    w = 2 * 3.141592653589793 * f
    print(f"f={f:.0e} Hz: P = {dipole_radiation(p0, w, mu0, c):.3e} W")
```

### Angular distribution has a sin²θ lobe

```python
from pysicrs import dipole_angular_distribution

for th in (0.0, 0.7854, 1.5708):          # 0°, 45°, 90°
    print(f"θ={th:.2f}: dP/dΩ = {dipole_angular_distribution(th, 1e-29, 2e9, mu0=1.25663706212e-6, c=2.99792458e8):.3e}")
# peaks at θ=90°, vanishes at θ=0° (along the dipole axis)
```

---

## Advantages & Limitations

✅ Exact constants in SI units — no ambiguity between unit systems

✅ Retarded/advanced propagation naturally causal

✅ Cover the textbook radiation `ω⁴` laws used in antenna and scattering theory

❌ Point charges only (finite-size and multipole corrections modeled by extra terms)

❌ No numerical solving of full wave/curl equations here — that's the PDE Maxwell FDTD

❌ No dielectric/magnetization response functions; those live in CASIMIR's Lifshitz machinery

---

## References

1. Jackson, J.D. (1999). *Classical Electrodynamics*, 3rd ed. Wiley.
2. Griffiths, D.J. (2017). *Introduction to Electrodynamics*, 4th ed. Cambridge.
3. Landau, L.D. & Lifshitz, E.M. (1975). *The Classical Theory of Fields*, 4th ed.
4. Zangwill, A. (2012). *Modern Electrodynamics*. Cambridge.

---

## Related Topics

- [PDE](pde.md) – Maxwell 3D FDTD complements the analytic formulas here
- [Casimir](casimir.md) – vacuum energy is a mode sum of EM modes between plates
- [Linalg](linalg.md) – field tensors $F_{\mu\nu}$ interactive with raise/lower