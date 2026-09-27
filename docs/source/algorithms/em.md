# Electromagnetism

Pysic-rs implements the classical results of electrodynamics — Green's functions and the
canonical radiation-formula family — in natural SI units with explicit physical constants.

> **Python binding status:** `green_fn_static`, `larmor_formula`, `dipole_radiation`,
> `dipole_angular_distribution`, `field_strength_point_charge`, and
> `compton_wavelength_shift` are bound. The retarded/advanced Green's functions and the
> remaining multipole helpers are Rust-only. See the [API page](../api/em.md).

**Every figure on this page was produced by calling these functions** — see
`docs/make_em_figures.py`. Closed forms, where shown, are dashed reference lines and
are labelled as such; the fitted log–log slopes are measured from the library's own
output.

---

## 1. Maxwell's Equations

```{figure} ../_static/figures/em_point_charge.svg
:alt: Coulomb field of a point charge, linear and log-log
:width: 100%

The canonical solution of Gauss's law, from `field_strength_point_charge`. The fitted
log–log slope is $-2.0000$, recovering $E_r\propto r^{-2}$ from the returned values
rather than asserting it.
```

### 1.1 Differential Form

$$
\nabla\cdot\mathbf{E} = \frac{\rho}{\varepsilon_0}, \qquad
\nabla\times\mathbf{E} = -\frac{\partial\mathbf{B}}{\partial t}
$$

$$
\nabla\cdot\mathbf{B} = 0, \qquad
\nabla\times\mathbf{B} = \mu_0\mathbf{J} + \mu_0\varepsilon_0\frac{\partial\mathbf{E}}{\partial t}
$$

### 1.2 Covariant Form (4D)

The electromagnetic field tensor $F^{\mu\nu}$ encodes both $\mathbf{E}$ and $\mathbf{B}$:

$$
F^{\mu\nu} = \begin{pmatrix}0&-E_x/c&-E_y/c&-E_z/c\\
E_x/c&0&-B_z&B_y\\
E_y/c&B_z&0&-B_x\\
E_z/c&-B_y&B_x&0\end{pmatrix}
$$

Maxwell's equations become

$$
\partial_\mu F^{\mu\nu} = \mu_0 J^\nu, \qquad
\partial_{[\lambda}F_{\mu\nu]} = 0
$$

The second equation is the **Bianchi identity** — it is automatic from $F_{\mu\nu} =
\partial_\mu A_\nu - \partial_\nu A_\mu$ (see [Gauge](gauge.md) for the non-abelian
generalization).

### 1.3 Gauge Invariance

Under $A_\mu \to A_\mu + \partial_\mu\Lambda$ (abelian gauge transformation):

$$
F_{\mu\nu} \to F_{\mu\nu} + \partial_\mu\partial_\nu\Lambda - \partial_\nu\partial_\mu\Lambda = F_{\mu\nu}
$$

since partial derivatives commute. The physics is invariant — only gauge-variant quantities
($A_\mu$ itself) change.

---

## 2. Green's Functions

```{figure} ../_static/figures/em_green_static.svg
:alt: Static Green's function and its inverse-r slope
:width: 100%

`green_fn_static` in the convention $\nabla^2G=-\delta^3(\mathbf r)$, so
$G=+1/(4\pi r)$ and a positive charge yields a positive potential — the same sign
convention as the retarded Green's function below. The measured log–log slope is
$-1.0000$.
```

### 2.1 Helmholtz/Poisson Operator

The Green function $G(\mathbf{r},\mathbf{r}')$ solves

$$
\nabla^2 G = -\delta(\mathbf{r}-\mathbf{r}')
$$

and builds solutions by convolution: $\phi(\mathbf{r}) = \int G(\mathbf{r},\mathbf{r}')
\rho(\mathbf{r}')/\varepsilon_0\,d^3r'$.

### 2.2 Static Green's Function (3D Poisson)

$$
\boxed{G_{\rm static}(\mathbf{r}) = \frac{1}{4\pi|\mathbf{r}-\mathbf{r}'|}
}
$$

*Proof.* In spherical coordinates centered at $\mathbf{r}'$, $G = G(r)$, and away from
the origin $\nabla^2 G = \frac{1}{r^2}\partial_r(r^2\partial_r G) = 0$, so
$r^2\partial_r G$ is a constant and $G = A/r + C$. Integrating $\nabla^2 G =
-\delta^3(\mathbf{r})$ over a ball of radius $r$ and applying the divergence theorem
gives $\oint\nabla G\cdot d\mathbf{S} = -1$; with $\nabla G = -(A/r^2)\hat{\mathbf r}$
this is $-4\pi A = -1$, hence $A = 1/(4\pi)$. Requiring $G\to0$ as $r\to\infty$ sets
$C = 0$. $\square$

### 2.3 Retarded Green's Function (Causal)

$$
\boxed{G_{\rm ret}(\mathbf{r},t;\,\mathbf{r}',t') = \frac{\delta\!\bigl(t'-t+|\mathbf{r}-\mathbf{r}'|/c\bigr)}{4\pi|\mathbf{r}-\mathbf{r}'|}
}
$$

This is the **causal** propagator: $G_{\rm ret}=0$ for $t'<t$ (no effect before the
cause). The delta function enforces that the signal propagates at speed $c$.

### 2.4 Advanced Green's Function

$$
G_{\rm adv}(\mathbf{r},t;\,\mathbf{r}',t') = \frac{\delta\!\bigl(t'-t-|\mathbf{r}-\mathbf{r}'|/c\bigr)}{4\pi|\mathbf{r}-\mathbf{r}'|}
$$

The **anti-causal** propagator: vanishes for $t'>t$. The combination $G_{\rm ret} -
G_{\rm adv}$ is proportional to the commutator of the electromagnetic field operators in
quantum field theory.

### 2.5 Hertz Potential

The retarded potentials in the Lorenz gauge ($\nabla\cdot\mathbf{A}+\mu_0\varepsilon_0
\partial\phi/\partial t = 0$) are

$$
\phi(\mathbf{r},t) = \frac{1}{4\pi\varepsilon_0}\int\frac{\rho(\mathbf{r}',t_{\rm ret})}{|\mathbf{r}-\mathbf{r}'|}\,d^3r'
$$

$$
\mathbf{A}(\mathbf{r},t) = \frac{\mu_0}{4\pi}\int\frac{\mathbf{J}(\mathbf{r}',t_{\rm ret})}{|\mathbf{r}-\mathbf{r}'|}\,d^3r'
$$

where $t_{\rm ret} = t - |\mathbf{r}-\mathbf{r}'|/c$ is the retarded time.

---

## 3. Radiation

### 3.1 Larmor Formula

```{figure} ../_static/figures/em_larmor.svg
:alt: Larmor power scaling quadratically in acceleration and in charge
:width: 100%

`larmor_formula` is quadratic in both arguments: the fitted log–log slopes are $2.000$
against acceleration and $2.000$ against charge. The dotted marker locates the
elementary charge.
```

**Theorem (Larmor).** A point charge $q$ with acceleration $a$ radiates total power

$$
\boxed{P = \frac{2}{3}\frac{q^2 a^2}{4\pi\varepsilon_0 c^3}}
$$

*Proof.* The angular distribution of radiated power is obtained from the Poynting vector
integrated over a sphere. In the far zone, the辐射 electric field is $\mathbf{E}_{\rm rad}
= \frac{q}{4\pi\varepsilon_0 c^2}\frac{\hat{n}\times(\hat{n}\times\mathbf{a})}{r}$, and
the Poynting flux $\mathbf{S} = \frac{1}{\mu_0}\mathbf{E}\times\mathbf{B}$. Integrating
$\mathbf{S}\cdot\hat{n}$ over the sphere and using the dipole pattern $\langle\sin^2\theta
\rangle = 2/3$ gives the Larmor formula. $\square$

### 3.2 Dipole Radiation

```{figure} ../_static/figures/em_dipole.svg
:alt: Dipole power against frequency and the sin-squared angular pattern
:width: 100%

Left: `dipole_radiation` gives a fitted slope of $4.000$, the $\omega^4$ law that makes
the sky blue. Right: `dipole_angular_distribution` traces the $\sin^2\theta$ torus with
its null along the dipole axis. The two are mutually consistent — integrating the
angular pattern over the sphere reproduces the total power to a relative
$5\times10^{-15}$, since $\int\sin^2\theta\,d\Omega=8\pi/3$.
```

For an oscillating electric dipole $\mathbf{p}(t) = p_0\cos(\omega t)\,\hat{z}$:

$$
\boxed{P = \frac{\mu_0\omega^4 p_0^2}{12\pi c}
      = \frac{\omega^4 p_0^2}{12\pi\varepsilon_0 c^3}}
$$

Both forms are the same number, since $\mu_0=1/(\varepsilon_0c^2)$. Note the **single**
power of $c$ alongside $\mu_0$: an electric dipole moment carries units of
$\mathrm{C\,m}$, and $[\mu_0p_0^2\omega^4/c]=\mathrm{W}$. The magnetic-dipole result in
§3.3 keeps $c^3$ because a magnetic moment has different dimensions — the two must not be
conflated.

The $\omega^4$ dependence is the **Rayleigh scattering law** — small particles scatter
shorter wavelengths much more efficiently (explaining the blue sky).

**Angular distribution:**

$$
\frac{dP}{d\Omega} = \frac{\mu_0\omega^4 p_0^2}{32\pi^2 c}\sin^2\theta
$$

This integrates to the boxed total power by construction, because
$\int\sin^2\theta\,d\Omega = 8\pi/3$ and $\tfrac{8\pi}{3}\cdot\tfrac1{32\pi^2}
= \tfrac1{12\pi}$.

The $\sin^2\theta$ pattern is the classic **dipole lobe** — maximum radiation perpendicular
to the dipole axis, zero along it.

### 3.3 Magnetic Dipole

For an oscillating magnetic dipole $m(t) = m_0\cos(\omega t)$:

$$
P = \frac{\mu_0 m_0^2\omega^4}{12\pi c^3}
$$

Same $\omega^4$ scaling but with magnetic moment replacing electric dipole moment.

### 3.4 Thomson Scattering

A free electron oscillating in an incident wave radiates with cross-section

$$
\sigma_T = \frac{8\pi}{3}r_e^2, \qquad r_e = \frac{e^2}{4\pi\varepsilon_0 m_e c^2} \approx 2.82\times 10^{-15}\,\text{m}
$$

where $r_e$ is the classical electron radius. This is frequency-independent (Rayleigh
regime for $\lambda \gg r_e$).

### 3.5 Compton Scattering

```{figure} ../_static/figures/em_compton.svg
:alt: Compton wavelength shift against scattering angle
:width: 70%
:align: center

`compton_wavelength_shift` normalised by the Compton wavelength
$\lambda_C=h/m_ec=2.4263$ pm: no shift in the forward direction, and exactly
$2\lambda_C$ on backscatter.
```

The wavelength shift in Compton scattering (photon off free electron) is

$$
\Delta\lambda = \lambda' - \lambda = \frac{h}{m_e c}(1-\cos\theta) = \lambda_C(1-\cos\theta)
$$

where $\lambda_C = h/(m_e c) \approx 2.43\times 10^{-12}$ m is the Compton wavelength.
This is the `compton_wavelength_shift` function.

---

## 4. Routines

| Routine | Formula |
|---------|---------|
| `green_fn_static(r)` | $1/(4\pi r)$ |
| `green_fn_retarded(r, t)` | causal $\delta$-propagator |
| `green_fn_advanced(r, t)` | anti-causal |
| `larmor_formula(q, a, eps0, c)` | $\frac{2}{3}\frac{q^2a^2}{4\pi\varepsilon_0c^3}$ |
| `dipole_radiation(p0, omega, mu0, c)` | $\frac{\mu_0\omega^4p_0^2}{12\pi c}$ |
| `radiation_power_oscillating(p0, omega, ...)` | dipole with cosine time dependence |
| `dipole_angular_distribution(theta, p0, ...)` | $dP/d\Omega$ |
| `magnetic_dipole_radiation(m, omega, ...)` | $\frac{\mu_0 m^2\omega^4}{12\pi c^3}$ |
| `field_strength_point_charge(...)` | $\mathbf{E}$, $\mathbf{B}$ of a point charge |
| `compton_wavelength_shift(...)` | $\Delta\lambda = \lambda_C(1-\cos\theta)$ |

---

## 5. Usage Examples

### Larmor radiation of an accelerating electron

```python
from pysicrs import larmor_formula

e = 1.602176634e-19     # C
eps0 = 8.8541878128e-12 # F/m
c = 2.99792458e8        # m/s

g_accel = 9.81
a = g_accel / 2
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

for th in (0.0, 0.7854, 1.5708):  # 0°, 45°, 90°
    dP = dipole_angular_distribution(th, 1e-29, 2e9, mu0=1.25663706212e-6, c=2.99792458e8)
    print(f"θ={th:.2f}: dP/dΩ = {dP:.3e}")
# peaks at θ=90°, vanishes at θ=0° (along the dipole axis)
```

---

## 6. Advantages & Limitations

✅ Exact constants in SI units — no ambiguity between unit systems

✅ Retarded/advanced propagation naturally causal

✅ Cover the textbook radiation $\omega^4$ laws used in antenna and scattering theory

❌ Point charges only (finite-size and multipole corrections modeled by extra terms)

❌ No numerical solving of full wave/curl equations — that's the PDE Maxwell FDTD

❌ No dielectric/magnetization response functions; those live in CASIMIR's Lifshitz machinery

---

## 7. References

1. Jackson, J.D. (1999). *Classical Electrodynamics*, 3rd ed. Wiley.
2. Griffiths, D.J. (2017). *Introduction to Electrodynamics*, 4th ed. Cambridge.
3. Landau, L.D. & Lifshitz, E.M. (1975). *The Classical Theory of Fields*, 4th ed.
4. Zangwill, A. (2012). *Modern Electrodynamics*. Cambridge.

---

## 8. Related Topics

- [PDE](pde.md) — Maxwell 3D FDTD complements the analytic formulas here
- [Casimir](casimir.md) — vacuum energy is a mode sum of EM modes between plates
- [Linear Algebra](linalg.md) — field tensors $F_{\mu\nu}$ interact with raise/lower
- [Gauge](gauge.md) — electrodynamics is the abelian U(1) limit of Yang–Mills
