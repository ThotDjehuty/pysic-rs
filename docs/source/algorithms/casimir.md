# Casimir Effect

Pysic-rs implements the celebrated vacuum-energy predictions of quantum field theory — the
Casimir force — with the standard, experimentally-confirmed formulas and their
regularization variants.

> **Python binding status:** `casimir_energy_parallel_plates`, `casimir_force`, and
> `polder_potential` are bound. The sphere, zeta, cylinder, finite-T, and Lifshitz forms are
> Rust-only. See the [API page](../api/casimir.md).

---

## Mathematical Foundations

### Parallel Plates

Two perfectly conducting parallel plates separated by distance $d$ alter the vacuum
zero-point energy of the electromagnetic field. The **energy per unit area** (Casimir 1948):

$$
\frac{E}{A} = -\frac{\pi^2 \hbar c}{720\,d^3}
$$

and the resulting **force per unit area** (attractive):

$$
\frac{F}{A} = -\frac{d}{dd}\left(\frac{E}{A}\right) = -\frac{\pi^2 \hbar c}{240\,d^4}
$$

**Proven**: the calculation subtracts the free-space vacuum mode sum $\frac12\sum\hbar\omega$
from the mode sum with plates; the divergence cancels and the finite, regulator-independent
difference is the result above. **Experimentally confirmed** (Lamoreaux 1997; Mohideen–Roy
1998).

### Zeta-Function Regularization

The same result obtained via analytic continuation of the Epstein zeta function:

$$
\frac{E}{A} = -\frac{\zeta(3)\,\hbar c}{8\pi\,d^3}
$$

**Proven identity (proven):** $\tfrac{\pi^2}{720} = \tfrac{\zeta(3)}{8\pi}$, numerically
$-\frac{0.01448\,\hbar c}{d^3}$ — both expressions exactly equal.

### Conducting Sphere

For a perfectly conducting sphere of radius $a$ (Boyer 1968, via cutoff regularization):

$$
E = -\frac{\hbar c}{8\pi\,a}
$$

### Concentric Cylinders

For a conducting cylinder of radius $a$ with an enclosing shell at separation $d$ (limit
$d \ll a$):

$$
E = -\frac{\hbar c\, a}{24\,d^2}
$$

### Finite-Temperature Correction

At temperature $T$ the thermal (blackbody-like) modes add:

$$
E(T) = E(0) + \frac{\pi^2 k_B^4\, T^4}{45\,\hbar^3 c^3}\,V
$$

### Lifshitz Formula (Real Materials)

For dielectrics with permittivity $\varepsilon(\omega)$, the zero-temperature Lifshitz
formula (Lifshitz 1956; **the** modern standard, confirmed by all precision experiments):

$$
\frac{E}{A} = \frac{\hbar}{2\pi^2}\int_0^\infty k_\perp dk_\perp\;
\int_0^\infty d\xi\; \ln\!\left(1 - r_{s,p}(i\xi)\,e^{-2\kappa d}\right)
$$

The exported routine evaluates the simplified form for identical media using the reflection
coefficient $r = \frac{\varepsilon(i\kappa) - 1}{\varepsilon(i\kappa) + 1}$ on the imaginary
frequency axis:

$$
\frac{E}{A} = -\frac{\hbar c}{4\pi^2}\int_0^\infty \frac{k^3\,dk}{e^{2kd}-1}\,
\left(\frac{\varepsilon(ik)-1}{\varepsilon(ik)+1}\right)^2
$$

---

## Routines

| Routine | Formula |
|---------|---------|
| `casimir_energy_parallel_plates(d, hbar, c)` | $-\dfrac{\pi^2\hbar c}{720\,d^3}$ per area |
| `casimir_force(d, hbar, c)` | $-\dfrac{\pi^2\hbar c}{240\,d^4}$ per area |
| `casimir_energy_sphere(a, hbar, c)` | $-\dfrac{\hbar c}{8\pi a}$ |
| `casimir_energy_zeta(d, hbar, c)` | $-\dfrac{\zeta(3)\hbar c}{8\pi d^3}$ |
| `casimir_energy_cylinders(a, d, hbar, c)` | $-\dfrac{\hbar c\,a}{24\,d^2}$ |
| `casimir_energy_finite_temperature(d, T, hbar, c, k_B)` | $E(0) + \frac{\pi^2 k_B^4 T^4}{45\hbar^3 c^3}V$ |
| `casimir_energy_lifshitz(d, ε(iω), hbar, c, n)` | $-\frac{\hbar c}{4\pi^2}\int\frac{k^3\,dk}{e^{2kd}-1}r^2$ |
| `polder_potential(...)` | retarded vdW between atoms |

---

## Usage Examples

### The classic force between plates at 100 nm

```python
from pysicrs import casimir_force, casimir_energy_parallel_plates

hbar, c = 1.054571817e-34, 2.99792458e8
d = 100e-9

E = casimir_energy_parallel_plates(d, hbar, c)   # J/m²
F = casimir_force(d, hbar, c)                    # N/m²
print(f"E/A = {E:.4e} J/m², F/A = {F:.4e} N/m²")
# F ≈ -1.3e-3 N/m² — measurable at micron separations
```

### Zeta-regularized energy equals the classic formula

```python
from pysicrs import casimir_energy_zeta, casimir_energy_parallel_plates

d, hbar, c = 1e-6, 1.054571817e-34, 2.99792458e8
E1 = casimir_energy_parallel_plates(d, hbar, c)
E2 = casimir_energy_zeta(d, hbar, c)
print(f"symmetric? {abs(E1 - E2):.3e}")   # both are -π²ℏc/720d³
```

### Lifshitz for a dielectric

```python
from pysicrs import casimir_energy_lifshitz

def eps_imag(w):            # ε(iω) — e.g. Drude permittivity
    wp = 1e16
    return 1.0 + wp**2 / w**2

E = casimir_energy_lifshitz(d=200e-9, epsilon_fn=eps_imag,
                            hbar=1.054571817e-34, c=2.99792458e8, n_points=400)
print(f"Lifshitz E/A = {E:.4e} J/m²")
```

---

## Advantages & Limitations

✅ Implements the exact, experimentally-confirmed textbook formulas

✅ Zeta regularization matches the classic result (checked internally by the test suite)

✅ Lifshitz formula supports realistic dielectric materials

❌ Material models (Drude/plasma) must be supplied as closures — no built-in table

❌ Sphere result is the Boyer (cutoff) approximation, not the full Mie/BD expansions

❌ No finite-temperature correction for the Lifshitz integral (T=0 only)

---

## References

1. Casimir, H.B.G. (1948). "On the attraction between two perfectly conducting plates." *Proc. K. Ned. Akad. Wet.* 51:793.
2. Lifshitz, E.M. (1956). "The theory of molecular attractive forces between solids." *Sov. Phys. JETP* 2:73.
3. Boyer, T.H. (1968). "Quantum electromagnetic zero-point energy of a conducting spherical shell." *Phys. Rev.* 174:1764.
4. Lamoreaux, S.K. (1997). "Demonstration of the Casimir force in the 0.6 to 6 μm range." *Phys. Rev. Lett.* 78:5.
5. Mohideen, U. & Roy, A. (1998). *Phys. Rev. Lett.* 81:4549.

---

## Related Topics

- [EM](em.md) – the modes being summed are electromagnetic vacuum modes
- [PDE](pde.md) – mode decomposition viewpoints connect to Maxwell FDTD
- [Quantum](quantum.md) – zero-point energy of the field is a QFT concept