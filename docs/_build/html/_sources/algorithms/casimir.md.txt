# Casimir Effect

Pysic-rs implements the celebrated vacuum-energy predictions of quantum field theory — the
Casimir force — with the standard, experimentally-confirmed formulas and their
regularisation variants.

> **Python binding status:** `casimir_energy_parallel_plates`, `casimir_force`, and
> `polder_potential` are bound. The sphere, zeta, cylinder, finite-T, and Lifshitz forms are
> Rust-only. See the [API page](../api/casimir.md).

---

## 1. Parallel Plates

### 1.1 Setup & Regularisation

Two perfectly conducting parallel plates separated by distance $d$ alter the vacuum
zero-point energy of the electromagnetic field. The naively divergent mode sum
$\frac{1}{2}\sum\hbar\omega$ is regularised by subtracting the free-space contribution.

### 1.2 Casimir's Original Calculation

The electromagnetic modes between the plates have wavevectors $k_n = n\pi/d$ ($n=1,2,
\ldots$) in the direction perpendicular to the plates. The energy per unit area is

$$
\frac{E}{A} = \frac{1}{2}\sum_{n=1}^\infty\int\frac{d^2k_\perp}{(2\pi)^2}\,2\sqrt{k_\perp^2+\left(\frac{n\pi}{d}\right)^2}\,\hbar
$$

(the factor 2 is for polarisations). The sum diverges quartically. Regularise with
$e^{-\alpha\omega}$ and let $\alpha\to 0$ after subtraction of the free-space term:

$$
\boxed{\frac{E}{A} = -\frac{\pi^2\hbar c}{720\,d^3}}
$$

### 1.3 Force

$$
\boxed{\frac{F}{A} = -\frac{d}{dd}\left(\frac{E}{A}\right) = -\frac{\pi^2\hbar c}{240\,d^4}}
$$

The force is **attractive** (the plates are pushed together). Experimentally confirmed:
Lamoreaux (1997) measured $F/A$ to 5% accuracy at separations $0.6$–$6\,\mu\text{m}$;
Mohideen & Roy (1998) confirmed to 1% at $100$–$900\,\text{nm}$.

### 1.4 Physical Interpretation

The Casimir force arises because the plates **restrict** the vacuum fluctuations between
them: only modes with $k_n = n\pi/d$ are allowed. Outside, all modes contribute. The
resulting imbalance in zero-point pressure pushes the plates together.

---

## 2. Zeta-Function Regularisation

### 2.1 Epstein Zeta Function

Define $Z(s) = \sum_{n=1}^\infty(n\pi/d)^{-s} = (d/\pi)^s\zeta(s)$. The regularised
energy is

$$
\frac{E}{A} = -\frac{1}{2}\frac{\sqrt{\pi}}{(4\pi)^{3/2}}\frac{\Gamma(3/2)}{\Gamma(3/2)}\,Z(3)
= -\frac{\zeta(3)\,\hbar c}{8\pi\,d^3}
$$

### 2.2 Equivalence

**Proven identity:** $\frac{\pi^2}{720} = \frac{\zeta(3)}{8\pi}$ (numerically both
$\approx 0.01370\dots$). Both expressions are exactly equal — the zeta-regularised and
the cutoff-regularised results agree.

---

## 3. Conducting Sphere

For a perfectly conducting sphere of radius $a$ (Boyer 1968, via cutoff regularisation):

$$
\boxed{E = -\frac{\hbar c}{8\pi\,a}}
$$

The sign is **attractive** (the sphere is pulled inward), but this is controversial — the
full Mie/BD (Bender–Deutsch) expansion gives corrections that depend on the precise
boundary conditions.

---

## 4. Concentric Cylinders

For a conducting cylinder of radius $a$ with an enclosing shell at separation $d$ (limit
$d\ll a$):

$$
E = -\frac{\hbar c\,a}{24\,d^2}
$$

---

## 5. Finite-Temperature Correction

At temperature $T$, thermal (blackbody-like) modes add to the zero-point energy:

$$
E(T) = E(0) + \frac{\pi^2 k_B^4\,T^4}{45\,\hbar^3 c^3}\,V
$$

The thermal correction dominates at large separations ($d\gtrsim\hbar c/(k_BT)$, the
thermal wavelength). At room temperature, the crossover is at $d\approx 7\,\mu\text{m}$.

---

## 6. Lifshitz Formula (Real Materials)

### 6.1 General Form

For dielectrics with permittivity $\varepsilon(\omega)$, the zero-temperature Lifshitz
formula (Lifshitz 1956) is the **modern standard**, confirmed by all precision experiments:

$$
\frac{E}{A} = \frac{\hbar}{2\pi^2}\int_0^\infty k_\perp\,dk_\perp\int_0^\infty d\xi\;\ln\!\left(1-r_{s,p}(i\xi)\,e^{-2\kappa d}\right)
$$

where $\kappa = \sqrt{k_\perp^2+\xi^2/c^2}$, and $r_{s,p}$ are the Fresnel reflection
coefficients for s- and p-polarisation evaluated on the imaginary frequency axis.

### 6.2 Simplified Form (Identical Media)

For two identical media with reflection coefficient
$r = \frac{\varepsilon(i\kappa)-1}{\varepsilon(i\kappa)+1}$:

$$
\frac{E}{A} = -\frac{\hbar c}{4\pi^2}\int_0^\infty\frac{k^3\,dk}{e^{2kd}-1}\,r^2
$$

### 6.3 Drude Model

For metals, the Drude permittivity on the imaginary axis is $\varepsilon(i\xi) = 1+
\omega_p^2/(\xi(\xi+\gamma))$ where $\omega_p$ is the plasma frequency and $\gamma$ the
damping rate. In the clean limit ($\gamma\to 0$): $\varepsilon(i\xi) = 1+\omega_p^2/\xi^2$.

---

## 7. Polder Potential (Atom–Wall)

The retarded van der Waals potential between an atom and a conducting wall at distance $z$
is (Casimir & Polder 1948):

$$
V(z) = -\frac{3\hbar c\,\alpha_0}{8\pi z^4}\,f(\lambda/z)
$$

where $\alpha_0$ is the static polarisability, $\lambda$ is the atomic transition
wavelength, and $f$ is a crossover function: $f(x)\to 1$ for $x\gg 1$ (retarded regime)
and $f(x)\to\frac{3}{4}x$ for $x\ll 1$ (non-retarded regime).

---

## 8. Routines

| Routine | Formula |
|---------|---------|
| `casimir_energy_parallel_plates(d, hbar, c)` | $-\pi^2\hbar c/(720\,d^3)$ per area |
| `casimir_force(d, hbar, c)` | $-\pi^2\hbar c/(240\,d^4)$ per area |
| `casimir_energy_sphere(a, hbar, c)` | $-\hbar c/(8\pi a)$ |
| `casimir_energy_zeta(d, hbar, c)` | $-\zeta(3)\hbar c/(8\pi d^3)$ |
| `casimir_energy_cylinders(a, d, hbar, c)` | $-\hbar c\,a/(24\,d^2)$ |
| `casimir_energy_finite_temperature(d, T, hbar, c, k_B)` | $E(0) + \pi^2 k_B^4 T^4/(45\hbar^3 c^3)V$ |
| `casimir_energy_lifshitz(d, ε(iω), hbar, c, n)` | $-\hbar c/(4\pi^2)\int k^3 dk\,r^2/(e^{2kd}-1)$ |
| `polder_potential(...)` | retarded vdW between atom and wall |

---

## 9. Usage Examples

### The classic force between plates at 100 nm

```python
from pysicrs import casimir_force, casimir_energy_parallel_plates

hbar, c = 1.054571817e-34, 2.99792458e8
d = 100e-9

E = casimir_energy_parallel_plates(d, hbar, c)
F = casimir_force(d, hbar, c)
print(f"E/A = {E:.4e} J/m², F/A = {F:.4e} N/m²")
```

### Zeta-regularised energy equals the classic formula

```python
from pysicrs import casimir_energy_zeta, casimir_energy_parallel_plates

d, hbar, c = 1e-6, 1.054571817e-34, 2.99792458e8
E1 = casimir_energy_parallel_plates(d, hbar, c)
E2 = casimir_energy_zeta(d, hbar, c)
print(f"symmetric? {abs(E1-E2):.3e}")
```

### Lifshitz for a dielectric

```python
from pysicrs import casimir_energy_lifshitz

def eps_imag(w):
    wp = 1e16
    return 1.0 + wp**2/w**2

E = casimir_energy_lifshitz(d=200e-9, epsilon_fn=eps_imag,
                            hbar=1.054571817e-34, c=2.99792458e8, n_points=400)
print(f"Lifshitz E/A = {E:.4e} J/m²")
```

---

## 10. Advantages & Limitations

✅ Implements the exact, experimentally-confirmed textbook formulas

✅ Zeta regularisation matches the classic result (checked internally)

✅ Lifshitz formula supports realistic dielectric materials

❌ Material models (Drude/plasma) must be supplied as closures — no built-in table

❌ Sphere result is the Boyer (cutoff) approximation, not the full Mie/BD expansions

❌ No finite-temperature correction for the Lifshitz integral ($T=0$ only)

---

## 11. References

1. Casimir, H.B.G. (1948). "On the attraction between two perfectly conducting plates." *Proc. K. Ned. Akad. Wet.* 51:793.
2. Lifshitz, E.M. (1956). "The theory of molecular attractive forces between solids." *Sov. Phys. JETP* 2:73.
3. Boyer, T.H. (1968). "Quantum electromagnetic zero-point energy of a conducting spherical shell." *Phys. Rev.* 174:1764.
4. Lamoreaux, S.K. (1997). "Demonstration of the Casimir force in the 0.6 to 6 μm range." *Phys. Rev. Lett.* 78:5.
5. Mohideen, U. & Roy, A. (1998). *Phys. Rev. Lett.* 81:4549.
6. Bordag, M., Mohideen, U. & Mostepanenko, V. (2001). *New Developments in the Casimir Effect*. Springer.

---

## 12. Related Topics

- [Electromagnetism](em.md) — the modes being summed are electromagnetic vacuum modes
- [PDE](pde.md) — mode decomposition viewpoints connect to Maxwell FDTD
- [Quantum](quantum.md) — zero-point energy of the field is a QFT concept
- [Special Functions](special_functions.md) — $\zeta(3)$ in the regularised result
