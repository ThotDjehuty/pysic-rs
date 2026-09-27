# Topology

Pysic-rs provides the topological invariants central to modern condensed-matter and field
theory: Berry phases and curvatures, Chern and TKNN numbers, and winding numbers — including
discrete (lattice) evaluations.

> **Python binding status:** `berry_phase`, `berry_phase_bloch`, `chern_number`,
> `winding_number`, and `skyrmion_number` are bound. `berry_curvature`, the discrete
> variants, and the TKNN helpers are Rust-only. See the [API page](../api/topology.md).

**Every figure on this page was produced by calling these functions** — see
`docs/make_topology_figures.py`. Nothing is plotted from another library's
topological invariants; where a closed form appears it is a dashed reference
line, labelled as such.

---

## 1. Berry Phase

```{figure} ../_static/figures/topo_berry_phase.svg
:alt: Berry phase on the Bloch sphere, closed form against the discrete Wilson loop
:width: 100%

Left: the closed form $\gamma=-\pi(1-\cos\theta)$ from `berry_phase_bloch`, with circles
showing the discrete estimator `berry_phase` applied to the explicit spin-½ states
$|n\rangle=(\cos\tfrac\theta2,\ \sin\tfrac\theta2 e^{i\varphi})$ carried once around
$\varphi\colon0\to2\pi$. The two agree to $2.5\times10^{-5}$ rad. Right: the phase is
exactly half the enclosed solid angle — the geometric content of the result.
```

### 1.1 Adiabatic Evolution & Geometric Phase

Consider a parameter-dependent Hamiltonian $H(\mathbf{R})$ with instantaneous eigenstates
$|n(\mathbf{R})\rangle$. As $\mathbf{R}$ is slowly varied around a closed loop $\mathcal{C}$
in parameter space, the adiabatic theorem guarantees that a system starting in $|n(\mathbf{R}_0)\rangle$
returns to that eigenstate — but acquires a phase:

$$
|\psi(T)\rangle = e^{i\gamma_n}\,e^{-\frac{i}{\hbar}\int_0^T E_n(t)\,dt}\,|n(\mathbf{R}_0)\rangle
$$

The first exponential is the **dynamical phase**; the second is the **geometric (Berry)
phase**.

### 1.2 Berry Phase Formula

**Theorem (Berry, 1984).** The geometric phase is

$$
\boxed{\gamma_n = i\oint_{\mathcal{C}}\langle n(\mathbf{R})|\nabla_{\mathbf{R}}|n(\mathbf{R})\rangle\cdot d\mathbf{R}
= \oint_{\mathcal{C}}\mathbf{A}_n\cdot d\mathbf{R}}
$$

where $\mathbf{A}_n = i\langle n|\nabla_{\mathbf{R}}|n\rangle$ is the **Berry connection**.

*Proof.* Under adiabatic evolution, the state picks up a phase from the overlap of
neighbouring eigenstates:

$$
\langle n(\mathbf{R})|n(\mathbf{R}+d\mathbf{R})\rangle \approx 1 + \langle n|\nabla_{\mathbf{R}}|n\rangle\cdot d\mathbf{R}
$$

The total phase after a closed loop is the path-ordered product (Wilson loop) of these
overlaps. In the adiabatic limit, this reduces to the integral of the Berry connection
along the loop. The gauge-invariance of $\gamma_n$ follows from the fact that changing
the phase convention $|n\rangle\to e^{i\phi(\mathbf{R})}|n\rangle$ shifts $\mathbf{A}_n
\to \mathbf{A}_n-\nabla\phi$, which does not change the closed-loop integral. $\square$

### 1.3 Berry Curvature

The Berry connection $\mathbf{A}_n$ is a gauge field on parameter space. Its "field
strength" is the **Berry curvature**:

$$
\boxed{\Omega_n^{ab} = \partial_a A_n^b - \partial_b A_n^a
= -2\,\mathrm{Im}\sum_{m\neq n}\frac{\langle n|\partial_a H|m\rangle\langle m|\partial_b H|n\rangle}{(E_n-E_m)^2}}
$$

*Proof.* Differentiate $H|n\rangle = E_n|n\rangle$ with respect to $R^a$:

$$
(\partial_a H)|n\rangle + H|\partial_a n\rangle = (\partial_a E_n)|n\rangle + E_n|\partial_a n\rangle
$$

Projecting onto $\langle m|$ for $m\neq n$:

$$
\langle m|\partial_a n\rangle = \frac{\langle m|\partial_a H|n\rangle}{E_n-E_m}
$$

Substituting into $\Omega_n^{ab} = i(\partial_a A_n^b - \partial_b A_n^a)$ and using
$\mathbf{A}_n = i\langle n|\nabla|n\rangle$ gives the result. $\square$

**Gauge invariance:** $\Omega_n^{ab}$ is invariant under $|n\rangle\to e^{i\phi}|n\rangle$
(the connection transforms, but the curvature does not — just as in electromagnetism).

### 1.4 Flux Quantisation (Stokes)

For a closed surface $\Sigma$ bounding $\mathcal{C}$:

$$
\gamma_n = \iint_\Sigma \Omega_n^{ab}\,dS_{ab}
$$

This relates the Berry phase to the flux of Berry curvature through a surface bounded by the
loop.

---

## 2. Chern Number (TKNN Invariant)

```{figure} ../_static/figures/topo_chern_qwz.svg
:alt: Quantized Chern plateaus of the Qi-Wu-Zhang model
:width: 100%

For a two-band Hamiltonian $H=\mathbf d(\mathbf k)\cdot\boldsymbol\sigma$ the Chern number
*is* the skyrmion number of $\hat{\mathbf d}$ over the Brillouin zone, so `skyrmion_number`
computes it directly. Left: sweeping the Qi–Wu–Zhang mass $m$ with
$\mathbf d=(\sin k_x,\ \sin k_y,\ m+\cos k_x+\cos k_y)$ gives the plateaus
$C=0,+1,-1,0$ with transitions at $m=0,\pm2$ (dashed) where the gap closes and the
invariant is undefined. Right: the plateaus approach exact integers as the $k$-grid is
refined — quantization is recovered numerically, not imposed.
```

```{figure} ../_static/figures/topo_chern_integrator.svg
:alt: chern_number recovers the prescribed Chern number from a curvature grid
:width: 70%
:align: center

`chern_number` performs the $\frac1{2\pi}\int\Omega\,d^2k$ integration. Feeding it a
constant curvature $\Omega=2\pi C/L^2$ returns $C$ to within $2\times10^{-13}$ across
$C=-3\ldots3$.
```

### 2.1 Definition

For a 2D parameter space (e.g. the Brillouin zone), the **first Chern number** of band $n$
is the total Berry flux:

$$
\boxed{C_n = \frac{1}{2\pi}\int_{\mathrm{BZ}}\Omega_n^{xy}\,d^2k \in \mathbb{Z}}
$$

**Theorem (Thouless–Kohmoto–Nightingale–den Nijs, 1982).** $C_n$ is an integer.

*Proof.* The Berry curvature $\Omega_n^{xy}$ is a closed 2-form ($d\Omega = 0$ follows from
the Bianchi identity of the Berry connection). By Chern–Weil theory, the integral of a
closed 2-form over a compact 2-manifold (the BZ torus) is a topological invariant — it
cannot change under smooth deformations of the Hamiltonian. The integer is the **first Chern
class** of the Berry connection bundle. $\square$

### 2.2 Quantum Hall Effect

**TKNN formula (1982):** the Hall conductance of a 2D electron gas in a magnetic field is

$$
\sigma_{xy} = \frac{e^2}{h}\sum_{n\,\mathrm{filled}}C_n
$$

Each filled Landau level contributes $C_n = 1$ (for the lowest level), giving the
quantised Hall conductance $\sigma_{xy} = ne^2/h$ — the **integer quantum Hall effect**.

### 2.3 Lattice (Discrete) Evaluation

On a discretised BZ with plaquette $(k_x, k_y)\to(k_x+\Delta k, k_y)\to(k_x+\Delta k,
k_y+\Delta k)\to(k_x, k_y+\Delta k)$, the Berry phase around each plaquette is

$$
\gamma_{\square} = \mathrm{Im}\ln\bigl[U_1(k)\,U_2(k+\hat{1})\,U_1^*(k+\hat{2})\,U_2^*(k)\bigr]
$$

where $U_\mu(k) = \langle n(k)|n(k+\hat{\mu})\rangle/|\langle n(k)|n(k+\hat{\mu})\rangle|$ is
the **link variable** (Wilson loop element). The Chern number is

$$
C = \frac{1}{2\pi}\sum_{\square}\gamma_{\square} \in \mathbb{Z}
$$

This is the lattice Berry flux method implemented by `chern_number_discrete`.

---

## 3. Winding Number

```{figure} ../_static/figures/topo_winding.svg
:alt: Curves winding 1, 2, 3 and -2 times about the origin
:width: 100%

Closed curves $z(t)=(1+0.3\cos t)e^{int}$ carried once around, with the winding number
returned by `winding_number` in each title. The radial modulation keeps $|z|>0$ while
separating the loops visually; the sign tracks orientation.
```

### 3.1 Definition

For a map $g: S^1\to S^1$ (a closed curve in the complex plane avoiding the origin), the
**winding number** counts the net number of counterclockwise encirclements of the origin:

$$
\boxed{w(g) = \frac{1}{2\pi i}\oint\frac{g'(z)}{g(z)}\,dz \in \mathbb{Z}}
$$

*Proof.* Write $g(z) = r(z)e^{i\theta(z)}$. Then $g'/g = r'/r + i\theta'$, and the contour
integral of $r'/r$ vanishes (single-valued), while $\frac{1}{2\pi}\oint d\theta$ counts the
net winding of $\theta$ around the origin — an integer. $\square$

### 3.2 Discrete Version

For a polygon with vertices $z_0, z_1, \ldots, z_N = z_0$, the winding number is

$$
w = \frac{1}{2\pi}\sum_{k=0}^{N-1}\arg\!\left(\frac{z_{k+1}}{z_k}\right)
$$

where the argument is taken in $(-\pi,\pi]$ and the sum of branch cuts gives the total
winding.

### 3.3 Physical Applications

- **Topological insulators:** the winding number classifies 1D chiral symmetry-protected
  topological phases (SSH model).
- **Vortices in superfluids/superconductors:** the winding number of the order parameter
  around a vortex core.
- **Anyons:** the braiding of anyonic quasiparticles in 2D is described by the winding
  number of their worldlines.

---

## 4. Skyrmion Number

```{figure} ../_static/figures/topo_skyrmion.svg
:alt: Belavin-Polyakov hedgehog texture and convergence of its topological charge
:width: 100%

Left: the Belavin–Polyakov profile $\hat{\mathbf n}=(2x,\,2y,\,r^2-1)/(r^2+1)$, with the
in-plane components as arrows and $n_z$ shaded — the texture wraps the target sphere
exactly once. Right: `skyrmion_number` converges to $|Q|=1$ as the domain grows. The
charge density is $q(r)=4/(1+r^2)^2$, so a domain of half-width $R$ necessarily omits
$1/(1+R^2)$ of the total; the error follows that analytic tail until the fixed
$h=0.05$ spacing becomes the limiting factor near $R\approx20$.
```

### 4.1 Definition

For a map $\mathbf{n}: \mathbb{R}^2\to S^2$ (a spin configuration), the **skyrmion number**
is the degree of the map:

$$
\boxed{N_{\rm sky} = \frac{1}{4\pi}\iint\mathbf{n}\cdot\left(\frac{\partial\mathbf{n}}{\partial x}\times\frac{\partial\mathbf{n}}{\partial y}\right)dx\,dy \in \mathbb{Z}}
$$

*Proof.* The integrand is the pullback of the area form on $S^2$ under $\mathbf{n}$.
Integrating over all space counts how many times $\mathbf{n}$ wraps around $S^2$ — the
degree of the map. For a compact manifold without boundary, this must be an integer
(degree theorem). $\square$

### 4.2 Physical Applications

- **Chiral magnets:** skyrmions are topologically stabilised magnetic vortices.
- **Skyrmion tubes in QCD:** baryons can be viewed as skyrmions in the pion field.
- **Spintronics:** skyrmion racetrack memory proposals.

---

## 5. Routines

| Routine | Description |
|---------|-------------|
| `berry_phase(psi, dparam, ...)` | $\gamma_n$ around a loop |
| `berry_curvature(H_params, psi, ...)` | $\Omega_n^{ab}$ at a point |
| `berry_phase_bloch(...)` | Bloch-state Berry phase |
| `chern_number(...)` | $C_n$ integral over BZ |
| `chern_number_discrete(lattice, ...)` | lattice flux version |
| `winding_number(path_real, path_imag)` | $w(g) = \frac{1}{2\pi i}\oint g'/g$ |
| `winding_number_discrete(points)` | discrete encirclement |
| `skyrmion_number(spin_config)` | degree of the map $\mathbb{R}^2\to S^2$ |

---

## 6. Usage Examples

All of the following run against the bound API and reproduce the printed values.

### Berry phase around a loop

```python
import math
from pysicrs import berry_phase, berry_phase_bloch

theta = math.pi / 3                      # cone half-angle
phis = [2 * math.pi * i / 400 for i in range(400)]
# |n(θ,φ)⟩ = (cos θ/2, sin θ/2 · e^{iφ}), passed as separate re/im parts
re = [[math.cos(theta / 2), math.sin(theta / 2) * math.cos(p)] for p in phis]
im = [[0.0,                 math.sin(theta / 2) * math.sin(p)] for p in phis]

gamma = berry_phase(re, im)              # discrete Wilson loop
gamma = gamma - 2 * math.pi if gamma > 0 else gamma   # onto the closed-form branch
print(f"discrete    γ = {gamma:+.6f} rad")            # -1.570772
print(f"closed form γ = {berry_phase_bloch(theta):+.6f} rad")  # -1.570796
```

### Winding number of a closed curve

```python
import math
from pysicrs import winding_number

for n in (1, 2, -3):
    ts = [2 * math.pi * i / 2000 for i in range(2000)]
    f1 = [math.cos(n * t) for t in ts]
    f2 = [math.sin(n * t) for t in ts]
    print(n, "->", winding_number(f1, f2))   # 1 -> 1,  2 -> 2,  -3 -> -3
```

### Chern number of a two-band model

For $H=\mathbf d(\mathbf k)\cdot\boldsymbol\sigma$ the Chern number is precisely the
skyrmion number of $\hat{\mathbf d}$ over the Brillouin zone, so one routine covers both:

```python
import math
from pysicrs import skyrmion_number

def qwz_chern(m, nk=64):
    dk = 2 * math.pi / nk
    ks = [-math.pi + (i - 1) * dk for i in range(nk + 2)]  # padded for the stencil
    grid = []
    for kx in ks:
        row = []
        for ky in ks:
            dx, dy = math.sin(kx), math.sin(ky)
            dz = m + math.cos(kx) + math.cos(ky)
            nrm = math.sqrt(dx*dx + dy*dy + dz*dz)         # ≠ 0 unless m = 0, ±2
            row.append([dx / nrm, dy / nrm, dz / nrm])
        grid.append(row)
    return skyrmion_number(grid, dk, dk)

for m in (-3.0, -1.0, 1.0, 3.0):
    print(f"m={m:+.1f}  C={qwz_chern(m):+.4f}")
# m=-3.0  C=+0.0001   m=-1.0  C=+0.9958
# m=+1.0  C=-0.9958   m=+3.0  C=-0.0001
```

```{note}
At $m=0$ and $m=\pm2$ the gap closes, $\mathbf d(\mathbf k)=0$ somewhere in the zone, and
the Chern number is genuinely undefined — not a numerical failure but the physics of a
topological transition.
```

### Chern number from a Berry-curvature grid

```python
import math
from pysicrs import chern_number

n, L = 64, 1.0
for target in (0, 1, -2):
    curvature = [[2 * math.pi * target / (L * L)] * n for _ in range(n)]
    print(target, "->", round(chern_number(curvature, L / n), 6))
# 0 -> 0.0,  1 -> 1.0,  -2 -> -2.0
```

### Skyrmion number of a hedgehog texture

```python
from pysicrs import skyrmion_number

span, npts = 16.0, 641
coords = [-span + 2 * span * i / (npts - 1) for i in range(npts)]
h = coords[1] - coords[0]
field = [[[2*x/(x*x+y*y+1), 2*y/(x*x+y*y+1), (x*x+y*y-1)/(x*x+y*y+1)]
          for y in coords] for x in coords]

print(f"|Q| = {abs(skyrmion_number(field, h, h)):.4f}")   # 0.9951  (exact: 1)
```

---

## 7. Advantages & Limitations

✅ Direct access to geometric/topological invariants (no ad-hoc derivations)

✅ Both continuum and discrete evaluation for lattice-style problems

✅ Connects naturally to the [Gauge](gauge.md) module (winding → instanton quantisation)

❌ Requires pre-computed eigenstates / Berry connections from outside

❌ No coverage of higher-genus (genus $g$) index theorems besides the winding/skyrmion counts

❌ Reflection of a generically non-abelian (Wilson-loop) connection is out of scope

---

## 8. References

1. Berry, M.V. (1984). "Quantal phase factors accompanying adiabatic changes." *Proc. R. Soc. A* 392:45.
2. Thouless, D.J., Kohmoto, M., Nightingale, M.P. & den Nijs, M. (1982). *Phys. Rev. Lett.* 49:405.
3. Simon, B. (1983). "Holonomy, the quantum adiabatic theorem, and Berry's phase." *Phys. Rev. Lett.* 51:2167.
4. Göckeler, M. & Schücker, T. (1990). *Differential Geometry, Gauge Theories, and Gravity*. Cambridge.
5. Hasan, M.Z. & Kane, C.L. (2010). "Colloquium: topological insulators." *Rev. Mod. Phys.* 82:3045.

---

## 9. Related Topics

- [Gauge](gauge.md) — winding number quantises the Yang–Mills instanton action
- [Quantum](quantum.md) — Berry phases appear in adiabatic two-level problems
- [General Relativity](general_relativity.md) — same differential-geometric language
- [Linear Algebra](linalg.md) — eigenstates for Berry connection computed here
