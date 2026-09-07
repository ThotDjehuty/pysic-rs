# Topology

Pysic-rs provides the topological invariants central to modern condensed-matter and field
theory: Berry phases and curvatures, Chern and TKNN numbers, and winding numbers — including
a discrete (lattice) evaluation.

> **Python binding status:** `berry_phase`, `berry_phase_bloch`, `chern_number`,
> `winding_number`, and `skyrmion_number` are bound. `berry_curvature`, the discrete variants,
> and the TKNN helpers are Rust-only. See the [API page](../api/topology.md).

---

## Mathematical Foundations

### Berry Phase & Curvature

For a parameter-dependent Hamiltonian $H(\mathbf R)$ with eigenstates $|n(\mathbf R)\rangle$,
the adiabatic evolution acquires a geometric phase (Berry 1984; **proven** — it is gauge
invariant):

$$
\gamma_n = i\oint \langle n(\mathbf R)|\nabla_{\mathbf R} n(\mathbf R)\rangle \cdot d\mathbf R
= \oint \mathbf A_n \cdot d\mathbf R
$$

with **Berry connection** $\mathbf A_n = i\langle n|\nabla_{\mathbf R}|n\rangle$, and the
**Berry curvature** (a gauge-invariant 2-form):

$$
\Omega_n^{ab} = i\left(\partial_{a}A_n^b - \partial_b A_n^a\right)
= -2\operatorname{Im}\sum_{m\neq n}
\frac{\langle n|\partial_a H|m\rangle\langle m|\partial_b H|n\rangle}{(E_n - E_m)^2}
$$

**Proven properties:** $\gamma_n$ is real, gauge-invariant mod $2\pi$, and — crucially — a
**phase** under parameter loops, which explains the Aharonov–Bohm phase and (with spin-orbit
coupling) gives rise to anomalous Hall / spin Hall responses.

### Chern Number (TKNN)

For a 2D band structure, the **first Chern number** of band $n$ is the integral of the Berry
curvature over the Brillouin zone:

$$
C_n = \frac{1}{2\pi}\int_{\text{BZ}} \Omega_n^{xy}\,d^2k \in \mathbb{Z}
$$

**Proven: integer, topological invariant** — unchanged under any smooth deformation of the
Hamiltonian, and the hallmark of the integer quantum Hall effect (Thouless–Kohmoto–
Nightingale–den Nijs 1982). The **lattice (discrete)** evaluation is the flux of Berry
curvature over plaquettes summing to the integer.

### Winding Number

For a map $g: S^1 \to S^1$ (a closed curve in the complex plane avoiding $0$), the winding
number counts the net number of counterclockwise encirclements of the origin:

$$
w(g) = \frac{1}{2\pi i}\oint \frac{g'(z)}{g(z)}\,dz \in \mathbb{Z}
$$

**Discrete version:** computed from the accumulated argument change of $g(z_k)$ around the
closed polygon. **Skyrmion number** (2D map $\mathbb R^2 \to S^2$, used in chiral magnets) is
the counting of sphere windings — analogous to degree of the map.

---

## Routines

| Routine | Description |
|---------|-------------|
| `berry_phase(psi, dparam, ...)` | $\gamma_n$ around a loop |
| `berry_curvature(H_params, psi, ...)` | $\Omega_n^{ab}$ at a point |
| `chern_number(...)` | $C_n$ integral over BZ |
| `chern_number_discrete(lattice, ...)` | lattice flux version |
| `winding_number(path_real, path_imag)` | $w(g) = \frac1{2\pi i}\oint g'/g$ |
| `winding_number_discrete(points)` | discrete encirclement |
| `skyrmion_number(spin_config)` | degree of the map $\mathbb R^2\to S^2$ |

---

## Usage Examples

### Discrete Chern number of a Landau-level band

```python
from pysicrs import chern_number_discrete

# Lattice Berry flux per plaquette on a 12×12 BZ sample
# Construction: build |u(k)> eigenvectors from H(k) over mesh
C = chern_number_discrete(plaq_flux)
print(int(round(C)))   # 1 for lowest Landau level
```

### Winding number of a closed curve

```python
from pysicrs import winding_number_discrete
import math, cmath

# z(t) = e^{i·2t}: winds twice around the origin
points = [cmath.exp(2j * t) for t in (i / 100 * 2 * math.pi for i in range(101))]
w = winding_number_discrete([p.real for p in points], [p.imag for p in points])
print(w)   # 2
```

### Berry phase around a loop

```python
from pysicrs import berry_phase
# Spin-1/2 in a varying magnetic field: loop that encloses a monopole
# gives γ = 2π × (solid angle)/4π
```

---

## Advantages & Limitations

✅ Direct access to geometric/topological invariants (no ad-hoc derivations)

✅ Both continuum and discrete evaluation for lattice-style problems

✅ Connects naturally to the Gauge module (winding → instanton quantization)

❌ Requires pre-computed eigenstates / Berry connections from outside

❌ No coverage of higher-genus (genus `g`) index theorems besides the winding/skyrmion counts

❌ Reflection of a generically non-abelian (Wilson-loop) connection is out of scope

---

## References

1. Berry, M.V. (1984). "Quantal phase factors accompanying adiabatic changes." *Proc. R. Soc. A* 392:45.
2. Thouless, D.J., Kohmoto, M., Nightingale, M.P. & den Nijs, M. (1982). *Phys. Rev. Lett.* 49:405.
3. Simon, B. (1983). "Holonomy, the quantum adiabatic theorem, and Berry's phase." *Phys. Rev. Lett.* 51:2167.
4. Göckeler, M. & Schücker, T. (1990). *Differential Geometry, Gauge Theories, and Gravity*. Cambridge.

---

## Related Topics

- [Gauge](gauge.md) – winding number quantizes the Yang–Mills instanton action
- [Quantum](quantum.md) – Berry phases appear in adiabatic two-level problems
- [General Relativity](general_relativity.md) – same differential-geometric language