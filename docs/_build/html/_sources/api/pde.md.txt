# API Reference: PDE

The Python surface exposes PDE solvers at top level (`from pysicrs import ...`). See also
[Algorithms: PDE Solvers](../algorithms/pde.md).

```python
import pysicrs
```

## Quantum

### `schrodinger_split_step_1d`

```python
schrodinger_split_step_1d(psi0_real, psi0_imag, v, dx, dt, n_steps,
                          mass=1.0, hbar=1.0)
    -> (real: list[float], imag: list[float])
```

Split-step Fourier time evolution of the 1D TDSE. `psi0_real`, `psi0_imag`, `v` are
equal-length arrays on the spatial grid.

### `schrodinger_eigen_1d`

```python
schrodinger_eigen_1d(v, dx, n_states, mass=1.0, hbar=1.0)
    -> (energies: list[float], states: list[list[float]])
```

Imaginary-time propagation to bound-state energies/eigenstates.

### `dirac_split_step_1d`

```python
dirac_split_step_1d(psi0_real, psi0_imag, v, dx, dt, n_steps,
                    c_speed=1.0, mass=1.0, hbar=1.0)
    -> (real: list[float], imag: list[float])
```

Split-step Fourier 1D Dirac equation with scalar potential.

## Parabolic / hyperbolic / elliptic

### `heat_crank_nicolson_1d` / `heat_crank_nicolson_2d`

```python
heat_crank_nicolson_1d(u0, dx, dt, alpha, n_steps) -> list[float]
heat_crank_nicolson_2d(u0, nx, ny, dx, dy, dt, alpha, n_steps) -> list[list[float]]
```

A-stable Crank–Nicolson (1D) / ADI (2D).

### `wave_fdtd_1d` / `wave_fdtd_2d`

```python
wave_fdtd_1d(u0, v0, dx, dt, v, n_steps) -> list[float]
wave_fdtd_2d(u0, v0, nx, ny, dx, dy, dt, v, n_steps) -> list[list[float]]
```

Leapfrog FDTD. Respect the CFL ratio $\frac{v\Delta t}{\Delta x} \le 1$.

### `poisson_fft_2d`

```python
poisson_fft_2d(f, nx, ny) -> list[list[float]]
```

Spectral (FFT) Poisson solve on a periodic $nx\times ny$ grid.

## Rust-only (not yet bound)

| function | purpose |
|----------|---------|
| `poisson_sor_2d` | SOR Poisson solve |
| `maxwell_step_3d` / `maxwell_max_dt` | 3D Yee-grid FDTD |