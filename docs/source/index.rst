.. Pysic-rs documentation master file

Pysic-rs Documentation
======================

**High-performance mathematical physics engine in Rust with Python bindings**

.. image:: https://img.shields.io/badge/version-0.1.0-blue.svg
   :target: https://github.com/ThotDjehuty/pysic-rs/releases
   :alt: Version

.. image:: https://img.shields.io/badge/license-MIT-green.svg
   :target: https://github.com/ThotDjehuty/pysic-rs/blob/main/LICENSE
   :alt: License

Pysic-rs provides CPU-only, textbook-validated implementations of standard mathematical
physics, built entirely in Rust and exposed to Python through PyO3. Every algorithm is a
canonical result from the established literature — nothing experimental, nothing
research-pending. This is the physics companion to
`Optimiz-rs <https://optimiz-rs.readthedocs.io/>`_, which supplies the probabilistic and
statistical machinery (MCMC, HMM, SDE, Kalman filtering). Where you need physics *and*
statistics, the two libraries interoperate.

**🗂 Module Map** — Pysic-rs covers thirteen foundational areas:

.. raw:: html

   <div class="grid cards" style="margin: 1.5rem 0;">
     <div style="padding: 1rem; border-radius: 0.75rem; border: 1px solid #c7d2fe; background: #eef2ff; margin-bottom: 0.5rem;">
       <b>🧮 Special Functions</b> — Γ, B, J_n, P_ℓ, Y_ℓ^m, erf, ζ · <b>Linear Algebra</b> — Cholesky, LU, inverse, tensor ops · <b>Calculus</b> — gradient, Hessian, Simpson, Romberg
     </div>
     <div style="padding: 1rem; border-radius: 0.75rem; border: 1px solid #c7d2fe; background: #eef2ff; margin-bottom: 0.5rem;">
       <b>⏱ Equations</b> — RK4/RK45, symplectic integrators · Schrödinger/Dirac split-step, heat, wave FDTD, Poisson, Maxwell 3D FDTD
     </div>
     <div style="padding: 1rem; border-radius: 0.75rem; border: 1px solid #c7d2fe; background: #eef2ff; margin-bottom: 0.5rem;">
       <b>⚛️ Quantum</b> — Pauli algebra, density matrices, coherent states, path integrals, propagators · <b>🌌 GR</b> — Schwarzschild/Kerr/FLRW, Christoffel, geodesics, ADM
     </div>
     <div style="padding: 1rem; border-radius: 0.75rem; border: 1px solid #c7d2fe; background: #eef2ff; margin-bottom: 0.5rem;">
       <b>🎯 Gauge & Topology</b> — SU(2)/SU(3), Yang–Mills, Berry phase, Chern numbers · <b>⚙️ Classical</b> — Hamiltonian, Lagrangian, rigid body, Euler angles
     </div>
     <div style="padding: 1rem; border-radius: 0.75rem; border: 1px solid #c7d2fe; background: #eef2ff; margin-bottom: 0.5rem;">
       <b>⚡ EM & Casimir</b> — Green's functions, radiation formulas, parallel-plate/sphere/Lifshitz Casimir energy
     </div>
   </div>

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   getting-started
   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: Examples

   examples

.. toctree::
   :maxdepth: 2
   :caption: Algorithms

   algorithms/special_functions
   algorithms/linalg
   algorithms/calculus
   algorithms/ode
   algorithms/pde
   algorithms/fourier
   algorithms/quantum
   algorithms/general_relativity
   algorithms/gauge
   algorithms/classical
   algorithms/topology
   algorithms/em
   algorithms/casimir

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/special_functions
   api/linalg
   api/calculus
   api/ode
   api/pde
   api/fourier
   api/quantum
   api/general_relativity
   api/gauge
   api/classical
   api/topology
   api/em
   api/casimir

.. toctree::
   :maxdepth: 1
   :caption: Theory

   theory/mathematical_foundations

.. toctree::
   :maxdepth: 1
   :caption: Advanced

   contributing
   changelog

Features
--------

**Modules Included:**

- **Special Functions**: Gamma (Lanczos), Beta, Bessel J/Y, Legendre P/P_ℓ^m, spherical harmonics, Chebyshev T/U, Airy, erf/erfc, exponential integral, Riemann zeta
- **Linear Algebra**: Cholesky, LU with pivoting, Gauss–Jordan inverse, determinant, tensor raise/lower, metric signature, Lie bracket
- **Numerical Calculus**: gradient, Hessian, Jacobian, trapezoid/Simpson/Gauss–Legendre/Romberg integration, linear & cubic-spline interpolation
- **ODE Solvers**: RK4, RK45 (Dormand–Prince), backward Euler, Crank–Nicolson, leapfrog, velocity Verlet, Yoshida symplectic, Hamiltonian integration
- **PDE Solvers**: Schrödinger (split-step Fourier), Dirac, heat (Crank–Nicolson + ADI), wave (FDTD 1D/2D), Poisson (FFT + SOR), Maxwell (3D Yee)
- **Fourier Analysis**: FFT/IFFT, power spectral density, Welch PSD, spectral derivatives
- **Quantum Mechanics**: Pauli matrices, spin operators, density matrices, partial trace, coherent/number/squeeze states, Rayleigh–Schrödinger perturbation, path integrals, Feynman propagators
- **General Relativity**: Schwarzschild, Kerr, Kerr–Newman, FLRW, Minkowski; Christoffel, Riemann, Ricci, Einstein tensors; geodesics; ADM 3+1; stress–energy
- **Gauge Theory**: SU(2)/SU(3) structure constants, Gell-Mann matrices, connections, covariant derivative, Yang–Mills action/EOM, instantons ('t Hooft)
- **Classical Mechanics**: Hamiltonian, Lagrange/Euler–Lagrange, Poisson bracket, rigid-body (Euler equations, quaternions, Euler angles)
- **Topology**: Berry phase & curvature, Chern/TKNN numbers, winding numbers, skyrmion number
- **Electromagnetism**: Green's functions (static, retarded, Helmholtz, dyadic), radiation formulas (Larmor, dipole, cyclotron, synchrotron, Thomson, Compton)
- **Casimir Effect**: parallel plates, sphere, cylinders, ζ-regularization, finite-T, Lifshitz, Polder, van der Waals

**Performance:**

- Pure Rust core with zero-cost abstractions
- Python bindings via PyO3 (abi3) — no NumPy dependency at the boundary
- Optional Rayon parallelism
- 50–100× faster than pure Python physics code

Quick Example
-------------

.. code-block:: python

    from pysicrs import constants, schwarzschild_metric, rk4_solve

    # Physical constants
    print(constants()["c"])      # 2.997925e+08

    # Schwarzschild metric at r = 10M
    g = schwarzschild_metric(10.0, 1.0)
    print(g[0, 0])     # -(1 - 2/10) = -0.8

    # Integrate the harmonic oscillator with RK4
    def osc(t, y):
        return [y[1], -y[0]]

    times, traj = rk4_solve(osc, [1.0, 0.0], 0.0, 10.0, 10000)
    print(traj[-1][0])  # cos(10) ≈ -0.8391

Installation
------------

From source:

.. code-block:: bash

    # Clone repository
    git clone https://github.com/ThotDjehuty/pysic-rs.git
    cd pysic-rs

    # Build and install
    pip install maturin
    maturin develop --release

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`