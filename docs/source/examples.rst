Example Notebooks
=================

Fifteen pre-executed notebooks covering the library's mathematical-physics
domains. Each serves as both a validation suite and a worked example: every
code cell prints labelled output, produces at least one figure, and compares
the library against a reference known **independently** of the routine under
test.

Each code cell is wrapped in the sandwich structure ``Theorem / Pivot Equation
/ Demonstration`` (PRE) and ``Expected Result / Graph Reading / Conclusion``
(POST), so the physics being checked is stated before the number appears.

.. hint::
   Notebooks are rendered here via nbsphinx. To re-run them locally:

   .. code-block:: bash

      pip install pysic-rs numpy scipy matplotlib jupyter
      jupyter nbconvert --to notebook --inplace --execute notebooks/*.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Notebooks

   notebooks/00_quickstart_and_constants.ipynb
   notebooks/01_ode_solvers.ipynb
   notebooks/02_pde_heat.ipynb
   notebooks/03_pde_wave.ipynb
   notebooks/04_pde_poisson.ipynb
   notebooks/05_pde_quantum.ipynb
   notebooks/06_fourier_analysis.ipynb
   notebooks/07_quantum_states.ipynb
   notebooks/08_propagators.ipynb
   notebooks/09_gauge_theory.ipynb
   notebooks/10_topology.ipynb
   notebooks/11_general_relativity.ipynb
   notebooks/12_classical_mechanics.ipynb
   notebooks/13_electromagnetism_casimir.ipynb
   notebooks/15_dyson_schwinger.ipynb

How these notebooks grade the library
-------------------------------------

Three classes of invariant do most of the work, because each is falsifiable
*without knowing the answer*:

**Unitarity.**
  A propagator built as a product of exponentials of Hermitian operators is
  unitary, so :math:`\lVert\psi(t)\rVert^2` is an exact invariant of the
  discrete scheme — not merely to truncation order. This is a binary test with
  no tolerance to tune: either the norm holds to round-off or the operator is
  wrong.

**Vanishing on a known solution.**
  Schwarzschild is a vacuum solution, so its Ricci scalar is identically zero
  at every radius *and every polar angle*. Any non-zero result is a defect,
  with no discretisation excuse available beyond the finite-difference
  residual.

**Integrality.**
  Topological invariants are integers. A winding number that returns
  :math:`2.9997` is not imprecise; it is broken.

Where a convergence order is estimated from a finite set of grid refinements,
the notebooks report a percentile bootstrap interval on that estimate rather
than quoting a single number.

Catalogue
---------

.. list-table::
   :widths: 12 52 36
   :header-rows: 1

   * - Notebook
     - Content
     - Verified against
   * - ``00_quickstart_and_constants``
     - API overview, CODATA 2018 constants
     - 15 constants, max relative error < 1e-6
   * - ``01_ode_solvers``
     - RK4, backward Euler, leapfrog, RK45 (Dormand–Prince)
     - RK4 order 4; RK45 reaches :math:`e^{-1}` with the error tracking ``rtol``
   * - ``02_pde_heat``
     - Heat 1D (Crank–Nicolson), 2D (ADI)
     - positivity, boundedness, conservation under the stated boundary condition
   * - ``03_pde_wave``
     - Wave FDTD 1D/2D
     - independent NumPy reference (difference 0.0); CFL :math:`\le 1/\sqrt2`
   * - ``04_pde_poisson``
     - Poisson FFT 2D
     - substitution into :math:`\nabla^2 u = f`
   * - ``05_pde_quantum``
     - Schrödinger split-step, eigenstates, Dirac
     - norm conserved to 1e-14; spectrum vs :math:`\frac{2\hbar^2}{m\Delta x^2}\sin^2\frac{k\pi}{2(N+1)}`
   * - ``06_fourier_analysis``
     - FFT/IFFT, Parseval, spectral derivative, Welch PSD
     - Parseval identity; spectral derivative of a known mode
   * - ``07_quantum_states``
     - Pauli / su(2), Fock and coherent states, density matrices
     - Clifford algebra :math:`\{\sigma_i,\sigma_j\}=2\delta_{ij}`; coherent-state Poisson statistics
   * - ``08_propagators``
     - Feynman propagator, free propagator
     - :math:`D = i/(p^2-m^2+i\varepsilon)`; decay :math:`e^{-mr}` and :math:`r^{-2}` at :math:`m=0`
   * - ``09_gauge_theory``
     - su(2)/su(3) structure constants, instanton action
     - Jacobi identity; instanton action :math:`8\pi^2/g^2` exact
   * - ``10_topology``
     - Berry phase, Chern numbers, winding, skyrmion charge
     - :math:`\gamma=-\pi` at the equator; winding exactly integer
   * - ``11_general_relativity``
     - Schwarzschild metric, Christoffel, Ricci/Einstein, ADM
     - vacuum :math:`R \equiv 0`; flat-slice Hamiltonian constraint :math:`H=6c^2`
   * - ``12_classical_mechanics``
     - Euler rotations, Euler equations, inertia tensor
     - orthogonality and :math:`\det=1`; axial :math:`\dot\omega = 0`
   * - ``13_electromagnetism_casimir``
     - Green's functions, Coulomb, Larmor, dipole, Compton, Casimir, Polder
     - Larmor, Compton, Casimir and Polder closed forms
   * - ``15_dyson_schwinger``
     - Lippmann–Schwinger / Dyson equation, rainbow gap equation
     - algebraic residual at round-off; Born limit; geometric Neumann decay

.. note::
   These notebooks were the instrument that exposed eight incorrect routines in
   release 0.1.0 — among them a Bessel function returning :math:`+50.6` for a
   quantity bounded by 1, an Airy function returning ``NaN`` at the origin, and
   a vacuum spacetime reporting non-zero Ricci curvature. All eight are
   repaired in 0.2.0 and each is now pinned by a test that fails on the old
   code. See the :doc:`changelog` for the post-mortem.

Where this physics is used
--------------------------

1. **Quantum computing** — Pauli matrices and density matrices underpin
   quantum error correction and algorithm simulation.

2. **General relativity** — Schwarzschild and Kerr metrics support black-hole
   and gravitational-wave modelling; the ADM constraints are the starting
   point for numerical-relativity evolution.

3. **Electromagnetism and photonics** — Green's functions and radiation
   formulas are the basis of antenna modelling and electromagnetic
   compatibility analysis.

4. **Condensed matter** — gauge theory and topological invariants describe
   topological insulators, quantum Hall systems and superconductors. The
   library's logo is the Hofstadter butterfly it computes.

5. **Particle and nuclear physics** — su(3) structure constants, instantons
   and the Dyson–Schwinger equations are central to QCD and non-perturbative
   phenomena.

6. **Quantum gravity** — Wheeler–DeWitt minisuperspace, Ashtekar variables and
   Chern–Simons formulations of 2+1 gravity, including the BTZ black hole.
