Example Notebooks
=================

All notebooks are **pre-executed** (kernel ``rhftlab``, mode synthetic —
closed-form formulas / analytic validations), validated end-to-end with
**0 errors** and **≥1 figure** per code cell. Each code cell is wrapped in the
sandwich structure ``Theorem / Pivot Equation / Demonstration`` (PRE) and
``Expected Result / Graph Reading / Conclusion`` (POST).

.. hint::
   Notebooks are rendered here via nbsphinx (executed notebook rendering).
   To re-run locally:

   .. code-block:: bash

      conda activate rhftlab
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

Catalogue
---------

.. list-table::
   :widths: 12 55 33
   :header-rows: 1

   * - Notebook
     - Content
     - Verified Points / CONSTAT
   * - ``00_quickstart_and_constants.ipynb``
     - API overview, CODATA constants
     - 15 constants dict, error < 1e-6
   * - ``01_ode_solvers.ipynb``
     - RK4, backward Euler, leapfrog, RK45
     - RK4 order 4; backward_euler CONSTAT (fixed-point iteration, h·λ<1 strict);
       leapfrog CONSTAT (2nd half-kick defective); rk45 stub None
   * - ``02_pde_heat.ipynb``
     - Heat 1D (Crank–Nicolson), 2D (ADI)
     - Conservation, positivity, boundedness (mass not conserved in Dirichlet → assert adapted)
   * - ``03_pde_wave.ipynb``
     - Wave FDTD 1D/2D
     - Matches independent numpy reference (diff 0.0); CFL ≤ 1/√2
   * - ``04_pde_poisson.ipynb``
     - Poisson FFT 2D
     - CONSTAT: sol/u ratio not constant (233–641)
   * - ``05_pde_quantum.ipynb``
     - Schrödinger split-step, eigenstates, Dirac
     - Norm conserved (1±1e-9); eigen CONSTAT (energy shifted); Dirac CONSTAT (norm ×4/step)
   * - ``06_fourier_analysis.ipynb``
     - FFT/IFFT, Parseval, spectral derivative, Welch
     - Even orders exact; order 1 ≈ 0 (CONSTAT); Welch normalization documented
   * - ``07_quantum_states.ipynb``
     - Pauli/su(2), Fock/coherent states, density
     - CONSTAT: σ_y zero (imaginary parts lost); f_123 = 1; coherent state Poisson
   * - ``08_propagators.ipynb``
     - Feynman propagator, free propagator
     - D = i/(p²−m²+iε); exponential decay e^{-mr} and r⁻² (m=0)
   * - ``09_gauge_theory.ipynb``
     - su(3) structure constants, instanton action
     - CONSTAT su(3): table not antisymmetric (f_147 = −1/2); instanton 8π²/g² exact
   * - ``10_topology.ipynb``
     - Berry, Chern, winding, skyrmion
     - γ=−π equator; chern = ΣF dk²/2π; W integer; skyrmion stub (CONSTAT)
   * - ``11_general_relativity.ipynb``
     - Schwarzschild, Christoffel, Ricci/Einstein, ADM
     - Metric/Christoffel exact; CONSTAT vs literature: Ricci −2/r² in "vacuum";
       ADM constraints (16πGρ convention)
   * - ``12_classical_mechanics.ipynb``
     - Euler rotations, Euler equations, inertia
     - Orthogonality/det = 1; ω̇ axial zero; 3D inertia tensor exact
   * - ``13_electromagnetism_casimir.ipynb``
     - Green, Coulomb, Larmor, dipole, Compton, Casimir, Polder
     - Larmor/Compton/Casimir/Polder exact; dipole CONSTAT (1/c² factor on E/A, F/A matches formulas)