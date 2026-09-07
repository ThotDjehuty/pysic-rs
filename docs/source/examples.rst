Example Notebooks
=================

All notebooks are **pre-executed** (kernel ``rhftlab``, mode synthétique —
formules fermées / vérifications analytiques), validated end-to-end with
**0 erreur** and **≥1 figure** per code cell. Each code cell is wrapped in the
sandwich structure ``Théorème / Équation pivot / Démonstration`` (PRE) and
``Résultat attendu / Lecture du graphique / Conclusion`` (POST).

.. hint::
   Les notebooks sont présentés ici via nbsphinx (rendu du notebook exécuté).
   Pour les relancer localement :

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
     - Contenu
     - Points vérifiés / CONSTAT
   * - ``00_quickstart_and_constants.ipynb``
     - Panorama de l'API, constantes CODATA
     - dict 15 constantes, écart < 1e-6
   * - ``01_ode_solvers.ipynb``
     - RK4, backward Euler, leapfrog, RK45
     - RK4 ordre 4 ; backward_Euler CONSTAT (itération de point fixe, h·λ·<1 strict) ;
       leapfrog CONSTAT (2e demi-coup défectueux) ; rk45 stub None
   * - ``02_pde_heat.ipynb``
     - Chaleur 1D (Crank–Nicolson), 2D (ADI)
     - Conservation, positivité, bornitude (masse non conservée en Dirichlet → assert adapté)
   * - ``03_pde_wave.ipynb``
     - Ondes FDTD 1D/2D
     - Identiques à une référence numpy indépendante (écart 0.0) ; CFL ≤ 1/√2
   * - ``04_pde_poisson.ipynb``
     - Poisson FFT 2D
     - CONSTAT : rapport sol/u non constant (233–641)
   * - ``05_pde_quantum.ipynb``
     - Schrödinger split-step, états propres, Dirac
     - Norme conservée (1±1e-9) ; eigen CONSTAT (énergie décalées) ; Dirac CONSTAT (norme ×4/pas)
   * - ``06_fourier_analysis.ipynb``
     - FFT/IFFT, Parseval, dérivée spectrale, Welch
     - Ordres pairs exacts ; ordre 1 ≈ 0 (CONSTAT) ; normalisation Welch documentée
   * - ``07_quantum_states.ipynb``
     - Pauli/su(2), états nombre/cohérents, densité
     - CONSTAT : σ_y nulle (parties imaginaires perdues) ; f_123 = 1 ; Poisson des états cohérents
   * - ``08_propagators.ipynb``
     - Propagateur de Feynman, propagateur libre
     - D = i/(p²−m²+iε) ; décroissance e^{-mr} et r⁻² (m=0)
   * - ``09_gauge_theory.ipynb``
     - Constantes su(3), action instanton
     - CONSTAT su(3) : table non antisymétrique (f_147 = −1/2) ; instanton 8π²/g² exact
   * - ``10_topology.ipynb``
     - Berry, Chern, enroulement, skyrmion
     - γ=−π équateur ; chern = ΣF dk²/2π ; W entier ; skyrmion stub (CONSTAT)
   * - ``11_general_relativity.ipynb``
     - Schwarzschild, Christoffel, Ricci/Einstein, ADM
     - Métrique/Christoffel exactes ; CONSTAT selon la littérature : Ricci −2/r² en « vide » ;
       contraintes ADM (convention 16πGρ)
   * - ``12_classical_mechanics.ipynb``
     - Rotations d'Euler, équations d'Euler, inertie
     - Orthogonalité/det = 1 ; ω̇ axial nul ; tenseur d'inertie 3D exact
   * - ``13_electromagnetism_casimir.ipynb``
     - Green, Coulomb, Larmor, dipôle, Compton, Casimir, Polder
     - Larmor/Compton/Casimir/Polder exacts ; dipôle CONSTAT (facteur 1/c² sur E/A, F/A conforme aux formules)