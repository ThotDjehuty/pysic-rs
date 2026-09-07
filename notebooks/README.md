# Pysic-rs — Exercices notebookés

Suite générique de 14 notebooks (kernel `rhftlab`), exécutés de bout en bout
(0 erreur, ≥ 1 figure par cellule de code), couvrant l'ensemble de l'API
`pysicrs` (~48 fonctions).

Chaque cellule de code est enveloppée dans une structure sandwich :
`Théorème / Modèle utilisé` + `Équation pivot` + `Démonstration` (PRE) puis
`Résultat attendu` + `Lecture du graphique` + `Conclusion` (POST). Les écarts
binding vs littérature sont marqués **CONSTAT** (constat d'écart, sans assert
failing).

| # | Notebook | Contenu | Statut |
|---|----------|---------|--------|
| 00 | `00_quickstart_and_constants.ipynb` | Panorama API, constantes CODATA | ✓ 7 cells · 0 err · 1 fig |
| 01 | `01_ode_solvers.ipynb` | RK4, backward Euler, leapfrog, RK45 | ✓ 19 · 0 · 4 |
| 02 | `02_pde_heat.ipynb` | Chaleur 1D CN / 2D ADI | ✓ 10 · 0 · 2 |
| 03 | `03_pde_wave.ipynb` | Ondes FDTD 1D/2D | ✓ 10 · 0 · 2 |
| 04 | `04_pde_poisson.ipynb` | Poisson FFT 2D | ✓ 7 · 0 · 1 |
| 05 | `05_pde_quantum.ipynb` | Schrödinger split-step, eigen, Dirac | ✓ 16 · 0 · 4 |
| 06 | `06_fourier_analysis.ipynb` | FFT, Parseval, dérivée spectrale, Welch | ✓ 13 · 0 · 3 |
| 07 | `07_quantum_states.ipynb` | Pauli/su(2), états nombre/cohérents, densité | ✓ 13 · 0 · 3 |
| 08 | `08_propagators.ipynb` | Propagateur de Feynman, propagateur libre | ✓ 10 · 0 · 2 |
| 09 | `09_gauge_theory.ipynb` | Constantes su(3), action instanton | ✓ 10 · 0 · 2 |
| 10 | `10_topology.ipynb` | Berry, Chern, enroulement, skyrmion | ✓ 16 · 0 · 4 |
| 11 | `11_general_relativity.ipynb` | Schwarzschild, Christoffel, Ricci/Einstein, ADM | ✓ 16 · 0 · 4 |
| 12 | `12_classical_mechanics.ipynb` | Rotations d'Euler, équations d'Euler, inertie | ✓ 13 · 0 · 3 |
| 13 | `13_electromagnetism_casimir.ipynb` | Green, Coulomb, Larmor, dipôle, Compton, Casimir, Polder | ✓ 19 · 0 · 5 |

Total : **208 cellules, 0 erreur, 49 figures**.

## Constats documentés (binding ≠ littérature)

- **01 ODE** : `backward_euler` — itération de point fixe, converge seulement si
  `h·λ < 1` (pas inconditionnellement A-stable) ; `leapfrog_integrate` —
  énergie qui dérive (2ᵉ demi-impulsion défectueuse, cf. `src/ode/symplectic.rs`) ;
  `rk45_solve` — stub renvoyant `None`.
- **04 Poisson** : rapport `sol/u` non constant (233–641) en intérieur.
- **05 Quantique** : `schrodinger_eigen_1d` énergies décalées (e.g. ∞ well :
  1301 vs 0.049) ; `dirac_split_step_1d` norme ×4 par pas (probabilité non
  unitaire).
- **06 Fourier** : dérivée spectrale d'ordre 1 ≈ 0 (seuls les ordres pairs
  marchent) ; fréquence Welch en cycles/échantillon.
- **07 États** : `pauli_matrices` renvoie parties réelles uniquement → σ_y nulle ;
  `su2_structure_constants` f_123 = 1 (générateurs σᵢ/2).
- **09 Jauge** : table su(3) non totalement antisymétrique (`f_147` = −1/2 au
  lieu de ±1/2) — bug de remplissage Rust.
- **10 Topologie** : `skyrmion_number` = stub (0.0).
- **11 GR** : `ricci_scalar`/`einstein_tensor` Schwartzschild codés en dur M=1
  → R ≈ −2/r² en « vide » (au lieu de 0) ; contraintes ADM en convention
  16πGρ / 8πG j.
- **13 EM** : `dipole_radiation` divise par `12πc³` au lieu de `12πc`
  (facteur c²).

## Exécution

```bash
conda activate rhftlab
jupyter nbconvert --to notebook --inplace --execute notebooks/*.ipynb
```

Rendu dans la doc : voir `docs/source/examples.rst` (nbsphinx, notebooks
pré-exécutés).