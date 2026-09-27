//! Correctness tests for the gauge module, written against physics-mandated
//! invariants rather than against the current implementation.
//!
//! These follow the same discipline as the solver fixes in 09bbb47 and the
//! special-function repairs in 4d95e5a: assert the closed-form/structural
//! property first, then make the implementation satisfy it.

#[cfg(test)]
mod gauge_invariants {
    use pysic::gauge::curvature::field_strength;
    use pysic::gauge::yang_mills::{instanton_action, thooft_symbol, yang_mills_action, yang_mills_eom};

    /// F_μν is a 2-form: F_μν = -F_νμ for ANY gauge field, always.
    #[test]
    #[ignore = "KNOWN BUG (gauge repair pending): field_strength differentiates \
                along stride n_components for μ but stride 1 for ν, so it is not \
                antisymmetric. Run with `cargo test -- --ignored`."]
    fn field_strength_is_antisymmetric() {
        let n_components = 2;
        let n_points = 8;
        let dx = 0.1;

        // An arbitrary, non-symmetric A_μ(x) so no accidental cancellation.
        let a_field: Vec<Vec<f64>> = (0..n_points)
            .map(|p| (0..n_components).map(|c| 0.3 * p as f64 + 1.7 * c as f64 + 0.11 * (p * c) as f64).collect())
            .collect();

        let f = field_strength(&a_field, dx, n_points, n_components);

        for p in 0..n_points {
            for mu in 0..n_components {
                for nu in 0..n_components {
                    let sum = f[p][mu][nu] + f[p][nu][mu];
                    assert!(
                        sum.abs() < 1e-12,
                        "F_{mu}{nu} + F_{nu}{mu} = {sum:.6e} at point {p} (must vanish: F is antisymmetric)"
                    );
                }
            }
        }
    }

    /// The Yang-Mills action S = -¼∫F^a_μν F^{a μν} is QUADRATIC in F.
    /// Scaling F → λF must scale S → λ²S.
    #[test]
    #[ignore = "KNOWN BUG (gauge repair pending): the action contracts a single F \
                instead of F·F, so S(2F)/S(F) = 2 rather than 4."]
    fn yang_mills_action_is_quadratic_in_field_strength() {
        let n = 2;
        let dx = 0.1;
        let g_inv: Vec<Vec<f64>> = (0..n)
            .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();

        // f_tensor is indexed [mu*n+nu][rho*n+sigma], so it needs to be n² × n².
        let f: Vec<Vec<f64>> = (0..n * n)
            .map(|i| (0..n * n).map(|j| 0.5 + 0.25 * i as f64 - 0.125 * j as f64).collect())
            .collect();
        let f2: Vec<Vec<f64>> = f.iter().map(|row| row.iter().map(|v| 2.0 * v).collect()).collect();

        let s1 = yang_mills_action(&f, &g_inv, n, dx);
        let s2 = yang_mills_action(&f2, &g_inv, n, dx);

        // Guard against a degenerate all-zero action making the test vacuous.
        assert!(s1.abs() > 1e-30, "action vanished identically; test would be vacuous");

        let ratio = s2 / s1;
        assert!(
            (ratio - 4.0).abs() < 1e-9,
            "S(2F)/S(F) = {ratio:.6} but must be 4 (action is quadratic in F). \
             A ratio of 2 means the action is linear in F."
        );
    }

    /// The non-abelian term ig[A_μ, F^μν] must actually contribute.
    /// Changing the coupling g on a non-commuting configuration must change the EOM.
    #[test]
    #[ignore = "KNOWN BUG (gauge repair pending): the commutator is built from f64 \
                scalars so it is identically zero, and the index arithmetic panics \
                out of bounds on a plain n=2 input."]
    fn yang_mills_eom_depends_on_coupling() {
        let n = 2;
        let dx = 0.1;
        let f: Vec<Vec<f64>> = (0..n * n)
            .map(|i| (0..n * n).map(|j| 0.7 + 0.3 * i as f64 - 0.2 * j as f64).collect())
            .collect();
        let a: Vec<Vec<f64>> = (0..n).map(|p| (0..n).map(|c| 1.0 + 0.5 * p as f64 + 0.25 * c as f64).collect()).collect();

        let weak = yang_mills_eom(&f, &a, 0.0, dx, n);
        let strong = yang_mills_eom(&f, &a, 10.0, dx, n);

        let differs = weak.iter().zip(&strong).any(|(w, s)| (w - s).abs() > 1e-12);
        assert!(
            differs,
            "EOM is identical at g=0 and g=10 ({weak:?} vs {strong:?}) — the commutator \
             term is identically zero, so the theory is abelian by construction."
        );
    }

    /// 't Hooft symbol η^a_{μν}: antisymmetric in μν, and defined over the FULL
    /// 4D index range (a = 1..3, μ,ν = 0..3), since instantons live in 4D.
    #[test]
    #[ignore = "KNOWN BUG (gauge repair pending): η^a_{μν} has no time components \
                and index 3 is unreachable, so it cannot build a 4D self-dual \
                instanton."]
    fn thooft_symbol_is_antisymmetric_and_covers_four_dimensions() {
        for a in 0..3 {
            for mu in 0..4 {
                for nu in 0..4 {
                    let s = thooft_symbol(a, mu, nu) + thooft_symbol(a, nu, mu);
                    assert!(
                        s.abs() < 1e-12,
                        "η^{a}_{mu}{nu} + η^{a}_{nu}{mu} = {s} (must vanish)"
                    );
                }
            }
        }

        // η^a_{0i} = -δ_{ai} (up to convention) must be non-zero: these are exactly
        // the components that make the instanton self-dual.
        let mut nonzero_with_time_index = 0;
        for a in 0..3 {
            for i in 1..4 {
                if thooft_symbol(a, 0, i).abs() > 1e-12 {
                    nonzero_with_time_index += 1;
                }
            }
        }
        assert!(
            nonzero_with_time_index >= 3,
            "all η^a_{{0i}} vanish — the symbol has no time components, so it cannot \
             build a self-dual instanton in 4D"
        );

        // Index 3 must be reachable at all.
        let mut touches_index_three = false;
        for a in 0..3 {
            for mu in 0..4 {
                if thooft_symbol(a, mu, 3).abs() > 1e-12 || thooft_symbol(a, 3, mu).abs() > 1e-12 {
                    touches_index_three = true;
                }
            }
        }
        assert!(touches_index_three, "η^a_{{μν}} is identically zero whenever an index is 3");
    }

    /// f^{abc} must be totally antisymmetric in all three indices, for every
    /// swap — not merely cyclic. The previous su(3) table filled only cyclic
    /// permutations (and overwrote some entries with contradictory values),
    /// leaving 27 non-zero entries instead of 54.
    #[test]
    fn structure_constants_are_totally_antisymmetric() {
        use pysic::gauge::bundles::{su2_structure_constants, su3_structure_constants};

        let f2 = su2_structure_constants();
        for a in 0..3 {
            for b in 0..3 {
                for c in 0..3 {
                    assert!((f2[a][b][c] + f2[b][a][c]).abs() < 1e-12, "su(2) a<->b at {a}{b}{c}");
                    assert!((f2[a][b][c] + f2[a][c][b]).abs() < 1e-12, "su(2) b<->c at {a}{b}{c}");
                }
            }
        }

        let f3 = su3_structure_constants();
        for a in 0..8 {
            for b in 0..8 {
                for c in 0..8 {
                    assert!(
                        (f3[a][b][c] + f3[b][a][c]).abs() < 1e-12,
                        "su(3) not antisymmetric under a<->b at ({a},{b},{c}): \
                         {} vs {}", f3[a][b][c], f3[b][a][c]
                    );
                    assert!(
                        (f3[a][b][c] + f3[a][c][b]).abs() < 1e-12,
                        "su(3) not antisymmetric under b<->c at ({a},{b},{c})"
                    );
                }
            }
        }

        let nz2 = (0..3).flat_map(|a| (0..3).flat_map(move |b| (0..3).map(move |c| (a, b, c))))
            .filter(|&(a, b, c)| f2[a][b][c].abs() > 1e-12).count();
        let nz3 = (0..8).flat_map(|a| (0..8).flat_map(move |b| (0..8).map(move |c| (a, b, c))))
            .filter(|&(a, b, c)| f3[a][b][c].abs() > 1e-12).count();
        assert_eq!(nz2, 6, "su(2) should have 1 triple × 3! = 6 non-zero entries");
        assert_eq!(nz3, 54, "su(3) should have 9 triples × 3! = 54 non-zero entries");
    }

    /// Every tabulated su(3) value must be present with the right magnitude.
    #[test]
    fn su3_structure_constants_match_the_gell_mann_table() {
        use pysic::gauge::bundles::su3_structure_constants;
        let f = su3_structure_constants();
        let s = 3.0_f64.sqrt() / 2.0;
        // (a, b, c) 1-indexed as usually tabulated, with the expected value.
        let table = [
            (1, 2, 3, 1.0),
            (1, 4, 7, 0.5),
            (1, 5, 6, -0.5),
            (2, 4, 6, 0.5),
            (2, 5, 7, 0.5),
            (3, 4, 5, 0.5),
            (3, 6, 7, -0.5),
            (4, 5, 8, s),
            (6, 7, 8, s),
        ];
        for (a, b, c, want) in table {
            let got = f[a - 1][b - 1][c - 1];
            assert!(
                (got - want).abs() < 1e-12,
                "f^{{{a}{b}{c}}} = {got}, expected {want}"
            );
        }
    }

    /// Jacobi identity: f^{ade}f^{bcd} + f^{bde}f^{cad} + f^{cde}f^{abd} = 0.
    /// This is the statement that the algebra actually closes.
    #[test]
    fn su3_structure_constants_satisfy_jacobi() {
        use pysic::gauge::bundles::su3_structure_constants;
        let f = su3_structure_constants();
        let mut worst = 0.0f64;
        for a in 0..8 {
            for b in 0..8 {
                for c in 0..8 {
                    for e in 0..8 {
                        let mut sum = 0.0;
                        for d in 0..8 {
                            sum += f[a][d][e] * f[b][c][d]
                                + f[b][d][e] * f[c][a][d]
                                + f[c][d][e] * f[a][b][d];
                        }
                        worst = worst.max(sum.abs());
                    }
                }
            }
        }
        assert!(worst < 1e-12, "Jacobi identity violated, worst residual {worst:e}");
    }

    /// Closed form: S_inst = 8π²/g². This one should already hold.
    #[test]
    fn instanton_action_matches_closed_form() {
        for g in [0.5, 1.0, 2.0, 3.3] {
            let expected = 8.0 * std::f64::consts::PI.powi(2) / (g * g);
            let got = instanton_action(g);
            assert!(
                (got - expected).abs() < 1e-12,
                "S_inst({g}) = {got}, expected {expected}"
            );
        }
    }
}
