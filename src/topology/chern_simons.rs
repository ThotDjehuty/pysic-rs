//! Chern–Simons theory: the cubic functional, Gauss linking numbers (knot
//! invariants), and 2+1 gravity as an SL(2,ℝ)×SL(2,ℝ) Chern–Simons theory
//! (the BTZ black hole).
//!
//! Every routine here is graded against a quantity known independently of the
//! routine itself — an exact integer (linking number), an exact algebraic
//! identity (antisymmetry of the cubic form), unimodularity of the holonomy
//! (det = 1, i.e. the flow stays in the group), or a closed form for the BTZ
//! mass and angular momentum. See the unit tests at the foot of this file.

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────
// 1. The Chern–Simons cubic functional on a simplicial 2-skeleton
// ─────────────────────────────────────────────────────────────────────────

/// Discrete Chern–Simons cubic functional on a connection matrix `a`:
///
/// ```text
/// CS(A) = Σ_{i<j<k} ( A_ij A_jk A_ki − A_ik A_kj A_ji ).
/// ```
///
/// This is the lattice image of the cubic term `(2/3) tr(A∧A∧A)` of the
/// continuum Chern–Simons action, summed over the oriented triangles of a
/// simplicial 2-skeleton whose edges carry the scalar connection `A_ij`.
///
/// It is a *heuristic* summary statistic, not a quantised topological
/// invariant, so it is graded only by the exact algebraic identities it must
/// satisfy for every input:
///
///   * `CS(Aᵀ) = −CS(A)`  — reversing every edge orientation flips the sign;
///   * `CS(A) = 0` whenever `A` is symmetric (`A = Aᵀ`).
///
/// Both are falsifiable without knowing the value of `CS(A)` itself.
pub fn chern_simons_triple(a: &[Vec<f64>]) -> f64 {
    let n = a.len();
    let mut total = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            for k in (j + 1)..n {
                total += a[i][j] * a[j][k] * a[k][i] - a[i][k] * a[k][j] * a[j][i];
            }
        }
    }
    total
}

// ─────────────────────────────────────────────────────────────────────────
// 2. Gauss linking number — a knot/link invariant
// ─────────────────────────────────────────────────────────────────────────

/// Gauss linking number of two closed space curves:
///
/// ```text
/// lk = (1/4π) ∮∮ (r₁ − r₂) · (dr₁ × dr₂) / |r₁ − r₂|³.
/// ```
///
/// Each curve is a closed polyline (the last point is joined back to the
/// first). The double integral is evaluated at segment midpoints. The result
/// is a topological invariant and must be an integer: the Hopf link returns
/// ±1, an unlink returns 0, and a partner curve wrapped `n` times returns `n`.
/// We return the raw value, not a rounded one, precisely so that the distance
/// to the nearest integer is visible as evidence.
pub fn linking_number(curve1: &[[f64; 3]], curve2: &[[f64; 3]]) -> f64 {
    let n1 = curve1.len();
    let n2 = curve2.len();
    let mut sum = 0.0;

    for i in 0..n1 {
        let a0 = curve1[i];
        let a1 = curve1[(i + 1) % n1];
        let dr1 = [a1[0] - a0[0], a1[1] - a0[1], a1[2] - a0[2]];
        let m1 = [
            0.5 * (a0[0] + a1[0]),
            0.5 * (a0[1] + a1[1]),
            0.5 * (a0[2] + a1[2]),
        ];

        for j in 0..n2 {
            let b0 = curve2[j];
            let b1 = curve2[(j + 1) % n2];
            let dr2 = [b1[0] - b0[0], b1[1] - b0[1], b1[2] - b0[2]];
            let m2 = [
                0.5 * (b0[0] + b1[0]),
                0.5 * (b0[1] + b1[1]),
                0.5 * (b0[2] + b1[2]),
            ];

            let r = [m1[0] - m2[0], m1[1] - m2[1], m1[2] - m2[2]];
            let rn = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
            if rn < 1e-12 {
                continue;
            }

            // dr1 × dr2
            let cross = [
                dr1[1] * dr2[2] - dr1[2] * dr2[1],
                dr1[2] * dr2[0] - dr1[0] * dr2[2],
                dr1[0] * dr2[1] - dr1[1] * dr2[0],
            ];
            let triple = r[0] * cross[0] + r[1] * cross[1] + r[2] * cross[2];
            sum += triple / (rn * rn * rn);
        }
    }

    sum / (4.0 * PI)
}

// ─────────────────────────────────────────────────────────────────────────
// 3. 2+1 gravity as SL(2,ℝ)×SL(2,ℝ) Chern–Simons — the BTZ black hole
// ─────────────────────────────────────────────────────────────────────────

/// Matrix exponential of a traceless 2×2 real matrix `X = [[a, b], [c, −a]]`,
/// i.e. a general element of the Lie algebra sl(2,ℝ), returned as
/// `[[m00, m01], [m10, m11]]`.
///
/// Closed form via Cayley–Hamilton. With `δ = a² + bc` one has `X² = δ I`, so
///
/// ```text
/// exp X = cosh(√δ) I + (sinh(√δ)/√δ) X        (δ > 0, hyperbolic)
/// exp X = cos(√−δ) I + (sin(√−δ)/√−δ) X       (δ < 0, elliptic)
/// exp X = I + X                                (δ = 0, parabolic)
/// ```
///
/// The result lies in SL(2,ℝ): `det(exp X) = 1` exactly, because
/// `tr X = 0`. That identity is the routine's falsifiable grade.
pub fn sl2_exp(a: f64, b: f64, c: f64) -> [[f64; 2]; 2] {
    let delta = a * a + b * c;
    let (ch, sh_over) = if delta > 1e-15 {
        let s = delta.sqrt();
        (s.cosh(), s.sinh() / s)
    } else if delta < -1e-15 {
        let s = (-delta).sqrt();
        (s.cos(), s.sin() / s)
    } else {
        (1.0, 1.0)
    };
    [
        [ch + sh_over * a, sh_over * b],
        [sh_over * c, ch - sh_over * a],
    ]
}

/// Trace of `exp X` for `X = [[a, b], [c, −a]] ∈ sl(2,ℝ)`.
/// Equals `2 cosh √δ` (hyperbolic), `2 cos √−δ` (elliptic), `2` (parabolic),
/// with `δ = a² + bc`.
pub fn sl2_exp_trace(a: f64, b: f64, c: f64) -> f64 {
    let m = sl2_exp(a, b, c);
    m[0][0] + m[1][1]
}

/// Traces of the two BTZ holonomies around the angular circle.
///
/// Anti-de Sitter 2+1 gravity is a Chern–Simons theory for
/// SL(2,ℝ)_L × SL(2,ℝ)_R. The BTZ black hole is a *flat* connection whose only
/// gauge-invariant data is the holonomy around the non-contractible angular
/// circle. In each sector the holonomy is hyperbolic with
///
/// ```text
/// tr H_± = 2 cosh( π (r₊ ∓ r₋) / ℓ ),
/// ```
///
/// where `r₊, r₋` are the outer/inner horizon radii and `ℓ` the AdS radius.
/// We build each holonomy as a genuine, non-diagonal sl(2,ℝ) exponential (so
/// the 2×2 matrix exp is exercised, not short-circuited) and return
/// `(tr H₊, tr H₋)`.
pub fn btz_holonomy_traces(r_plus: f64, r_minus: f64, ell: f64) -> (f64, f64) {
    let lam_p = PI * (r_plus - r_minus) / ell;
    let lam_m = PI * (r_plus + r_minus) / ell;
    // X = [[0, λ], [λ, 0]] ⇒ δ = λ² ⇒ hyperbolic ⇒ tr exp X = 2 cosh λ.
    (sl2_exp_trace(0.0, lam_p, lam_p), sl2_exp_trace(0.0, lam_m, lam_m))
}

/// BTZ mass and angular momentum from the horizon radii, in units `8G = 1`:
///
/// ```text
/// M = (r₊² + r₋²) / ℓ² ,     J = 2 r₊ r₋ / ℓ.
/// ```
pub fn btz_mass_spin(r_plus: f64, r_minus: f64, ell: f64) -> (f64, f64) {
    let m = (r_plus * r_plus + r_minus * r_minus) / (ell * ell);
    let j = 2.0 * r_plus * r_minus / ell;
    (m, j)
}

/// Recover the horizon radii from the two holonomy traces — the inverse of
/// [`btz_holonomy_traces`] — and return `(r₊, r₋)`. Composing the two is a
/// round trip against a closed form: feed `(r₊, r₋)`, read it back unchanged.
pub fn btz_horizons_from_traces(tr_plus: f64, tr_minus: f64, ell: f64) -> (f64, f64) {
    let lam_p = (tr_plus / 2.0).acosh(); // = π(r₊ − r₋)/ℓ
    let lam_m = (tr_minus / 2.0).acosh(); // = π(r₊ + r₋)/ℓ
    let r_plus = ell * (lam_m + lam_p) / (2.0 * PI);
    let r_minus = ell * (lam_m - lam_p) / (2.0 * PI);
    (r_plus, r_minus)
}

// ─────────────────────────────────────────────────────────────────────────
// Tests — each graded against an independently known quantity
// ─────────────────────────────────────────────────────────────────────────
#[cfg(test)]
mod tests {
    use super::*;

    fn circle_xy(radius: f64, center: [f64; 3], n: usize) -> Vec<[f64; 3]> {
        (0..n)
            .map(|i| {
                let t = 2.0 * PI * i as f64 / n as f64;
                [center[0] + radius * t.cos(), center[1] + radius * t.sin(), center[2]]
            })
            .collect()
    }

    /// Partner circle in the x–z plane, optionally wrapped `wrap` times.
    fn circle_xz_wrapped(radius: f64, center: [f64; 3], n: usize, wrap: usize) -> Vec<[f64; 3]> {
        (0..n)
            .map(|i| {
                let t = 2.0 * PI * wrap as f64 * i as f64 / n as f64;
                [center[0] + radius * t.cos(), center[1], center[2] + radius * t.sin()]
            })
            .collect()
    }

    #[test]
    fn linking_hopf_is_plus_minus_one() {
        // Two interlocked unit circles: the Hopf link, lk = ±1.
        let c1 = circle_xy(1.0, [0.0, 0.0, 0.0], 400);
        let c2 = circle_xz_wrapped(1.0, [1.0, 0.0, 0.0], 400, 1);
        let lk = linking_number(&c1, &c2);
        assert!((lk.abs() - 1.0).abs() < 0.03, "Hopf link lk = {lk}, expected ±1");
    }

    #[test]
    fn linking_unlink_is_zero() {
        // Two far-apart circles: the unlink, lk = 0.
        let c1 = circle_xy(1.0, [0.0, 0.0, 0.0], 400);
        let c2 = circle_xz_wrapped(1.0, [12.0, 0.0, 0.0], 400, 1);
        let lk = linking_number(&c1, &c2);
        assert!(lk.abs() < 0.01, "unlink lk = {lk}, expected 0");
    }

    #[test]
    fn linking_scales_with_wrap_number() {
        // Partner wrapped n times through the disk: lk = n (degree of the Gauss map).
        let c1 = circle_xy(1.0, [0.0, 0.0, 0.0], 800);
        for n in 2..=3usize {
            let c2 = circle_xz_wrapped(1.0, [1.0, 0.0, 0.0], 400 * n, n);
            let lk = linking_number(&c1, &c2);
            assert!(
                (lk.abs() - n as f64).abs() < 0.05,
                "wrap {n}: lk = {lk}, expected {n}"
            );
        }
    }

    #[test]
    fn cubic_vanishes_on_symmetric_connection() {
        // CS(A) = 0 for A = Aᵀ.
        let a = vec![
            vec![0.0, 0.7, -0.2, 0.5],
            vec![0.7, 0.0, 0.9, -0.3],
            vec![-0.2, 0.9, 0.0, 0.4],
            vec![0.5, -0.3, 0.4, 0.0],
        ];
        assert!(chern_simons_triple(&a).abs() < 1e-14);
    }

    #[test]
    fn cubic_is_antisymmetric_under_transpose() {
        // CS(Aᵀ) = −CS(A).
        let a = vec![
            vec![0.0, 0.3, 0.8, -0.1],
            vec![-0.5, 0.0, 0.2, 0.6],
            vec![0.4, -0.7, 0.0, 0.9],
            vec![0.1, -0.2, 0.3, 0.0],
        ];
        let n = a.len();
        let at: Vec<Vec<f64>> = (0..n)
            .map(|i| (0..n).map(|j| a[j][i]).collect())
            .collect();
        let cs = chern_simons_triple(&a);
        let cst = chern_simons_triple(&at);
        assert!(cs.abs() > 1e-6, "test connection should be non-degenerate");
        assert!((cst + cs).abs() < 1e-14, "CS(Aᵀ) = {cst}, −CS(A) = {}", -cs);
    }

    #[test]
    fn sl2_exp_is_unimodular() {
        // det exp X = 1 for every traceless X (the flow stays in SL(2,ℝ)).
        for &(a, b, c) in &[(0.3, 1.1, -0.4), (0.0, 0.9, 0.9), (0.5, -0.2, 0.7), (1.2, 0.0, 0.0)] {
            let m = sl2_exp(a, b, c);
            let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
            assert!((det - 1.0).abs() < 1e-12, "det = {det} for X=({a},{b},{c})");
        }
    }

    #[test]
    fn btz_holonomy_matches_cosh_closed_form() {
        let (ell, r_plus, r_minus) = (1.0, 1.7, 0.6);
        let (tr_p, tr_m) = btz_holonomy_traces(r_plus, r_minus, ell);
        let exp_p = 2.0 * (PI * (r_plus - r_minus) / ell).cosh();
        let exp_m = 2.0 * (PI * (r_plus + r_minus) / ell).cosh();
        assert!((tr_p - exp_p).abs() < 1e-10, "tr H₊ = {tr_p}, closed form {exp_p}");
        assert!((tr_m - exp_m).abs() < 1e-10, "tr H₋ = {tr_m}, closed form {exp_m}");
    }

    #[test]
    fn btz_mass_spin_round_trip_through_holonomy() {
        // (r₊, r₋) → holonomy traces → (r₊, r₋) → (M, J), graded against the
        // closed form evaluated on the original radii.
        let (ell, r_plus, r_minus) = (1.3, 2.1, 0.8);
        let (tr_p, tr_m) = btz_holonomy_traces(r_plus, r_minus, ell);
        let (rp, rm) = btz_horizons_from_traces(tr_p, tr_m, ell);
        assert!((rp - r_plus).abs() < 1e-9, "r₊ recovered = {rp}");
        assert!((rm - r_minus).abs() < 1e-9, "r₋ recovered = {rm}");

        let (m, j) = btz_mass_spin(rp, rm, ell);
        let m_exact = (r_plus * r_plus + r_minus * r_minus) / (ell * ell);
        let j_exact = 2.0 * r_plus * r_minus / ell;
        assert!((m - m_exact).abs() < 1e-9, "M = {m}, closed form {m_exact}");
        assert!((j - j_exact).abs() < 1e-9, "J = {j}, closed form {j_exact}");

        // Cosmic-censorship / extremality sanity: |J| ≤ M ℓ.
        assert!(j.abs() <= m * ell + 1e-12, "|J| = {} exceeds Mℓ = {}", j.abs(), m * ell);
    }
}
