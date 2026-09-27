//! Closed-form tests for the topology and electromagnetism modules.
//!
//! Same discipline as tests/gauge_correctness.rs: assert the analytic result or
//! the structural invariant, not whatever the implementation happens to return.

const PI: f64 = std::f64::consts::PI;
const C: f64 = 2.997_924_58e8;
const EPS0: f64 = 8.854_187_817e-12;
const MU0: f64 = 4e-7 * PI;

#[cfg(test)]
mod em_radiation {
    use super::*;
    use pysic::em::radiation::{
        compton_wavelength_shift, dipole_angular_distribution, dipole_radiation, larmor_formula,
        magnetic_dipole_radiation, thomson_cross_section,
    };

    /// Larmor: P = q²a²/(6π ε₀ c³).
    #[test]
    fn larmor_matches_closed_form() {
        let (q, a) = (1.602_176_634e-19, 1e20);
        let expected = q * q * a * a / (6.0 * PI * EPS0 * C.powi(3));
        let got = larmor_formula(q, a, EPS0, C);
        assert!(
            (got - expected).abs() / expected < 1e-12,
            "Larmor: got {got:e}, expected {expected:e}"
        );
    }

    /// Electric dipole: P = μ₀|p|²ω⁴/(12π c). One power of c, not three.
    #[test]
    fn electric_dipole_power_matches_closed_form() {
        let (p, omega) = (1e-30f64, 1e15f64);
        let expected = MU0 * p * p * omega.powi(4) / (12.0 * PI * C);
        let got = dipole_radiation(p, omega, MU0, C);
        assert!(
            (got - expected).abs() / expected < 1e-12,
            "dipole power: got {got:e}, expected {expected:e}"
        );
    }

    /// The same result written with ε₀: P = |p|²ω⁴/(12π ε₀ c³).
    #[test]
    fn electric_dipole_power_agrees_in_both_unit_systems() {
        let (p, omega) = (2.5e-30f64, 3e14f64);
        let via_mu0 = dipole_radiation(p, omega, MU0, C);
        let via_eps0 = p * p * omega.powi(4) / (12.0 * PI * EPS0 * C.powi(3));
        assert!(
            (via_mu0 - via_eps0).abs() / via_eps0 < 1e-9,
            "μ₀ form {via_mu0:e} disagrees with ε₀ form {via_eps0:e}"
        );
    }

    /// ∫ (dP/dΩ) dΩ must reproduce the total power, since ∫sin²θ dΩ = 8π/3.
    #[test]
    fn angular_distribution_integrates_to_total_power() {
        let (p, omega) = (1e-30, 1e15);
        let n = 20_000;
        let mut integral = 0.0;
        for k in 0..n {
            let theta = PI * (k as f64 + 0.5) / n as f64;
            let dtheta = PI / n as f64;
            // dΩ = 2π sinθ dθ for an azimuthally symmetric distribution.
            integral += dipole_angular_distribution(theta, p, omega, MU0, C)
                * 2.0
                * PI
                * theta.sin()
                * dtheta;
        }
        let total = dipole_radiation(p, omega, MU0, C);
        assert!(
            (integral - total).abs() / total < 1e-6,
            "∫dP/dΩ = {integral:e} but total power = {total:e}"
        );
    }

    /// A magnetic moment has different dimensions: P = μ₀m²ω⁴/(12π c³).
    /// Guards against "fixing" this one to match the electric case.
    #[test]
    fn magnetic_dipole_keeps_c_cubed() {
        let (m, omega) = (1e-24f64, 1e12f64);
        let expected = MU0 * m * m * omega.powi(4) / (12.0 * PI * C.powi(3));
        let got = magnetic_dipole_radiation(m, omega, MU0, C);
        assert!(
            (got - expected).abs() / expected < 1e-12,
            "magnetic dipole: got {got:e}, expected {expected:e}"
        );
    }

    /// Δλ = (h/m_e c)(1 − cos θ); at θ = π this is twice the Compton wavelength.
    #[test]
    fn compton_shift_matches_closed_form() {
        let (h, m_e) = (6.626_070_15e-34, 9.109_383_7015e-31);
        let lambda_c = h / (m_e * C);
        for theta in [0.0, PI / 3.0, PI / 2.0, PI] {
            let expected = lambda_c * (1.0 - theta.cos());
            let got = compton_wavelength_shift(theta, h, m_e, C);
            assert!(
                (got - expected).abs() <= 1e-18 + 1e-12 * expected.abs(),
                "Δλ({theta}): got {got:e}, expected {expected:e}"
            );
        }
    }

    #[test]
    fn thomson_cross_section_matches_closed_form() {
        let r_e = 2.817_940_3262e-15;
        let expected = (8.0 / 3.0) * PI * r_e * r_e;
        assert!((thomson_cross_section(r_e) - expected).abs() / expected < 1e-12);
    }
}

#[cfg(test)]
mod em_green_functions {
    use super::*;
    use pysic::em::green_fn::green_fn_static;

    /// Convention ∇²G = −δ³(r) gives G = +1/(4π|r|) — the standard EM choice,
    /// so a positive charge produces a positive potential. This must agree with
    /// the boxed equation on the docs page and with `green_fn_retarded`.
    #[test]
    fn static_green_fn_follows_documented_convention() {
        for r in [0.25, 1.0, 3.5] {
            let expected = 1.0 / (4.0 * PI * r);
            let got = green_fn_static(r);
            assert!(
                (got - expected).abs() / expected.abs() < 1e-12,
                "G({r}) = {got}, expected {expected} (convention ∇²G = −δ³)"
            );
        }
    }

    /// The static and retarded Green's functions must carry the SAME sign
    /// convention. This is the inconsistency that flagged the original bug.
    #[test]
    fn static_and_retarded_green_fns_share_a_convention() {
        use pysic::em::green_fn::green_fn_retarded;
        let r = 2.0;
        // Sample the retarded G at the light-cone peak t = r/c, where its
        // regularized delta is positive, and compare signs with the static one.
        let g_ret = green_fn_retarded(r / C, r, C);
        let g_static = green_fn_static(r);
        assert!(
            g_ret > 0.0 && g_static > 0.0,
            "sign mismatch: G_ret = {g_ret:e}, G_static = {g_static:e} — both should be positive"
        );
    }

    /// 1/r falloff.
    #[test]
    fn static_green_fn_scales_as_inverse_r() {
        let ratio = green_fn_static(1.0) / green_fn_static(4.0);
        assert!((ratio - 4.0).abs() < 1e-9, "expected 1/r falloff, got ratio {ratio}");
    }

    #[test]
    fn static_green_fn_is_regularized_at_origin() {
        assert_eq!(green_fn_static(0.0), 0.0);
    }
}

#[cfg(test)]
mod topology_invariants {
    use super::*;
    use pysic::topology::berry::berry_phase_bloch;
    use pysic::topology::chern::chern_number;
    use pysic::topology::winding::{skyrmion_number, winding_number_real};

    /// Berry phase on a cone of half-angle θ is −½ of the enclosed solid angle.
    #[test]
    fn bloch_berry_phase_is_half_the_solid_angle() {
        for theta in [0.0, PI / 6.0, PI / 3.0, PI / 2.0, PI] {
            let solid_angle = 2.0 * PI * (1.0 - theta.cos());
            let expected = -0.5 * solid_angle;
            let got = berry_phase_bloch(theta);
            assert!(
                (got - expected).abs() < 1e-12,
                "γ({theta}) = {got}, expected −Ω/2 = {expected}"
            );
        }
    }

    /// A curve traversed n times has winding number exactly n.
    #[test]
    fn winding_number_counts_traversals() {
        for n in [-3i32, -1, 1, 2, 5] {
            let steps = 4000;
            let f1: Vec<f64> = (0..steps)
                .map(|k| (n as f64 * 2.0 * PI * k as f64 / steps as f64).cos())
                .collect();
            let f2: Vec<f64> = (0..steps)
                .map(|k| (n as f64 * 2.0 * PI * k as f64 / steps as f64).sin())
                .collect();
            let got = winding_number_real(&f1, &f2);
            assert!(
                (got - n as f64).abs() < 1e-9,
                "winding number for n={n} came out as {got}"
            );
        }
    }

    /// Constant Berry curvature Ω = 2π/L² over an L × L patch gives c₁ = 1.
    #[test]
    fn chern_number_of_constant_curvature() {
        let (n, l) = (64usize, 1.0f64);
        let dk = l / n as f64;
        let curvature = vec![vec![2.0 * PI / (l * l); n]; n];
        let got = chern_number(&curvature, dk);
        assert!(
            (got - 1.0).abs() < 1e-6,
            "c₁ = {got}, expected 1 for constant curvature 2π/L²"
        );
    }

    /// Belavin–Polyakov profile: n̂ = (2x, 2y, r²−1)/(r²+1) has |Q| = 1.
    /// This is the test the previous placeholder binding could never pass.
    ///
    /// The charge density is q(r) = 4/(1+r²)², which peaks at the origin and
    /// leaves a truncation tail of exactly 1/(1+R²) outside radius R. So the
    /// grid has to resolve the core *and* reach far enough out; we assert the
    /// error shrinks as the grid improves rather than trusting one grid.
    fn belavin_polyakov_charge(span: f64, n: usize) -> f64 {
        let coords: Vec<f64> = (0..n)
            .map(|k| -span + 2.0 * span * k as f64 / (n - 1) as f64)
            .collect();
        let h = coords[1] - coords[0];
        let field: Vec<Vec<[f64; 3]>> = coords
            .iter()
            .map(|&x| {
                coords
                    .iter()
                    .map(|&y| {
                        let r2 = x * x + y * y;
                        let den = r2 + 1.0;
                        [2.0 * x / den, 2.0 * y / den, (r2 - 1.0) / den]
                    })
                    .collect()
            })
            .collect();
        skyrmion_number(&field, h, h)
    }

    #[test]
    fn skyrmion_number_of_belavin_polyakov_hedgehog() {
        // (span, n, h): each row resolves the core better and reaches further.
        let coarse = belavin_polyakov_charge(8.0, 129).abs(); // h = 0.125
        let medium = belavin_polyakov_charge(12.0, 385).abs(); // h ≈ 0.0625
        let fine = belavin_polyakov_charge(16.0, 1025).abs(); // h ≈ 0.03125

        let e_coarse = (coarse - 1.0).abs();
        let e_medium = (medium - 1.0).abs();
        let e_fine = (fine - 1.0).abs();

        println!("|Q| coarse={coarse:.6} medium={medium:.6} fine={fine:.6}");
        println!("err  coarse={e_coarse:.6} medium={e_medium:.6} fine={e_fine:.6}");

        assert!(
            e_medium < e_coarse && e_fine < e_medium,
            "error must shrink with resolution: {e_coarse:.6} -> {e_medium:.6} -> {e_fine:.6}"
        );
        assert!(
            e_fine < 0.02,
            "|Q| = {fine:.6} on the finest grid, expected 1 (error {e_fine:.6})"
        );
    }

    /// A uniform field is topologically trivial: Q = 0.
    #[test]
    fn skyrmion_number_of_uniform_field_vanishes() {
        let field: Vec<Vec<[f64; 3]>> = (0..24).map(|_| (0..24).map(|_| [0.0, 0.0, 1.0]).collect()).collect();
        let q = skyrmion_number(&field, 0.1, 0.1);
        assert!(q.abs() < 1e-12, "uniform field gave Q = {q}, expected 0");
    }
}
