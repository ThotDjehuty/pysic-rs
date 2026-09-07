//! Integration tests for pysic-rs

#[cfg(test)]
mod special_functions_tests {
    use pysic::special_functions::*;

    #[test]
    fn gamma_factorial() {
        for n in 0..10 {
            let expected = (1..=n as u64).product::<u64>() as f64;
            assert!((gamma(n as f64 + 1.0) - expected).abs() < 1e-6,
                "Γ({}) = {}, expected {}", n+1, gamma(n as f64 + 1.0), expected);
        }
    }

    #[test]
    fn erf_symmetry() {
        for x in 0..20 {
            let xi = x as f64 * 0.1;
            assert!((erf(xi) + erf(-xi)).abs() < 1e-10,
                "erf({}) + erf({}) = {}", xi, -xi, erf(xi) + erf(-xi));
        }
    }

    #[test]
    fn bessel_at_zero() {
        assert!((bessel_j0(0.0) - 1.0).abs() < 1e-15);
        assert!(bessel_j1(0.0).abs() < 1e-15);
    }
}

#[cfg(test)]
mod linalg_tests {
    use pysic::linalg::*;
    use ndarray::array;

    #[test]
    fn matrix_inverse_consistency() {
        let a = array![[2.0, 1.0], [1.0, 3.0]];
        let inv = inverse(&a).unwrap();
        let product = a.dot(&inv);
        for i in 0..2 {
            for j in 0..2 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((product[[i, j]] - expected).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn cholesky_reconstruction() {
        let a = array![[4.0, 2.0], [2.0, 3.0]];
        let l = cholesky(&a).unwrap();
        let recon = l.dot(&l.t());
        assert!((recon - a).mapv(|x| x.abs()).sum() < 1e-10);
    }

    #[test]
    fn lie_bracket_antisymmetric() {
        let x_field = |_p: &[f64]| vec![1.0, 0.0, 0.0];
        let y_field = |_p: &[f64]| vec![0.0, 1.0, 0.0];
        let point = vec![0.0, 0.0, 0.0];
        let xy = lie_bracket(&x_field, &y_field, &point, 1e-5);
        let yx = lie_bracket(&y_field, &x_field, &point, 1e-5);
        for i in 0..3 {
            assert!((xy[i] + yx[i]).abs() < 1e-4,
                "[X,Y] + [Y,X] = {} (should be 0)", xy[i] + yx[i]);
        }
    }
}

#[cfg(test)]
mod ode_tests {
    use pysic::ode::rk4::*;

    #[test]
    fn rk4_harmonic_oscillator() {
        fn harmonic(t: f64, y: &[f64]) -> Vec<f64> {
            let _ = t;
            vec![y[1], -y[0]]
        }
        let result = rk4_solve(&harmonic, &[1.0, 0.0], (0.0, 10.0), 10000);
        let final_y = result.1.last().unwrap()[0];
        assert!((final_y - 10.0_f64.cos()).abs() < 0.01,
            "y(10) = {}, expected cos(10) = {}", final_y, 10.0_f64.cos());
    }
}

#[cfg(test)]
mod quantum_tests {
    use pysic::quantum::pauli::*;
    use pysic::quantum::hilbert::coherent_state;
    use num_complex::Complex64;

    #[test]
    fn pauli_x_squared() {
        let (sx, _sy, _sz) = pauli_matrices();
        let mut sx2 = [[Complex64::new(0.0, 0.0); 2]; 2];
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    sx2[i][j] += sx[i][k] * sx[k][j];
                }
            }
        }
        assert!((sx2[0][0] - 1.0).norm() < 1e-10);
        assert!(sx2[0][1].norm() < 1e-10);
        assert!(sx2[1][0].norm() < 1e-10);
        assert!((sx2[1][1] - 1.0).norm() < 1e-10);
    }

    #[test]
    fn coherent_state_normalization() {
        let alpha = Complex64::new(2.0, 1.0);
        let state = coherent_state(alpha, 50);
        let norm: f64 = state.iter().map(|c: &Complex64| c.norm_sqr()).sum();
        assert!((norm - 1.0).abs() < 1e-6,
            "Coherent state norm: {}", norm);
    }
}

#[cfg(test)]
mod gr_tests {
    use pysic::general_relativity::metrics::*;

    #[test]
    fn schwarzschild_limit() {
        let g = schwarzschild_metric(1e10, 1.0);
        assert!((g[[0, 0]] + 1.0).abs() < 1e-7);
        assert!((g[[1, 1]] - 1.0).abs() < 1e-7);
    }

    #[test]
    fn kerr_reduces_to_schwarzschild() {
        let g_kerr = kerr_metric(10.0, std::f64::consts::FRAC_PI_2, 1.0, 0.0);
        let g_schw = schwarzschild_metric(10.0, 1.0);
        for i in 0..4 {
            for j in 0..4 {
                assert!((g_kerr[[i, j]] - g_schw[[i, j]]).abs() < 1e-10,
                    "Kerr(a=0) ≠ Schwarzschild at [{},{}]: {} vs {}", i, j, g_kerr[[i,j]], g_schw[[i,j]]);
            }
        }
    }
}

#[cfg(test)]
mod casimir_tests {
    use pysic::casimir::energy::*;

    #[test]
    fn casimir_energy_scales_as_d_cubed() {
        let hbar = 1.0;
        let c = 1.0;
        let e1 = casimir_energy_parallel_plates(1.0, hbar, c);
        let e2 = casimir_energy_parallel_plates(2.0, hbar, c);
        assert!((e2 / e1 - 0.125).abs() < 1e-10,
            "E(2)/E(1) = {}, expected 0.125", e2 / e1);
    }
}
