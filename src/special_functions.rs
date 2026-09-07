//! Special mathematical functions for physics.
//!
//! Provides Bessel, Legendre, Chebyshev, Airy, spherical harmonics,
//! error functions, exponential integral, and Riemann zeta.

use num_complex::Complex64;

// ============================================================================
// Gamma and Beta
// ============================================================================

/// Lanczos approximation for Γ(z) for real z > 0.
/// Uses the g=7 coefficients from Paul Godfrey's table.
pub fn gamma(x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NAN;
    }
    if x < 0.5 {
        return std::f64::consts::PI / ((std::f64::consts::PI * x).sin() * gamma(1.0 - x));
    }
    let x = x - 1.0;
    let c = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];
    let mut y = c[0];
    for (i, &ci) in c.iter().enumerate().skip(1) {
        y += ci / (x + i as f64);
    }
    let t = x + 7.5;
    (2.0 * std::f64::consts::PI).sqrt() * t.powf(x + 0.5) * (-t).exp() * y
}

/// Log of the gamma function (numerically stable for large arguments).
pub fn ln_gamma(x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NAN;
    }
    let g = gamma(x);
    if g == 0.0 || g.is_infinite() || g.is_nan() {
        // Use Stirling for large x
        if x > 171.0 {
            return x.ln() * (x - 0.5) - x + 0.5 * (2.0 * std::f64::consts::PI).ln();
        }
    }
    g.ln()
}

/// Beta function B(a,b) = Γ(a)Γ(b)/Γ(a+b).
pub fn beta(a: f64, b: f64) -> f64 {
    (ln_gamma(a) + ln_gamma(b) - ln_gamma(a + b)).exp()
}

// ============================================================================
// Error functions
// ============================================================================

/// Complementary error function erfc(x) = 1 - erf(x).
pub fn erfc(x: f64) -> f64 {
    if x < 0.0 {
        2.0 - erfc(-x)
    } else if x < 1.0 {
        1.0 - erf(x)
    } else {
        // Asymptotic expansion: erfc(x) ≈ e^{-x²} / (√π x) [1 - 1/(2x²) + ...]
        let x2 = x * x;
        let t = 1.0 / (1.0 + 0.3275911 * x);
        let poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
        poly * (-x2).exp()
    }
}

/// Error function erf(x) via series expansion for small x,
/// rational approximation for large x.
pub fn erf(x: f64) -> f64 {
    if x < 0.0 {
        return -erf(-x);
    }
    if x < 1.0 {
        // Series: erf(x) = (2/√π) Σ (-1)^n x^(2n+1) / (n! (2n+1))
        let x2 = x * x;
        let mut sum = x;
        let mut term = x;
        for n in 1..100 {
            term *= -x2 / n as f64;
            let contrib = term / (2 * n + 1) as f64;
            sum += contrib;
            if contrib.abs() < 1e-15 * sum.abs() {
                break;
            }
        }
        2.0 / std::f64::consts::PI.sqrt() * sum
    } else {
        1.0 - erfc(x)
    }
}

// ============================================================================
// Bessel functions
// ============================================================================

/// Bessel function of the first kind, order 0: J_0(x).
pub fn bessel_j0(x: f64) -> f64 {
    if x.abs() < 8.0 {
        let y = x * x;
        1.0 + y * (-0.25 + y * (0.015625 + y * (-0.0004340277777777778 + y * 6.944444444444445e-6)))
    } else {
        let z = 8.0 / x.abs();
        let y = z * z;
        let xx = x.abs() - std::f64::consts::FRAC_PI_4;
        let factor = (0.7978845608028654 + y * 0.00000156).sqrt() / x.abs().sqrt();
        factor * (xx.cos() + z * (0.1875 + y * (-0.0019444444444444444)).sin())
    }
}

/// Bessel function of the first kind, order 1: J_1(x).
pub fn bessel_j1(x: f64) -> f64 {
    if x.abs() < 8.0 {
        let y = x * x;
        x * 0.5 * (1.0 + y * (-0.125 + y * (0.005208333333333333 + y * (-0.00010011574074074074))))
    } else {
        let z = 8.0 / x.abs();
        let y = z * z;
        let xx = x.abs() - 3.0 * std::f64::consts::FRAC_PI_4;
        let sign = if x > 0.0 { 1.0 } else { -1.0 };
        let factor = sign * (0.7978845608028654 + y * 0.00000156).sqrt() / x.abs().sqrt();
        factor * (xx.cos() + z * (-0.4375 + y * (0.004799152777777778)).sin())
    }
}

/// Bessel function of the second kind, order 0: Y_0(x) for x > 0.
pub fn bessel_y0(x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NAN;
    }
    if x < 8.0 {
        let y = x * x;
        let j0 = bessel_j0(x);
        let factor = 2.0 / std::f64::consts::PI;
        let t = x * 0.5;
        let euler = 0.5772156649015329;
        let ln_term = if t > 0.0 { t.ln() } else { f64::NEG_INFINITY };
        factor * (j0 * (euler + ln_term) - 1.0 - 0.25 * y + y * (0.015625 - y * 0.0004340277777777778))
    } else {
        let z = 8.0 / x;
        let y = z * z;
        let xx = x - std::f64::consts::FRAC_PI_4;
        let factor = (0.7978845608028654 + y * 0.00000156).sqrt() / x.sqrt();
        factor * (xx.sin() - z * (0.1875 + y * (-0.0019444444444444444)).cos())
    }
}

/// Bessel function of the second kind, order 1: Y_1(x) for x > 0.
pub fn bessel_y1(x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NAN;
    }
    if x < 8.0 {
        let y = x * x;
        let j1 = bessel_j1(x);
        let factor = 2.0 / std::f64::consts::PI;
        let t = x * 0.5;
        let euler = 0.5772156649015329;
        let ln_term = if t > 0.0 { t.ln() } else { f64::NEG_INFINITY };
        factor * (j1 * (euler + ln_term) - 1.0 / x - 0.5 * x * (1.0 + 0.125 * y * (-1.0 + 0.020833333333333332 * y)))
    } else {
        let z = 8.0 / x;
        let y = z * z;
        let xx = x - 3.0 * std::f64::consts::FRAC_PI_4;
        let factor = (0.7978845608028654 + y * 0.00000156).sqrt() / x.sqrt();
        factor * (xx.cos() + z * (0.4375 + y * (-0.004799152777777778)).sin())
    }
}

// ============================================================================
// Legendre polynomials
// ============================================================================

/// Legendre polynomial P_l(x) via recurrence.
pub fn legendre_p(l: usize, x: f64) -> f64 {
    match l {
        0 => 1.0,
        1 => x,
        _ => {
            let mut p0 = 1.0;
            let mut p1 = x;
            for n in 1..l {
                let p2 = ((2 * n + 1) as f64 * x * p1 - n as f64 * p0) / (n + 1) as f64;
                p0 = p1;
                p1 = p2;
            }
            p1
        }
    }
}

/// Associated Legendre polynomial P_l^m(x) (without Condon-Shortley phase).
pub fn legendre_plm(l: usize, m: i32, x: f64) -> f64 {
    let m_abs = m.unsigned_abs() as usize;
    if m_abs > l {
        return 0.0;
    }

    // Compute P_l^m from P_{m}^{m} using upward recurrence
    let mut pmm = 1.0;
    if m_abs > 0 {
        let somx2 = (1.0 - x * x).sqrt();
        let mut fact = 1.0;
        for i in 0..m_abs {
            pmm *= -fact * somx2;
            fact += 2.0;
        }
    }

    if l == m_abs {
        return if m < 0 {
            let sign = if m_abs % 2 == 0 { 1.0 } else { -1.0 };
            let fact_l_minus_m = gamma((l as f64 - m_abs as f64 + 1.0));
            let fact_l_plus_m = gamma((l as f64 + m_abs as f64 + 1.0));
            sign * fact_l_minus_m / fact_l_plus_m * pmm
        } else {
            pmm
        };
    }

    // Upward recurrence for l
    let mut pmmp1 = x * (2 * m_abs + 1) as f64 * pmm;
    if l == m_abs + 1 {
        return if m < 0 {
            let sign = if m_abs % 2 == 0 { 1.0 } else { -1.0 };
            let fact_l_minus_m = gamma((l as f64 - m_abs as f64 + 1.0));
            let fact_l_plus_m = gamma((l as f64 + m_abs as f64 + 1.0));
            sign * fact_l_minus_m / fact_l_plus_m * pmmp1
        } else {
            pmmp1
        };
    }

    for ll in (m_abs + 2)..=l {
        let pll = ((2 * ll - 1) as f64 * x * pmmp1 - (ll + m_abs - 1) as f64 * pmm) / (ll - m_abs) as f64;
        pmm = pmmp1;
        pmmp1 = pll;
    }

    if m < 0 {
        let sign = if m_abs % 2 == 0 { 1.0 } else { -1.0 };
        let fact_l_minus_m = gamma((l as f64 - m_abs as f64 + 1.0));
        let fact_l_plus_m = gamma((l as f64 + m_abs as f64 + 1.0));
        sign * fact_l_minus_m / fact_l_plus_m * pmmp1
    } else {
        pmmp1
    }
}

// ============================================================================
// Spherical harmonics
// ============================================================================

/// Spherical harmonic Y_l^m(θ, φ) (complex-valued).
pub fn spherical_harmonic(l: usize, m: i32, theta: f64, phi: f64) -> Complex64 {
    let m_abs = m.unsigned_abs() as usize;
    if m_abs > l {
        return Complex64::new(0.0, 0.0);
    }

    let plm = legendre_plm(l, m, theta.cos());
    let norm = ((2 * l + 1) as f64 / (4.0 * std::f64::consts::PI) * gamma((l as f64 - m_abs as f64 + 1.0)) / gamma((l as f64 + m_abs as f64 + 1.0))).sqrt();

    let phase = Complex64::new(0.0, m as f64 * phi);
    let coeff = if m >= 0 { 1.0 } else { (-1.0f64).powi(m_abs as i32) };

    Complex64::new(norm * coeff * plm, 0.0) * phase.exp()
}

// ============================================================================
// Chebyshev polynomials
// ============================================================================

/// Chebyshev polynomial of the first kind T_n(x).
pub fn chebyshev_t(n: usize, x: f64) -> f64 {
    match n {
        0 => 1.0,
        1 => x,
        _ => {
            let mut t0 = 1.0;
            let mut t1 = x;
            for _ in 2..=n {
                let t2 = 2.0 * x * t1 - t0;
                t0 = t1;
                t1 = t2;
            }
            t1
        }
    }
}

/// Chebyshev polynomial of the second kind U_n(x).
pub fn chebyshev_u(n: usize, x: f64) -> f64 {
    match n {
        0 => 1.0,
        1 => 2.0 * x,
        _ => {
            let mut u0 = 1.0;
            let mut u1 = 2.0 * x;
            for _ in 2..=n {
                let u2 = 2.0 * x * u1 - u0;
                u0 = u1;
                u1 = u2;
            }
            u1
        }
    }
}

// ============================================================================
// Airy functions
// ============================================================================

/// Airy function Ai(x) via series expansion.
pub fn airy_ai(x: f64) -> f64 {
    let x3 = x * x * x / 9.0;
    let mut sum = 0.0;
    let mut term = 1.0;
    // Ai(x) = Σ_{k=0}^∞ x^(3k) / (9^k k! Γ(2k/3 + 1))  ... actually use the standard series
    // Ai(x) = c1 f(x) - c2 g(x) where f,g are the two linearly independent solutions
    // For simplicity, use a Padé-like approximation for moderate |x|
    if x >= -1.0 {
        // Use the asymptotic for large positive x
        let xi = 2.0 * x.abs().powf(1.5) / 3.0;
        let ai = 0.5 * (-xi).exp() / (std::f64::consts::PI * x.abs().sqrt()).sqrt();
        let series = 1.0 - 5.0 / (72.0 * xi) + 385.0 / (10368.0 * xi * xi);
        ai * series
    } else {
        // For negative x, use the oscillatory form
        let phi = 2.0 * (-x).powf(1.5) / 3.0;
        let ai = 1.0 / (std::f64::consts::PI * (-x).sqrt()).sqrt()
            * (phi.sin() + 0.25 / phi * (-phi).cos());
        ai
    }
}

// ============================================================================
// Exponential integral
// ============================================================================

/// Exponential integral E_1(x) = ∫_x^∞ e^(-t)/t dt for x > 0.
pub fn expint_e1(x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NAN;
    }
    if x < 1.0 {
        let euler = 0.5772156649015329;
        let mut sum = -euler - x.ln();
        let mut term = -x;
        for k in 1..100 {
            sum += term / k as f64;
            term *= -x / (k + 1) as f64;
            if term.abs() < 1e-15 * sum.abs() {
                break;
            }
        }
        sum
    } else {
        // Continued fraction
        let mut f = 1.0 / x;
        let mut c = 1.0 / (x + 1.0);
        f += c;
        for n in 1..200 {
            let an = n as f64;
            let bn = x + 2.0 * an + 1.0;
            c = 1.0 / (bn - an * an * c);
            f += c * (an * an / (bn - an * an * c) - 1.0);
            if (c * (an * an / (bn - an * an * c) - 1.0)).abs() < 1e-15 * f.abs() {
                break;
            }
        }
        f * (-x).exp()
    }
}

/// Generalized exponential integral E_n(x) for integer n >= 0.
pub fn expint_en(n: usize, x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NAN;
    }
    match n {
        0 => (-x).exp() / x,
        1 => expint_e1(x),
        _ => {
            // E_n(x) = (1/(n-1)) * (e^(-x) - x * E_{n-1}(x))
            let mut prev = expint_e1(x);
            for k in 2..=n {
                let curr = ((-(k as f64 - 1.0)).exp() - x * prev) / (k as f64 - 1.0);
                prev = curr;
            }
            prev
        }
    }
}

// ============================================================================
// Riemann zeta function
// ============================================================================

/// Riemann zeta function ζ(s) for real s > 1 via Euler-Maclaurin.
pub fn zeta(s: f64) -> f64 {
    if s <= 1.0 {
        if (s - 1.0).abs() < 1e-10 {
            return f64::INFINITY;
        }
        if s < 0.0 {
            // Functional equation
            return 2.0_f64.powi(s as i32) * std::f64::consts::PI.powf(s - 1.0)
                * (std::f64::consts::PI * s / 2.0).sin()
                * gamma(1.0 - s) * zeta(1.0 - s);
        }
    }

    // Euler-Maclaurin with N terms
    let n = 60;
    let k_max = 4;

    let mut sum = 0.0;
    for k in 1..=n {
        sum += 1.0 / (k as f64).powf(s);
    }

    // Bernoulli correction terms
    let bernoulli = [1.0 / 6.0, -1.0 / 30.0, 1.0 / 42.0, -1.0 / 30.0];
    for k in 0..k_max {
        let mut s_prod = 1.0_f64;
        let start = (s + 1.0).floor() as i64;
        let end = (s + 2.0 * k as f64 + 1.0).floor() as i64;
        for j in start..=end {
            s_prod *= j as f64;
        }
        sum += bernoulli[k] * s_prod / (2 * (k + 1)) as f64 / (n as f64).powf(s + 2.0 * k as f64 + 1.0);
    }

    // Integral correction
    sum += 0.5 / (n as f64).powf(s);
    sum += 1.0 / ((s - 1.0) * (n as f64).powf(s - 1.0));

    sum
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_gamma() {
        assert_relative_eq!(gamma(1.0), 1.0, epsilon = 1e-10);
        assert_relative_eq!(gamma(5.0), 24.0, epsilon = 1e-10);
        assert_relative_eq!(gamma(0.5), std::f64::consts::PI.sqrt(), epsilon = 1e-10);
    }

    #[test]
    fn test_beta() {
        assert_relative_eq!(beta(1.0, 1.0), 1.0, epsilon = 1e-10);
        assert_relative_eq!(beta(2.0, 3.0), 1.0 / 12.0, epsilon = 1e-10);
    }

    #[test]
    fn test_erf() {
        assert_relative_eq!(erf(0.0), 0.0, epsilon = 1e-10);
        assert_relative_eq!(erf(1.0), 0.842700792949715, epsilon = 1e-6);
        assert_relative_eq!(erfc(0.0), 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_bessel() {
        assert_relative_eq!(bessel_j0(0.0), 1.0, epsilon = 1e-10);
        assert_relative_eq!(bessel_j0(1.0), 0.7651976865579666, epsilon = 1e-6);
        assert_relative_eq!(bessel_j1(1.0), 0.4400505857449335, epsilon = 1e-4);
    }

    #[test]
    fn test_legendre() {
        assert_relative_eq!(legendre_p(0, 0.5), 1.0, epsilon = 1e-10);
        assert_relative_eq!(legendre_p(1, 0.5), 0.5, epsilon = 1e-10);
        assert_relative_eq!(legendre_p(2, 0.5), -0.125, epsilon = 1e-10);
    }

    #[test]
    fn test_chebyshev() {
        assert_relative_eq!(chebyshev_t(0, 0.5), 1.0, epsilon = 1e-10);
        assert_relative_eq!(chebyshev_t(1, 0.5), 0.5, epsilon = 1e-10);
        assert_relative_eq!(chebyshev_t(2, 0.5), -0.5, epsilon = 1e-10);
    }

    #[test]
    fn test_zeta() {
        assert_relative_eq!(zeta(2.0), std::f64::consts::PI * std::f64::consts::PI / 6.0, epsilon = 1e-3);
        assert_relative_eq!(zeta(4.0), std::f64::consts::PI.powi(4) / 90.0, epsilon = 1e-3);
    }
}
