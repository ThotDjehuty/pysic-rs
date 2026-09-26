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
///
/// Rational minimax approximations (Abramowitz & Stegun 9.4.1/9.4.3), accurate
/// to about 1e-8 absolute. A short Taylor series is not usable here: it is only
/// good for |x| of order 1, and the oscillation must be carried by an explicit
/// phase for larger arguments.
pub fn bessel_j0(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 1e-3 {
        // The rational fit carries ~1e-9 of error at the origin; near zero the
        // Taylor series is exact to machine precision and gives J_0(0) = 1.
        let y = x * x;
        return 1.0 + y * (-0.25 + y / 64.0);
    }
    if ax < 8.0 {
        let y = x * x;
        let p1 = 57568490574.0
            + y * (-13362590354.0
                + y * (651619640.7 + y * (-11214424.18 + y * (77392.33017 + y * -184.9052456))));
        let p2 = 57568490411.0
            + y * (1029532985.0 + y * (9494680.718 + y * (59272.64853 + y * (267.8532712 + y))));
        p1 / p2
    } else {
        let z = 8.0 / ax;
        let y = z * z;
        let xx = ax - 0.785398164;
        let p1 = 1.0
            + y * (-0.1098628627e-2
                + y * (0.2734510407e-4 + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
        let p2 = -0.1562499995e-1
            + y * (0.1430488765e-3
                + y * (-0.6911147651e-5 + y * (0.7621095161e-6 + y * -0.934935152e-7)));
        (0.636619772 / ax).sqrt() * (xx.cos() * p1 - z * xx.sin() * p2)
    }
}

/// Bessel function of the first kind, order 1: J_1(x).
pub fn bessel_j1(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 1e-3 {
        let y = x * x;
        return 0.5 * x * (1.0 + y * (-0.125 + y / 192.0));
    }
    let ans = if ax < 8.0 {
        let y = x * x;
        let p1 = ax
            * (72362614232.0
                + y * (-7895059235.0
                    + y * (242396853.1
                        + y * (-2972611.439 + y * (15704.48260 + y * -30.16036606)))));
        let p2 = 144725228442.0
            + y * (2300535178.0 + y * (18583304.74 + y * (99447.43394 + y * (376.9991397 + y))));
        p1 / p2
    } else {
        let z = 8.0 / ax;
        let y = z * z;
        let xx = ax - 2.356194491;
        let p1 = 1.0
            + y * (0.183105e-2
                + y * (-0.3516396496e-4 + y * (0.2457520174e-5 + y * -0.240337019e-6)));
        let p2 = 0.04687499995
            + y * (-0.2002690873e-3
                + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)));
        (0.636619772 / ax).sqrt() * (xx.cos() * p1 - z * xx.sin() * p2)
    };
    if x < 0.0 {
        -ans
    } else {
        ans
    }
}

/// Bessel function of the second kind, order 0: Y_0(x), for x > 0.
///
/// Singular at the origin: `Y_0(x) ~ (2/pi) ln x` as `x -> 0+`. Returns NaN for
/// non-positive arguments.
pub fn bessel_y0(x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NAN;
    }
    if x < 8.0 {
        let y = x * x;
        let p1 = -2957821389.0
            + y * (7062834065.0
                + y * (-512359803.6 + y * (10879881.29 + y * (-86327.92757 + y * 228.4622733))));
        let p2 = 40076544269.0
            + y * (745249964.8 + y * (7189466.438 + y * (47447.26470 + y * (226.1030244 + y))));
        p1 / p2 + 0.636619772 * bessel_j0(x) * x.ln()
    } else {
        let z = 8.0 / x;
        let y = z * z;
        let xx = x - 0.785398164;
        let p1 = 1.0
            + y * (-0.1098628627e-2
                + y * (0.2734510407e-4 + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
        let p2 = -0.1562499995e-1
            + y * (0.1430488765e-3
                + y * (-0.6911147651e-5 + y * (0.7621095161e-6 + y * -0.934935152e-7)));
        (0.636619772 / x).sqrt() * (xx.sin() * p1 + z * xx.cos() * p2)
    }
}

/// Bessel function of the second kind, order 1: Y_1(x), for x > 0.
pub fn bessel_y1(x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NAN;
    }
    if x < 8.0 {
        let y = x * x;
        let p1 = x
            * (-4900604943000.0
                + y * (1275274390000.0
                    + y * (-51534381390.0
                        + y * (734926455.1 + y * (-4237922.726 + y * 8511.937935)))));
        let p2 = 24995805700000.0
            + y * (424441966400.0
                + y * (3733650367.0
                    + y * (22459040.02 + y * (102042.605 + y * (354.9632885 + y)))));
        p1 / p2 + 0.636619772 * (bessel_j1(x) * x.ln() - 1.0 / x)
    } else {
        let z = 8.0 / x;
        let y = z * z;
        let xx = x - 2.356194491;
        let p1 = 1.0
            + y * (0.183105e-2
                + y * (-0.3516396496e-4 + y * (0.2457520174e-5 + y * -0.240337019e-6)));
        let p2 = 0.04687499995
            + y * (-0.2002690873e-3
                + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)));
        (0.636619772 / x).sqrt() * (xx.sin() * p1 + z * xx.cos() * p2)
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
    // Maclaurin series Ai(x) = c1 f(x) - c2 g(x). Ai is entire, so the series
    // converges everywhere; past |x| ~ 7 cancellation eats the precision and we
    // switch to the asymptotic expansions instead.
    const C1: f64 = 0.355028053887817239; //  Ai(0)  = 3^(-2/3)/Gamma(2/3)
    const C2: f64 = 0.258819403792806798; // -Ai'(0) = 3^(-1/3)/Gamma(1/3)
    const U1: f64 = 5.0 / 72.0;
    const U2: f64 = 385.0 / 10368.0;
    const U3: f64 = 85085.0 / 2239488.0;

    if x.abs() <= 9.0 {
        let x3 = x * x * x;
        // f = sum_k t_k with t_0 = 1 and t_{k+1}/t_k = x^3 / ((3k+2)(3k+3))
        let mut f = 1.0;
        let mut tf = 1.0;
        // g = sum_k u_k with u_0 = x and u_{k+1}/u_k = x^3 / ((3k+3)(3k+4))
        let mut g = x;
        let mut tg = x;
        for k in 0..60 {
            let kf = k as f64;
            tf *= x3 / ((3.0 * kf + 2.0) * (3.0 * kf + 3.0));
            tg *= x3 / ((3.0 * kf + 3.0) * (3.0 * kf + 4.0));
            f += tf;
            g += tg;
            if tf.abs() < 1e-18 * f.abs().max(1e-300)
                && tg.abs() < 1e-18 * g.abs().max(1e-300)
            {
                break;
            }
        }
        C1 * f - C2 * g
    } else if x > 0.0 {
        // Ai(x) ~ exp(-xi) / (2 sqrt(pi) x^(1/4)) * sum (-1)^k u_k / xi^k
        let xi = 2.0 * x.powf(1.5) / 3.0;
        // u_k coefficients of the Airy asymptotic expansion.
        let series = 1.0 - U1 / xi + U2 / (xi * xi) - U3 / (xi * xi * xi);
        (-xi).exp() / (2.0 * std::f64::consts::PI.sqrt() * x.powf(0.25)) * series
    } else {
        // Ai(-z) ~ sin(xi + pi/4) / (sqrt(pi) z^(1/4)), oscillatory branch
        let z = -x;
        let xi = 2.0 * z.powf(1.5) / 3.0;
        let ph = xi + std::f64::consts::FRAC_PI_4;
        // Even u_k ride the sine, odd u_k the cosine.
        let sin_part = 1.0 - U2 / (xi * xi);
        let cos_part = U1 / xi - U3 / (xi * xi * xi);
        (ph.sin() * sin_part - ph.cos() * cos_part)
            / (std::f64::consts::PI.sqrt() * z.powf(0.25))
    }
}

// ============================================================================
// Exponential integral
// ============================================================================

/// Generalized exponential integral `E_n(x) = \int_1^\infty e^{-xt}/t^n dt`,
/// for `x > 0` (and `x >= 0` when `n > 1`).
///
/// Series for small `x`, modified-Lentz continued fraction otherwise
/// (Numerical Recipes 6.3). Returns NaN on arguments outside the domain.
pub fn expint_en(n: usize, x: f64) -> f64 {
    const EULER: f64 = 0.5772156649015329;
    const MAXIT: usize = 200;
    const EPS: f64 = 1e-15;
    const FPMIN: f64 = 1e-300;

    if x < 0.0 || (x == 0.0 && (n == 0 || n == 1)) {
        return f64::NAN;
    }
    if n == 0 {
        return (-x).exp() / x;
    }
    if x == 0.0 {
        return 1.0 / (n as f64 - 1.0);
    }

    let nm1 = n - 1;
    if x > 1.0 {
        // Continued fraction, evaluated by modified Lentz.
        let mut b = x + n as f64;
        let mut c = 1.0 / FPMIN;
        let mut d = 1.0 / b;
        let mut h = d;
        for i in 1..=MAXIT {
            let a = -((i * (nm1 + i)) as f64);
            b += 2.0;
            d = 1.0 / (a * d + b);
            c = b + a / c;
            let del = c * d;
            h *= del;
            if (del - 1.0).abs() < EPS {
                return h * (-x).exp();
            }
        }
        h * (-x).exp()
    } else {
        // Power series about the origin.
        let mut ans = if nm1 != 0 {
            1.0 / nm1 as f64
        } else {
            -x.ln() - EULER
        };
        let mut fact = 1.0;
        for i in 1..=MAXIT {
            fact *= -x / i as f64;
            let del = if i != nm1 {
                -fact / (i as f64 - nm1 as f64)
            } else {
                let mut psi = -EULER;
                for ii in 1..=nm1 {
                    psi += 1.0 / ii as f64;
                }
                fact * (-x.ln() + psi)
            };
            ans += del;
            if del.abs() < ans.abs() * EPS {
                return ans;
            }
        }
        ans
    }
}

/// Exponential integral `E_1(x) = \int_x^\infty e^{-t}/t\,dt` for `x > 0`.
pub fn expint_e1(x: f64) -> f64 {
    expint_en(1, x)
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

    /// Bessel functions are bounded by 1 (J) or grow only logarithmically at the
    /// origin (Y); they must also hit their tabulated zeros. Before the rewrite
    /// J_0(7.98) returned +50.6 instead of +0.176, because a four-term Taylor
    /// series was being used all the way out to |x| = 8.
    #[test]
    fn test_bessel_reference_values() {
        // Abramowitz & Stegun table 9.1
        for &(x, j0v, j1v) in &[
            (1.0_f64, 0.7651976866, 0.4400505857),
            (2.0, 0.2238907791, 0.5767248078),
            (5.0, -0.1775967713, -0.3275791376),
            (8.0, 0.1716508071, 0.2346363469),
            (10.0, -0.2459357645, 0.0434727462),
            (15.0, -0.0142244728, 0.2051040386),
        ] {
            assert!((bessel_j0(x) - j0v).abs() < 1e-7, "J0({x}) = {}", bessel_j0(x));
            assert!((bessel_j1(x) - j1v).abs() < 1e-7, "J1({x}) = {}", bessel_j1(x));
        }
        for &(x, y0v, y1v) in &[
            (1.0_f64, 0.0882569642, -0.7812128213),
            (5.0, -0.3085176252, 0.1478631434),
            (10.0, 0.0556711673, 0.2490154242),
        ] {
            assert!((bessel_y0(x) - y0v).abs() < 1e-7, "Y0({x}) = {}", bessel_y0(x));
            assert!((bessel_y1(x) - y1v).abs() < 1e-7, "Y1({x}) = {}", bessel_y1(x));
        }
        // Boundedness and parity, including across the 8.0 branch switch.
        let mut x = -25.0;
        while x <= 25.0 {
            assert!(bessel_j0(x).abs() <= 1.0 + 1e-9, "J0({x}) = {}", bessel_j0(x));
            assert!(bessel_j1(x).abs() <= 0.6, "J1({x}) = {}", bessel_j1(x));
            assert!((bessel_j0(x) - bessel_j0(-x)).abs() < 1e-12, "J0 not even at {x}");
            assert!((bessel_j1(x) + bessel_j1(-x)).abs() < 1e-12, "J1 not odd at {x}");
            x += 0.05;
        }
        // First zeros of J0 and J1.
        assert!(bessel_j0(2.404825558).abs() < 1e-7);
        assert!(bessel_j1(3.831705970).abs() < 1e-7);
    }

    /// Ai(0) and Ai'(0) are known in closed form, and Ai must stay finite
    /// everywhere. The previous implementation applied the large-|x| asymptotic
    /// at x = 0, where xi = 0, and returned NaN.
    #[test]
    fn test_airy_reference_values() {
        assert!((airy_ai(0.0) - 0.3550280538878172).abs() < 1e-12, "Ai(0) = {}", airy_ai(0.0));
        // A&S table 10.11
        for &(x, v) in &[
            (1.0_f64, 0.1352924163),
            (2.0, 0.0349241304),
            (5.0, 1.0834442e-4),
            (-1.0, 0.5355608833),
            (-2.0, 0.2274074282),
            (-5.0, 0.3507610090),
        ] {
            assert!((airy_ai(x) - v).abs() < 1e-8, "Ai({x}) = {}", airy_ai(x));
        }
        // Finite and bounded across the series/asymptotic switch at |x| = 9.
        let mut x = -15.0;
        while x <= 12.0 {
            let a = airy_ai(x);
            assert!(a.is_finite(), "Ai({x}) is not finite");
            assert!(a.abs() < 1.0, "Ai({x}) = {a} out of range");
            x += 0.05;
        }
        // First zero of Ai.
        assert!(airy_ai(-2.338107410).abs() < 1e-8);
    }

    /// E_n satisfies the recurrence n E_{n+1}(x) = e^{-x} - x E_n(x), and
    /// E_1 has known values. The old code wrote exp(-(k-1)) where exp(-x) was
    /// meant, so E_n was wrong for every n > 1.
    #[test]
    fn test_expint_reference_values() {
        // A&S table 5.1
        for &(x, v) in &[
            (0.5_f64, 0.5597735947),
            (1.0, 0.2193839344),
            (2.0, 0.0489005107),
            (5.0, 1.1482955e-3),
        ] {
            assert!((expint_e1(x) - v).abs() < 1e-9, "E1({x}) = {}", expint_e1(x));
        }
        // E_n(0) = 1/(n-1) for n > 1.
        for n in 2..8 {
            assert!(
                (expint_en(n, 0.0) - 1.0 / (n as f64 - 1.0)).abs() < 1e-12,
                "E_{n}(0) = {}",
                expint_en(n, 0.0)
            );
        }
        // Recurrence, across both the series and continued-fraction branches.
        for &x in &[0.1_f64, 0.5, 0.9, 1.1, 2.0, 5.0, 10.0] {
            for n in 1..6 {
                let lhs = n as f64 * expint_en(n + 1, x);
                let rhs = (-x).exp() - x * expint_en(n, x);
                assert!(
                    (lhs - rhs).abs() < 1e-10 * lhs.abs().max(1e-10),
                    "recurrence broken at n={n}, x={x}: {lhs} vs {rhs}"
                );
            }
        }
    }
}
