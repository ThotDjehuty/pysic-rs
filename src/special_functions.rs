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
    // Regrouped into a single exponential: forming `t.powf(x + 0.5)` on its own
    // overflows near x = 143 even though Gamma(x) itself is finite out to
    // x ~ 171.6, so the naive product reports `inf` for arguments the type can
    // still represent.
    (2.0 * std::f64::consts::PI).sqrt() * y * ((x + 0.5) * t.ln() - t).exp()
}

/// Log of the gamma function (numerically stable for large arguments).
pub fn ln_gamma(x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NAN;
    }
    // Reflection for small arguments: ln Γ(x) = ln(π/sin πx) − ln Γ(1−x).
    if x < 0.5 {
        let s = (std::f64::consts::PI * x).sin();
        return (std::f64::consts::PI / s.abs()).ln() - ln_gamma(1.0 - x);
    }
    // Lanczos in logarithmic form — never forms Γ(x) itself, so it neither
    // overflows past x = 171 nor loses the mantissa to a large exponential.
    let z = x - 1.0;
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
        y += ci / (z + i as f64);
    }
    let t = z + 7.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (z + 0.5) * t.ln() - t + y.ln()
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

/// Euler--Mascheroni constant.
const EULER_GAMMA: f64 = 0.577_215_664_901_532_9;

/// Crossover between the ascending series and the Hankel asymptotic
/// expansion for the Bessel functions, chosen by measuring both branches
/// against tabulated values: the series holds near machine precision to
/// |x| ~ 14 and then loses digits to cancellation, while the asymptotic
/// series improves with |x| and is already ~1e-12 here, so 15 is close to
/// where the two error curves cross. The *previous* implementation switched
/// at $|x|=8$ with a four-term Taylor polynomial below it, which is the
/// `bessel` defect reported in the paper.
const BESSEL_SERIES_LIMIT: f64 = 15.0;

/// Ascending series for `J_nu(x)`, `nu` a small non-negative integer:
/// `J_nu(x) = (x/2)^nu / nu! * sum_k (-1)^k (x^2/4)^k / (k! (k+nu)!/nu!)`.
fn bessel_j_series(nu: u32, x: f64) -> f64 {
    let h = 0.25 * x * x;
    let mut term = 1.0f64;
    let mut sum = 1.0f64;
    for k in 1..200 {
        term *= -h / (k as f64 * (k as f64 + nu as f64));
        sum += term;
        if term.abs() <= 1e-19 * sum.abs().max(1e-300) {
            break;
        }
    }
    let mut pref = 1.0f64;
    for j in 1..=nu {
        pref *= (x * 0.5) / j as f64;
    }
    pref * sum
}

/// The two Hankel asymptotic amplitude series `(P, Q)` shared by `J_nu` and
/// `Y_nu`, with `a_k(nu) = prod_{j=1..k}(4nu^2-(2j-1)^2) / (k! 8^k)`:
/// `P = sum_k (-1)^k a_{2k}/x^{2k}`, `Q = sum_k (-1)^k a_{2k+1}/x^{2k+1}`.
fn bessel_hankel_pq(nu: f64, ax: f64) -> (f64, f64) {
    let mu = 4.0 * nu * nu;
    let mut a = 1.0f64;
    let mut p = 1.0f64;
    let mut q = 0.0f64;
    for k in 1..=12usize {
        a *= (mu - ((2 * k - 1) as f64).powi(2)) / (k as f64 * 8.0 * ax);
        let sgn = if (k / 2) % 2 == 0 { 1.0 } else { -1.0 };
        if k % 2 == 0 {
            p += sgn * a;
        } else {
            q += sgn * a;
        }
    }
    (p, q)
}

/// `J_nu` by the Hankel asymptotic expansion (large `|x|`), for `nu >= 0`.
fn bessel_j_asymptotic(nu: f64, ax: f64) -> f64 {
    let chi = ax - (0.5 * nu + 0.25) * std::f64::consts::PI;
    let (p, q) = bessel_hankel_pq(nu, ax);
    (2.0 / (std::f64::consts::PI * ax)).sqrt() * (chi.cos() * p - chi.sin() * q)
}

/// `Y_nu` by the Hankel asymptotic expansion (large `x > 0`).
fn bessel_y_asymptotic(nu: f64, x: f64) -> f64 {
    let chi = x - (0.5 * nu + 0.25) * std::f64::consts::PI;
    let (p, q) = bessel_hankel_pq(nu, x);
    (2.0 / (std::f64::consts::PI * x)).sqrt() * (chi.sin() * p + chi.cos() * q)
}

/// Bessel function of the first kind, order 0: J_0(x). Even in `x`.
pub fn bessel_j0(x: f64) -> f64 {
    let ax = x.abs();
    if ax <= BESSEL_SERIES_LIMIT {
        bessel_j_series(0, ax)
    } else {
        bessel_j_asymptotic(0.0, ax)
    }
}

/// Bessel function of the first kind, order 1: J_1(x). Odd in `x`.
pub fn bessel_j1(x: f64) -> f64 {
    let ax = x.abs();
    let v = if ax <= BESSEL_SERIES_LIMIT {
        bessel_j_series(1, ax)
    } else {
        bessel_j_asymptotic(1.0, ax)
    };
    if x < 0.0 {
        -v
    } else {
        v
    }
}

/// Bessel function of the second kind, order 0: Y_0(x) for x > 0.
///
/// Small-argument form (A&S 9.1.13):
/// `Y_0 = (2/pi)[(ln(x/2)+gamma) J_0(x) + sum_{k>=1} (-1)^{k+1} H_k (x^2/4)^k/(k!)^2]`.
pub fn bessel_y0(x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NAN;
    }
    if x > BESSEL_SERIES_LIMIT {
        return bessel_y_asymptotic(0.0, x);
    }
    let h = 0.25 * x * x;
    let mut term = 1.0f64; // (x^2/4)^k/(k!)^2 at k = 0
    let mut harmonic = 0.0f64;
    let mut sum = 0.0f64;
    for k in 1..200 {
        term *= -h / (k as f64 * k as f64); // carries the (-1)^k
        harmonic += 1.0 / k as f64;
        // (-1)^{k+1} H_k (x^2/4)^k/(k!)^2 = -term * H_k
        let contrib = -term * harmonic;
        sum += contrib;
        if contrib.abs() <= 1e-19 * sum.abs().max(1e-300) {
            break;
        }
    }
    (2.0 / std::f64::consts::PI) * ((x * 0.5).ln() + EULER_GAMMA) * bessel_j0(x)
        + (2.0 / std::f64::consts::PI) * sum
}

/// Bessel function of the second kind, order 1: Y_1(x) for x > 0.
///
/// Small-argument form reduced from A&S 9.1.11 at `n = 1`:
/// `Y_1 = (2/pi)(ln(x/2)+gamma) J_1 - 2/(pi x)
///        - (x/(2 pi)) sum_{k>=0} (-1)^k (H_k + H_{k+1}) (x^2/4)^k/(k!(k+1)!)`.
pub fn bessel_y1(x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NAN;
    }
    if x > BESSEL_SERIES_LIMIT {
        return bessel_y_asymptotic(1.0, x);
    }
    let h = 0.25 * x * x;
    let mut term = 1.0f64; // (x^2/4)^k/(k!(k+1)!) at k = 0
    let mut h_k = 0.0f64; // H_0
    let mut sum = (h_k + 1.0) * term; // k = 0: H_0 + H_1 = 1
    for k in 1..200 {
        term *= -h / (k as f64 * (k as f64 + 1.0)); // carries the (-1)^k
        h_k += 1.0 / k as f64; // H_k
        let h_k1 = h_k + 1.0 / (k as f64 + 1.0); // H_{k+1}
        let contrib = term * (h_k + h_k1);
        sum += contrib;
        if contrib.abs() <= 1e-19 * sum.abs().max(1e-300) {
            break;
        }
    }
    (2.0 / std::f64::consts::PI) * ((x * 0.5).ln() + EULER_GAMMA) * bessel_j1(x)
        - 2.0 / (std::f64::consts::PI * x)
        - (x / (2.0 * std::f64::consts::PI)) * sum
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
/// `Ai(0) = 3^{-2/3}/Gamma(2/3)`.
const AI_0: f64 = 0.355_028_053_887_817_24;
/// `Ai'(0) = -3^{-1/3}/Gamma(1/3)`.
const AI_PRIME_0: f64 = -0.258_819_403_792_806_8;
/// Below this `|x|` the Maclaurin series is used; beyond it, the asymptotics.
/// The old implementation applied the large-`|x|` asymptotic for *all*
/// `x >= -1`, including `x = 0`, where it divides by `xi = 0` and returns NaN
/// — the `airy` defect reported in the paper.
const AIRY_SERIES_LIMIT: f64 = 5.0;
/// On the oscillatory side Ai does not decay, so the series suffers far less
/// cancellation and can be carried further out; pushing the join to `-9` also
/// moves it to where `zeta = (2/3)|x|^{3/2}` is large enough for the
/// asymptotic series to be accurate.
const AIRY_SERIES_LIMIT_NEG: f64 = -9.0;

/// Asymptotic `u_k` coefficients: `u_0 = 1`,
/// `u_k = (6k-5)(6k-3)(6k-1) / ((2k-1) 216 k) * u_{k-1}`.
fn airy_u(k: usize) -> f64 {
    let mut u = 1.0f64;
    for j in 1..=k {
        let jf = j as f64;
        u *= ((6.0 * jf - 5.0) * (6.0 * jf - 3.0) * (6.0 * jf - 1.0))
            / ((2.0 * jf - 1.0) * 216.0 * jf);
    }
    u
}

/// Airy function of the first kind, Ai(x).
///
/// Maclaurin series for `|x| <= 5` from `y'' = x y`, whose coefficients obey
/// `a_{n+3} = a_n / ((n+3)(n+2))` seeded by `a_0 = Ai(0)`, `a_1 = Ai'(0)`,
/// `a_2 = 0`; the decaying asymptotic for `x > 5`; and the oscillatory
/// asymptotic (A&S 10.4.60) for `x < -5`.
pub fn airy_ai(x: f64) -> f64 {
    if x <= AIRY_SERIES_LIMIT && x >= AIRY_SERIES_LIMIT_NEG {
        // Three interleaved chains, one per residue of n mod 3.
        let mut a = [AI_0, AI_PRIME_0, 0.0];
        let mut xp = [1.0, x, x * x];
        let mut sum = a[0] * xp[0] + a[1] * xp[1] + a[2] * xp[2];
        let x3 = x * x * x;
        for block in 1..60usize {
            let mut delta = 0.0;
            for j in 0..3 {
                let n_prev = 3 * (block - 1) + j;
                a[j] /= ((n_prev + 3) as f64) * ((n_prev + 2) as f64);
                xp[j] *= x3;
                delta += a[j] * xp[j];
            }
            sum += delta;
            if delta.abs() <= 1e-19 * sum.abs().max(1e-300) {
                break;
            }
        }
        sum
    } else if x > 0.0 {
        // Ai(x) ~ e^{-zeta}/(2 sqrt(pi) x^{1/4}) sum_k (-1)^k u_k / zeta^k
        let zeta = 2.0 / 3.0 * x.powf(1.5);
        let mut series = 0.0f64;
        for k in 0..12usize {
            let t = airy_u(k) / zeta.powi(k as i32);
            series += if k % 2 == 0 { t } else { -t };
        }
        (-zeta).exp() / (2.0 * std::f64::consts::PI.sqrt() * x.powf(0.25)) * series
    } else {
        // DLMF 9.7.9:
        // Ai(-z) ~ pi^{-1/2} z^{-1/4}[ cos(zeta-pi/4) sum (-1)^k u_{2k}/zeta^{2k}
        //                            + sin(zeta-pi/4) sum (-1)^k u_{2k+1}/zeta^{2k+1} ]
        let z = -x;
        let zeta = 2.0 / 3.0 * z.powf(1.5);
        let mut even = 0.0f64;
        let mut odd = 0.0f64;
        for k in 0..6usize {
            let te = airy_u(2 * k) / zeta.powi((2 * k) as i32);
            let to = airy_u(2 * k + 1) / zeta.powi((2 * k + 1) as i32);
            if k % 2 == 0 {
                even += te;
                odd += to;
            } else {
                even -= te;
                odd -= to;
            }
        }
        let phase = zeta - std::f64::consts::FRAC_PI_4;
        (phase.cos() * even + phase.sin() * odd)
            / (std::f64::consts::PI.sqrt() * z.powf(0.25))
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
    expint_en(1, x)
}

/// Generalized exponential integral E_n(x) for integer n >= 0.
/// `E_n(x) = integral_1^inf e^{-xt}/t^n dt`, for integer `n >= 0`, `x > 0`.
///
/// Modified Lentz continued fraction for `x > 1`, ascending series otherwise
/// (A&S 5.1). The previous upward recurrence evaluated `exp(-(k-1))` where
/// `exp(-x)` was intended — the `expint` defect reported in the paper — which
/// made every order wrong, `n = 1` included.
pub fn expint_en(n: usize, x: f64) -> f64 {
    const FPMIN: f64 = 1e-300;
    const EPS: f64 = 1e-16;
    const MAXIT: usize = 200;

    if x < 0.0 || (x == 0.0 && (n == 0 || n == 1)) {
        return f64::NAN;
    }
    if n == 0 {
        return (-x).exp() / x;
    }
    if x == 0.0 {
        return 1.0 / (n as f64 - 1.0);
    }
    let nf = n as f64;

    if x > 1.0 {
        // Modified Lentz algorithm for the continued fraction.
        let mut b = x + nf;
        let mut c = 1.0 / FPMIN;
        let mut d = 1.0 / b;
        let mut h = d;
        for i in 1..=MAXIT {
            let a = -(i as f64) * (nf - 1.0 + i as f64);
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
        // Ascending series.
        let mut ans = if n == 1 {
            -x.ln() - EULER_GAMMA
        } else {
            1.0 / (nf - 1.0)
        };
        let mut fact = 1.0f64;
        for i in 1..=MAXIT {
            fact *= -x / i as f64;
            let del = if i != n - 1 {
                -fact / (i as f64 - (nf - 1.0))
            } else {
                let mut psi = -EULER_GAMMA;
                for ii in 1..n {
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

// ============================================================================
// Riemann zeta function
// ============================================================================

/// Riemann zeta function ζ(s) for real s > 1 via Euler-Maclaurin.
pub fn zeta(s: f64) -> f64 {
    if (s - 1.0).abs() < 1e-14 {
        return f64::INFINITY;
    }
    if s == 0.0 {
        return -0.5;
    }
    if s < 1.0 {
        // Functional equation: zeta(s) = 2^s pi^{s-1} sin(pi s/2) Gamma(1-s) zeta(1-s).
        // Negative even integers are the trivial zeros, where sin(pi s/2) = 0.
        if s < 0.0 && (s / 2.0 - (s / 2.0).round()).abs() < 1e-14 {
            return 0.0;
        }
        return 2.0_f64.powf(s)
            * std::f64::consts::PI.powf(s - 1.0)
            * (std::f64::consts::PI * s / 2.0).sin()
            * gamma(1.0 - s)
            * zeta(1.0 - s);
    }

    // Euler--Maclaurin:
    //   zeta(s) = sum_{n<N} n^{-s} + N^{1-s}/(s-1) + N^{-s}/2
    //             + sum_{k=1..M} B_{2k}/(2k)! * (s)_{2k-1} * N^{-(s+2k-1)}
    // where (s)_{2k-1} = s(s+1)...(s+2k-2) is a *real* rising factorial. The
    // previous version built that product from `floor(s)` integers and
    // truncated the exponent with `powi(s as i32)`.
    const N: usize = 20;
    // B_{2k}/(2k)! for k = 1..7
    const B_OVER_FACT: [f64; 7] = [
        1.0 / 12.0,
        -1.0 / 720.0,
        1.0 / 30240.0,
        -1.0 / 1_209_600.0,
        1.0 / 47_900_160.0,
        -691.0 / 1_307_674_368_000.0,
        1.0 / 74_724_249_600.0,
    ];

    let nf = N as f64;
    let mut sum = 0.0f64;
    for n in 1..N {
        sum += (n as f64).powf(-s);
    }
    sum += nf.powf(1.0 - s) / (s - 1.0);
    sum += 0.5 * nf.powf(-s);

    for (idx, &coef) in B_OVER_FACT.iter().enumerate() {
        let k = idx + 1;
        // (s)_{2k-1} = prod_{j=0}^{2k-2} (s + j)
        let mut rising = 1.0f64;
        for j in 0..(2 * k - 1) {
            rising *= s + j as f64;
        }
        sum += coef * rising * nf.powf(-(s + 2.0 * k as f64 - 1.0));
    }

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

    // ── Bessel ───────────────────────────────────────────────────────────
    //
    // The former test asserted only x = 0 and x = 1 at epsilon 1e-6 / 1e-4.
    // Both the series/asymptotic join and the asymptotic branch itself lay
    // outside the sampled points, so a routine returning J0(7.98) = +50.6
    // passed. These tests sample *across* the join and assert the bound
    // |J_nu| <= 1, which is what makes them able to fail.

    #[test]
    fn test_bessel_tabulated() {
        // A&S Table 9.1 reference values.
        assert_relative_eq!(bessel_j0(0.0), 1.0, epsilon = 1e-15);
        assert_relative_eq!(bessel_j0(1.0), 0.765_197_686_557_966_6, epsilon = 1e-13);
        assert_relative_eq!(bessel_j0(5.0), -0.177_596_771_314_338_3, epsilon = 1e-12);
        assert_relative_eq!(bessel_j0(10.0), -0.245_935_764_451_348_33, epsilon = 1e-10);
        assert_relative_eq!(bessel_j0(20.0), 0.167_024_664_340_583_26, epsilon = 1e-9);
        assert_relative_eq!(bessel_j1(1.0), 0.440_050_585_744_933_5, epsilon = 1e-13);
        assert_relative_eq!(bessel_j1(10.0), 0.043_472_746_168_861_44, epsilon = 1e-10);
        assert_relative_eq!(bessel_y0(1.0), 0.088_256_964_215_676_96, epsilon = 1e-12);
        assert_relative_eq!(bessel_y0(10.0), 0.055_671_167_283_599_39, epsilon = 1e-9);
        assert_relative_eq!(bessel_y1(1.0), -0.781_212_821_300_288_7, epsilon = 1e-12);
    }

    #[test]
    fn test_bessel_j_is_bounded_by_one() {
        // |J_0|, |J_1| <= 1 for all real x. This is the assertion that the
        // defective asymptotic branch violated, by a factor of fifty.
        let mut x = -40.0;
        while x <= 40.0 {
            assert!(
                bessel_j0(x).abs() <= 1.0 + 1e-12,
                "|J0({x})| = {} exceeds 1",
                bessel_j0(x).abs()
            );
            assert!(
                bessel_j1(x).abs() <= 1.0 + 1e-12,
                "|J1({x})| = {} exceeds 1",
                bessel_j1(x).abs()
            );
            x += 0.01;
        }
    }

    #[test]
    fn test_bessel_branches_agree_at_the_join() {
        // A branch switch must not introduce a step. Comparing f(j-eps) to
        // f(j+eps) would only measure the function's own slope, so instead
        // evaluate *both* branches at the same argument and require them to
        // agree: that is the statement a mangled branch actually violates.
        let j = BESSEL_SERIES_LIMIT;
        for nu in [0u32, 1u32] {
            let s = bessel_j_series(nu, j);
            let a = bessel_j_asymptotic(nu as f64, j);
            assert!(
                (s - a).abs() < 1e-10,
                "J{nu} branches disagree by {} at the join |x| = {j}",
                (s - a).abs()
            );
        }
    }

    #[test]
    fn test_bessel_wronskian() {
        // J_1(x) Y_0(x) - J_0(x) Y_1(x) = 2/(pi x) exactly, for every x > 0.
        // An identity relating all four routines, independent of any table.
        for &x in &[0.5, 1.0, 2.5, 7.98, 12.0, 17.4, 17.6, 25.0, 33.0] {
            let w = bessel_j1(x) * bessel_y0(x) - bessel_j0(x) * bessel_y1(x);
            let exact = 2.0 / (std::f64::consts::PI * x);
            assert_relative_eq!(w, exact, epsilon = 1e-9);
        }
    }

    #[test]
    fn test_airy_tabulated_and_finite() {
        // Ai(0) = 3^{-2/3}/Gamma(2/3); the old asymptotic-for-all-x returned NaN here.
        assert_relative_eq!(airy_ai(0.0), 0.355_028_053_887_817_24, epsilon = 1e-15);
        assert_relative_eq!(airy_ai(1.0), 0.135_292_416_312_881_4, epsilon = 1e-13);
        assert_relative_eq!(airy_ai(5.0), 1.083_444_281_401_7e-4, epsilon = 1e-9);
        assert_relative_eq!(airy_ai(-1.0), 0.535_560_883_292_352_2, epsilon = 1e-13);
        assert_relative_eq!(airy_ai(-15.0), 0.278_217_490_878_740_6, epsilon = 1e-8);
        let mut x = -20.0;
        while x <= 10.0 {
            assert!(airy_ai(x).is_finite(), "Ai({x}) is not finite");
            x += 0.01;
        }
    }

    #[test]
    fn test_expint_tabulated_and_recurrence() {
        // E_1(1) = 0.2193839344; the old code was wrong at every order.
        assert_relative_eq!(expint_en(1, 1.0), 0.219_383_934_395_520_3, epsilon = 1e-11);
        assert_relative_eq!(expint_en(2, 1.0), 0.148_495_506_775_922_0, epsilon = 1e-11);
        assert_relative_eq!(expint_en(3, 1.0), 0.109_691_967_197_270_9, epsilon = 1e-11);
        // Recurrence n E_{n+1}(x) = e^{-x} - x E_n(x) must hold identically.
        for &x in &[0.2, 0.9, 1.5, 4.0] {
            for n in 1..6usize {
                let lhs = n as f64 * expint_en(n + 1, x);
                let rhs = (-x).exp() - x * expint_en(n, x);
                assert_relative_eq!(lhs, rhs, epsilon = 1e-11);
            }
        }
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
        // The former tolerance was 1e-3 — three orders looser than the closed
        // forms justify, and loose enough to pass an implementation wrong in
        // the fourth digit. These are the exact values.
        assert_relative_eq!(
            zeta(2.0),
            std::f64::consts::PI * std::f64::consts::PI / 6.0,
            epsilon = 1e-14
        );
        assert_relative_eq!(zeta(4.0), std::f64::consts::PI.powi(4) / 90.0, epsilon = 1e-14);
        assert_relative_eq!(zeta(6.0), std::f64::consts::PI.powi(6) / 945.0, epsilon = 1e-14);
        // Negative-argument closed forms, reached through the functional
        // equation: zeta(-1) = -1/12, zeta(-3) = 1/120, and the trivial zeros.
        assert_relative_eq!(zeta(-1.0), -1.0 / 12.0, epsilon = 1e-12);
        assert_relative_eq!(zeta(-3.0), 1.0 / 120.0, epsilon = 1e-12);
        assert_relative_eq!(zeta(0.0), -0.5, epsilon = 1e-14);
        for &s in &[-2.0, -4.0, -6.0] {
            assert!(zeta(s).abs() < 1e-14, "zeta({s}) should be a trivial zero");
        }
    }

    #[test]
    fn test_gamma_does_not_overflow_prematurely() {
        // Gamma(x) is finite out to x ~ 171.6; forming t^(x+0.5) separately
        // overflowed near x = 143 and reported inf for representable values.
        for &x in &[100.0, 140.0, 150.0, 170.0] {
            assert!(gamma(x).is_finite(), "gamma({x}) overflowed early");
        }
        // Gamma(n+1) = n! exactly for moderate n.
        assert_relative_eq!(gamma(21.0), 2_432_902_008_176_640_000.0, max_relative = 1e-13);
        // ln Gamma must agree with ln(Gamma) where the latter is representable.
        for &x in &[0.7, 2.5, 9.0, 60.0] {
            assert_relative_eq!(ln_gamma(x), gamma(x).ln(), epsilon = 1e-12);
        }
    }
}
