//! Green's functions for the electromagnetic field.

use num_complex::Complex64;

/// Static Green's function for the Laplacian in 3D, in the convention
/// ∇²G = -δ³(r):
/// G(r) = 1/(4π|r|)
///
/// This is the standard electromagnetic convention, so that the potential of a
/// positive charge comes out positive via φ = ∫G ρ/ε₀ d³r'. It matches
/// `green_fn_retarded` below, which likewise carries a positive 1/(4πr).
pub fn green_fn_static(r: f64) -> f64 {
    if r < 1e-15 {
        0.0
    } else {
        1.0 / (4.0 * std::f64::consts::PI * r)
    }
}

/// Retarded Green's function for the wave equation:
/// G_R(t, r) = θ(t) δ(t - r/c) / (4πr)
pub fn green_fn_retarded(t: f64, r: f64, c: f64) -> f64 {
    if t < 0.0 || r < 1e-15 {
        return 0.0;
    }
    // Approximate delta function with Gaussian
    let sigma = 0.01 * r / c;
    let dt = t - r / c;
    let delta = (-(dt * dt) / (2.0 * sigma * sigma)).exp() / (sigma * (2.0 * std::f64::consts::PI).sqrt());
    delta / (4.0 * std::f64::consts::PI * r)
}

/// Advanced Green's function: G_A(t, r) = G_R(-t, r).
pub fn green_fn_advanced(t: f64, r: f64, c: f64) -> f64 {
    green_fn_retarded(-t, r, c)
}

/// Free-space Green's function for the Helmholtz equation:
/// G(k, r) = -e^{ikr}/(4πr) (outgoing wave)
pub fn green_fn_helmholtz(k: f64, r: f64) -> Complex64 {
    if r < 1e-15 {
        return Complex64::new(0.0, 0.0);
    }
    let ikr = Complex64::new(0.0, k * r);
    -ikr.exp() / (4.0 * std::f64::consts::PI * r)
}

/// Green's function for the Poisson equation on a sphere.
pub fn green_fn_sphere(l: usize, theta: f64) -> f64 {
    // G_l(cos θ) = -1/(2l+1) P_l(cos θ)
    let p_l = crate::special_functions::legendre_p(l, theta.cos());
    -p_l / (2.0 * l as f64 + 1.0)
}

/// Dyadic Green's function for Maxwell's equations (free space):
/// G_{ij}(r) = (δ_{ij} + ∇_i∇_j/k²) e^{ikr}/(4πr)
pub fn green_fn_maxwell_dyadic(
    k: f64,
    r: f64,
    delta_ij: bool,
) -> Complex64 {
    if r < 1e-15 {
        return Complex64::new(0.0, 0.0);
    }
    let ikr = Complex64::new(0.0, k * r);
    let scalar = ikr.exp() / (4.0 * std::f64::consts::PI * r);

    if delta_ij {
        // Simplified: diagonal component
        scalar * (1.0 + Complex64::new(0.0, -1.0 / (k * r)))
    } else {
        // Off-diagonal: involves angular derivatives, simplified here
        scalar * (0.0 - 1.0 / (k * r * k * r))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_static_gfn() {
        let g = green_fn_static(1.0);
        assert!((g - 1.0 / (4.0 * std::f64::consts::PI)).abs() < 1e-10);
    }

    /// Gauss's law fixes the coefficient: for ∇²G = -δ³, the flux of ∇G through
    /// any enclosing sphere is exactly -1.
    #[test]
    fn static_gfn_satisfies_gauss_law() {
        for r in [0.3, 1.0, 7.0] {
            // G = 1/(4πr) ⇒ ∂_r G = -1/(4πr²) ⇒ flux = ∂_r G · 4πr² = -1.
            let dr = 1e-6;
            let dgdr = (green_fn_static(r + dr) - green_fn_static(r - dr)) / (2.0 * dr);
            let flux = dgdr * 4.0 * std::f64::consts::PI * r * r;
            assert!((flux + 1.0).abs() < 1e-6, "flux at r={r} was {flux}, expected -1");
        }
    }

    #[test]
    fn test_helmholtz_gfn() {
        let g = green_fn_helmholtz(1.0, 1.0);
        assert!(g.norm() > 0.0);
    }
}
