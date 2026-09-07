//! Green's functions for the electromagnetic field.

use num_complex::Complex64;

/// Static Green's function for the Laplacian in 3D:
/// G(r) = -1/(4π|r|)
pub fn green_fn_static(r: f64) -> f64 {
    if r < 1e-15 {
        0.0
    } else {
        -1.0 / (4.0 * std::f64::consts::PI * r)
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
        assert!((g - (-1.0 / (4.0 * std::f64::consts::PI))).abs() < 1e-10);
    }

    #[test]
    fn test_helmholtz_gfn() {
        let g = green_fn_helmholtz(1.0, 1.0);
        assert!(g.norm() > 0.0);
    }
}
