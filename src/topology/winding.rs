//! Winding numbers and topological charges.

use num_complex::Complex64;

/// Winding number of a complex curve around the origin.
/// n = (1/2πi) ∮ dz/z
///
/// For a discrete set of points on the curve.
pub fn winding_number(
    curve: &[(f64, f64)],  // Points (Re, Im) on the closed curve
) -> f64 {
    let n = curve.len();
    let mut total_angle = 0.0;

    for i in 0..n {
        let j = (i + 1) % n;
        let (x1, y1) = curve[i];
        let (x2, y2) = curve[j];

        let angle1 = y1.atan2(x1);
        let angle2 = y2.atan2(x2);
        let mut d_angle = angle2 - angle1;

        // Normalize to [-π, π]
        while d_angle > std::f64::consts::PI {
            d_angle -= 2.0 * std::f64::consts::PI;
        }
        while d_angle < -std::f64::consts::PI {
            d_angle += 2.0 * std::f64::consts::PI;
        }

        total_angle += d_angle;
    }

    (total_angle / (2.0 * std::f64::consts::PI)).round()
}

/// Discrete winding number from complex overlap phases.
/// n = (1/2π) Σ arg(⟨z_i|z_{i+1}⟩)
pub fn winding_number_discrete(
    overlaps: &[Complex64],  // Overlaps between adjacent points
) -> f64 {
    let total_phase: f64 = overlaps.iter().map(|z| z.arg()).sum();
    (total_phase / (2.0 * std::f64::consts::PI)).round()
}

/// Winding number of a real function f(x) = (f₁(x), f₂(x)) around the origin.
pub fn winding_number_real(
    f1: &[f64],
    f2: &[f64],
) -> f64 {
    let curve: Vec<(f64, f64)> = f1.iter().zip(f2.iter()).map(|(&a, &b)| (a, b)).collect();
    winding_number(&curve)
}

/// Topological charge (skyrmion number) for a 2D field:
/// Q = (1/4π) ∫ n̂ · (∂_x n̂ × ∂_y n̂) d²x
/// where n̂ is a unit vector field.
///
/// `n_field` is indexed `[ix][iy]` and each entry is the 3-component unit
/// vector n̂ at that grid point, so the grid may be any nx × ny. (The previous
/// signature was `&[[[f64; 3]; 3]]`, which pinned the second axis to exactly
/// three columns and made the function uncallable for a real grid.)
pub fn skyrmion_number(
    n_field: &[Vec<[f64; 3]>],  // n̂(x, y) unit vectors on an nx × ny grid
    dx: f64,
    dy: f64,
) -> f64 {
    let nx = n_field.len();
    let ny = n_field[0].len();
    let mut charge = 0.0;

    for i in 1..nx - 1 {
        for j in 1..ny - 1 {
            // ∂_x n̂
            let dnx = [
                (n_field[i + 1][j][0] - n_field[i - 1][j][0]) / (2.0 * dx),
                (n_field[i + 1][j][1] - n_field[i - 1][j][1]) / (2.0 * dx),
                (n_field[i + 1][j][2] - n_field[i - 1][j][2]) / (2.0 * dx),
            ];
            // ∂_y n̂
            let dny = [
                (n_field[i][j + 1][0] - n_field[i][j - 1][0]) / (2.0 * dy),
                (n_field[i][j + 1][1] - n_field[i][j - 1][1]) / (2.0 * dy),
                (n_field[i][j + 1][2] - n_field[i][j - 1][2]) / (2.0 * dy),
            ];
            // Cross product
            let cross = [
                dnx[1] * dny[2] - dnx[2] * dny[1],
                dnx[2] * dny[0] - dnx[0] * dny[2],
                dnx[0] * dny[1] - dnx[1] * dny[0],
            ];
            // Dot with n̂
            let dot = n_field[i][j][0] * cross[0]
                + n_field[i][j][1] * cross[1]
                + n_field[i][j][2] * cross[2];
            charge += dot;
        }
    }

    charge / (4.0 * std::f64::consts::PI) * dx * dy
}

/// Instanton number in lattice gauge theory: n = (1/32π²) Σ_P ε_{μνρσ} tr(F_{μν} F_{ρσ}).
pub fn instanton_number_lattice(
    wilson_loops: &Vec<Vec<Vec<Vec<Vec<Vec<Complex64>>>>>>,
    n_lattice: usize,
) -> f64 {
    let mut charge = 0.0;

    for i in 0..n_lattice {
        for j in 0..n_lattice {
            for k in 0..n_lattice {
                for l in 0..n_lattice {
                    // Q = Σ_P Im tr(U_P)
                    for mu in 0..4 {
                        for nu in 0..4 {
                            if mu == nu { continue; }
                            charge += wilson_loops[i][j][k][l][mu][nu].arg();
                        }
                    }
                }
            }
        }
    }

    charge / (32.0 * std::f64::consts::PI * std::f64::consts::PI)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_winding_unit_circle() {
        // Unit circle: winding number should be 1
        let curve: Vec<(f64, f64)> = (0..360)
            .map(|i| {
                let angle = i as f64 * std::f64::consts::PI / 180.0;
                (angle.cos(), angle.sin())
            })
            .collect();
        let n = winding_number(&curve);
        assert!((n - 1.0).abs() < 0.1, "Expected 1, got {}", n);
    }

    #[test]
    fn test_winding_origin() {
        // Curve passing through origin
        let curve = vec![(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)];
        let n = winding_number(&curve);
        assert!(n.abs() >= 0.9);
    }
}
