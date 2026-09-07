//! Rigid body dynamics: Euler equations, moment of inertia.

/// Compute the inertia tensor for a set of point masses.
///
/// I_{ij} = Σ_m m_a (|r_a|² δ_{ij} - r_{a,i} r_{a,j})
pub fn inertia_tensor(
    masses: &[f64],
    positions: &[[f64; 3]],
) -> [[f64; 3]; 3] {
    let mut i = [[0.0; 3]; 3];

    for (m, r) in masses.iter().zip(positions.iter()) {
        let r_sq = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
        for j in 0..3 {
            for k in 0..3 {
                if j == k {
                    i[j][k] += m * (r_sq - r[j] * r[k]);
                } else {
                    i[j][k] -= m * r[j] * r[k];
                }
            }
        }
    }

    i
}

/// Angular momentum: L = I ω.
pub fn angular_momentum(
    inertia: &[[f64; 3]; 3],
    omega: &[f64],
) -> Vec<f64> {
    (0..3).map(|i| (0..3).map(|j| inertia[i][j] * omega[j]).sum()).collect()
}

/// Rotational kinetic energy: T = ½ ω · L = ½ ωᵀ I ω.
pub fn rotational_kinetic_energy(
    inertia: &[[f64; 3]; 3],
    omega: &[f64],
) -> f64 {
    0.5 * (0..3).map(|i| omega[i] * (0..3).map(|j| inertia[i][j] * omega[j]).sum::<f64>()).sum::<f64>()
}

/// Euler's equations for a rigid body (torque-free):
/// I₁ ω̇₁ = (I₂ - I₃) ω₂ ω₃
/// I₂ ω̇₂ = (I₃ - I₁) ω₃ ω₁
/// I₃ ω̇₃ = (I₁ - I₂) ω₁ ω₂
pub fn euler_equations(
    omega: &[f64],
    moments: &[f64],  // Principal moments of inertia [I₁, I₂, I₃]
    torque: &[f64],   // External torque [N₁, N₂, N₃]
) -> Vec<f64> {
    let mut omega_dot = vec![0.0; 3];
    omega_dot[0] = (torque[0] + (moments[1] - moments[2]) * omega[1] * omega[2]) / moments[0];
    omega_dot[1] = (torque[1] + (moments[2] - moments[0]) * omega[2] * omega[0]) / moments[1];
    omega_dot[2] = (torque[2] + (moments[0] - moments[1]) * omega[0] * omega[1]) / moments[2];
    omega_dot
}

/// Torque-free Euler equations (special case).
pub fn torque_free(
    omega: &[f64],
    moments: &[f64],
) -> Vec<f64> {
    euler_equations(omega, moments, &[0.0, 0.0, 0.0])
}

/// Precession frequency for a symmetric top (I₁ = I₂).
pub fn precession_frequency(
    omega3: f64,
    torque_magnitude: f64,
    moments: &[f64],
) -> f64 {
    let i1 = moments[0];
    let i3 = moments[2];
    torque_magnitude / (i3 * omega3 * (i3 - i1).abs())
}

/// Euler angles to rotation matrix.
pub fn euler_to_rotation(
    phi: f64,    // Precession
    theta: f64,  // Nutation
    psi: f64,    // Intrinsic rotation
) -> [[f64; 3]; 3] {
    let (sp, cp) = phi.sin_cos();
    let (st, ct) = theta.sin_cos();
    let (ss, cs) = psi.sin_cos();

    [
        [cp * cs - sp * ct * ss, -cp * ss - sp * ct * cs, sp * st],
        [sp * cs + cp * ct * ss, -sp * ss + cp * ct * cs, -cp * st],
        [st * ss, st * cs, ct],
    ]
}

/// Quaternion multiplication (for efficient rotations).
pub fn quaternion_multiply(q1: &[f64; 4], q2: &[f64; 4]) -> [f64; 4] {
    [
        q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3],
        q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2],
        q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1],
        q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0],
    ]
}

/// Quaternion to rotation matrix.
pub fn quaternion_to_rotation(q: &[f64; 4]) -> [[f64; 3]; 3] {
    let (w, x, y, z) = (q[0], q[1], q[2], q[3]);
    let (w2, x2, y2, z2) = (w * w, x * x, y * y, z * z);
    let (wx, wy, wz) = (w * x, w * y, w * z);
    let (xy, xz, yz) = (x * y, x * z, y * z);
    let norm = w2 + x2 + y2 + z2;
    let s = 2.0 / norm;

    [
        [1.0 - s * (y2 + z2), s * (xy - wz), s * (xz + wy)],
        [s * (xy + wz), 1.0 - s * (x2 + z2), s * (yz - wx)],
        [s * (xz - wy), s * (yz + wx), 1.0 - s * (x2 + y2)],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_inertia_point_mass() {
        let masses = [1.0];
        let positions = [[1.0, 0.0, 0.0]];
        let i = inertia_tensor(&masses, &positions);
        // For a point mass at (1,0,0): I = diag(0, 1, 1)
        assert!((i[0][0]).abs() < 1e-10);
        assert!((i[1][1] - 1.0).abs() < 1e-10);
        assert!((i[2][2] - 1.0).abs() < 1e-10);
    }
}
