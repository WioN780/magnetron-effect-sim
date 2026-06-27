// (PhaseSpace is defined in particles.rs)

/// Evaluates electric and magnetic fields at normalized spatial coordinates.
pub trait ElectroMagneticField {
    /// Evaluates normalized electric field [Ex, Ey, Ez] at normalized position [x, y, z]
    fn evaluate_e(&self, pos: &[f64; 3]) -> [f64; 3];
    /// Evaluates normalized magnetic field [Bx, By, Bz] at normalized position [x, y, z]
    fn evaluate_b(&self, pos: &[f64; 3]) -> [f64; 3];
}


// Vector algebra helper utilities for [f64; 3]
#[inline]
pub fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
pub fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
pub fn add(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
pub fn scale(a: &[f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
pub fn norm_sq(a: &[f64; 3]) -> f64 {
    dot(a, a)
}

/// Computes the time step dt from steps per gyroperiod and cyclotron frequency.
/// dt = 2 * pi * m / (steps_per_gyroperiod * |q| * B)
pub fn compute_dt(steps_per_gyroperiod: u32, b_magnitude: f64, q: f64, m: f64) -> f64 {
    let steps = steps_per_gyroperiod as f64;
    let omega_c = (q.abs() * b_magnitude) / m;
    if omega_c.abs() < 1e-12 {
        // Fallback if field is zero
        2.0 * std::f64::consts::PI / steps
    } else {
        2.0 * std::f64::consts::PI / (steps * omega_c)
    }
}

/// Relativistic Higuera-Cary Boris pusher for a single particle.
///
/// Updates position [x, y, z] and relativistic momentum u = gamma * v / c.
/// charge_m_ratio is q/m in normalized units.
pub fn push_single_higuera_cary(
    pos: &mut [f64; 3],
    u: &mut [f64; 3],
    field: &dyn ElectroMagneticField,
    dt: f64,
    charge_m_ratio: f64,
    pos_scale: f64,
) {
    let e = field.evaluate_e(pos);
    let b = field.evaluate_b(pos);

    // 1. Half-step electric acceleration: u_minus = u_i + epsilon
    let epsilon = scale(&e, 0.5 * charge_m_ratio * dt);
    let u_minus = add(u, &epsilon);

    // 2. Compute relativistic factor gamma_minus
    let u_minus_sq = norm_sq(&u_minus);
    let gamma_minus = (1.0 + u_minus_sq).sqrt();

    // 3. Compute rotation vector beta
    let beta = scale(&b, 0.5 * charge_m_ratio * dt);
    let beta_sq = norm_sq(&beta);
    let beta_dot_u_minus = dot(&beta, &u_minus);

    // 4. Compute gamma_new using the Higuera-Cary biquadratic formula
    let diff = gamma_minus * gamma_minus - beta_sq;
    let radical = (diff * diff + 4.0 * (beta_sq + beta_dot_u_minus * beta_dot_u_minus)).sqrt();
    let gamma_new_sq = 0.5 * (diff + radical);
    let gamma_new = gamma_new_sq.sqrt();

    // 5. Perform magnetic rotation
    let t_vec = scale(&beta, 1.0 / gamma_new);
    let t_sq = norm_sq(&t_vec);
    let s_vec = scale(&t_vec, 2.0 / (1.0 + t_sq));

    // u_plus = u_minus + (u_minus + u_minus x t) x s
    let u_minus_cross_t = cross(&u_minus, &t_vec);
    let u_prime = add(&u_minus, &u_minus_cross_t);
    let u_prime_cross_s = cross(&u_prime, &s_vec);
    let u_plus = add(&u_minus, &u_prime_cross_s);

    // 6. Final half-step electric acceleration
    let u_final = add(&u_plus, &epsilon);

    // 7. Update position using final velocity v_final = u_final / gamma_final
    let u_final_sq = norm_sq(&u_final);
    let gamma_final = (1.0 + u_final_sq).sqrt();
    let v_final = scale(&u_final, 1.0 / gamma_final);

    *u = u_final;
    *pos = add(pos, &scale(&v_final, pos_scale * dt));
}

/// Non-relativistic Boris pusher for a single particle (for reference comparison).
///
/// Updates position [x, y, z] and velocity v (which acts as momentum here, scaled by c).
/// charge_m_ratio is q/m in normalized units.
pub fn push_single_non_relativistic_boris(
    pos: &mut [f64; 3],
    v: &mut [f64; 3],
    field: &dyn ElectroMagneticField,
    dt: f64,
    charge_m_ratio: f64,
    pos_scale: f64,
) {
    let e = field.evaluate_e(pos);
    let b = field.evaluate_b(pos);

    // 1. Half-step electric acceleration: v_minus = v_i + epsilon
    let epsilon = scale(&e, 0.5 * charge_m_ratio * dt);
    let v_minus = add(v, &epsilon);

    // 2. Compute rotation vector t
    let t_vec = scale(&b, 0.5 * charge_m_ratio * dt);
    let t_sq = norm_sq(&t_vec);
    let s_vec = scale(&t_vec, 2.0 / (1.0 + t_sq));

    // 3. Perform magnetic rotation: v_plus = v_minus + (v_minus + v_minus x t) x s
    let v_minus_cross_t = cross(&v_minus, &t_vec);
    let v_prime = add(&v_minus, &v_minus_cross_t);
    let v_prime_cross_s = cross(&v_prime, &s_vec);
    let v_plus = add(&v_minus, &v_prime_cross_s);

    // 4. Final half-step electric acceleration
    let v_final = add(&v_plus, &epsilon);

    *v = v_final;
    *pos = add(pos, &scale(&v_final, pos_scale * dt));
}
// (Batch pushers are defined in particles.rs)


#[cfg(test)]
mod tests {
    use super::*;

    struct SimpleField {
        e_val: [f64; 3],
        b_val: [f64; 3],
    }

    impl ElectroMagneticField for SimpleField {
        fn evaluate_e(&self, _pos: &[f64; 3]) -> [f64; 3] { self.e_val }
        fn evaluate_b(&self, _pos: &[f64; 3]) -> [f64; 3] { self.b_val }
    }

    #[test]
    fn test_non_relativistic_limit() {
        let field = SimpleField {
            e_val: [1e-8, -2e-8, 3e-8],
            b_val: [1e-6, 2e-6, -3e-6],
        };

        // Initialize particle with very small velocity (v << c), so gamma ≈ 1
        let mut pos_hc = [0.0, 0.0, 0.0];
        let mut u_hc = [1e-6, -2e-6, 1.5e-6]; // relativistic momentum (u = gamma * v ≈ v)

        let mut pos_nr = [0.0, 0.0, 0.0];
        let mut v_nr = [1e-6, -2e-6, 1.5e-6]; // non-relativistic velocity

        let dt = 0.05;
        let q_over_m = -1.0;

        // Push both 100 steps
        for _ in 0..100 {
            push_single_higuera_cary(&mut pos_hc, &mut u_hc, &field, dt, q_over_m, 1.0);
            push_single_non_relativistic_boris(&mut pos_nr, &mut v_nr, &field, dt, q_over_m, 1.0);
        }

        // Verify that Higuera-Cary reduces to non-relativistic Boris within high precision
        for i in 0..3 {
            assert!((pos_hc[i] - pos_nr[i]).abs() < 1e-12, "pos[{}] differs: HC={}, NR={}", i, pos_hc[i], pos_nr[i]);
            assert!((u_hc[i] - v_nr[i]).abs() < 1e-12, "vel[{}] differs: HC={}, NR={}", i, u_hc[i], v_nr[i]);
        }
    }
}
