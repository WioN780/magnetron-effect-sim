// Diagnostics and automated parameter recovery pipeline

use crate::config::MagnetronConfig;
use crate::constants::{E, M_E, MU_0};
use crate::particles::run_coaxial_simulation;

/// Computes the central homogeneous magnetic field of the solenoid (T)
/// for a given solenoid current Ic, solenoid turns, length, and diameter.
pub fn compute_b_field(i_c: f64, turn_count: f64, length: f64, diameter: f64) -> f64 {
    (MU_0 * turn_count * i_c) / (length * length + diameter * diameter).sqrt()
}

/// Computes the solenoid current (A) for a given magnetic field B (T)
pub fn compute_i_c(b: f64, turn_count: f64, length: f64, diameter: f64) -> f64 {
    let term = (length * length + diameter * diameter).sqrt();
    b * term / (MU_0 * turn_count)
}

/// Find the inflection point of the I_a(I_c) curve.
/// The inflection point is the point of steepest descent (maximum absolute negative slope).
/// We look for the zero crossing of the second derivative in the neighborhood of the steepest descent.
pub fn find_inflection_point(i_c_vals: &[f64], i_a_vals: &[f64]) -> Option<f64> {
    if i_c_vals.len() < 3 || i_a_vals.len() < 3 {
        return None;
    }
    let h = i_c_vals[1] - i_c_vals[0];
    let mut min_slope = f64::MAX;
    let mut min_idx = 1;

    for i in 1..(i_c_vals.len() - 1) {
        let slope = (i_a_vals[i + 1] - i_a_vals[i - 1]) / (2.0 * h);
        if slope < min_slope {
            min_slope = slope;
            min_idx = i;
        }
    }

    let d2 = |idx: usize| -> f64 {
        (i_a_vals[idx + 1] - 2.0 * i_a_vals[idx] + i_a_vals[idx - 1]) / (h * h)
    };

    let d2_prev = d2(min_idx - 1);
    let d2_curr = d2(min_idx);
    let d2_next = d2(min_idx + 1);

    if d2_prev * d2_curr <= 0.0 {
        let x_prev = i_c_vals[min_idx - 1];
        let x_curr = i_c_vals[min_idx];
        let t = d2_prev / (d2_prev - d2_curr);
        Some(x_prev + t * (x_curr - x_prev))
    } else if d2_curr * d2_next <= 0.0 {
        let x_curr = i_c_vals[min_idx];
        let x_next = i_c_vals[min_idx + 1];
        let t = d2_curr / (d2_curr - d2_next);
        Some(x_curr + t * (x_next - x_curr))
    } else {
        Some(i_c_vals[min_idx])
    }
}

/// Perform least-squares linear regression on y = alpha * x + beta
/// Returns (alpha, beta)
pub fn linear_regression(x: &[f64], y: &[f64]) -> (f64, f64) {
    let n = x.len() as f64;
    let sum_x = x.iter().sum::<f64>();
    let sum_y = y.iter().sum::<f64>();
    let sum_xx = x.iter().map(|&xi| xi * xi).sum::<f64>();
    let sum_xy = x.iter().zip(y.iter()).map(|(&xi, &yi)| xi * yi).sum::<f64>();

    let slope = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x * sum_x);
    let intercept = (sum_y - slope * sum_x) / n;
    (slope, intercept)
}

/// Represents the result of a sweep for a single U_a voltage value.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct UaSweep {
    pub u_a: f64,
    pub i_c_values: Vec<f64>,
    pub i_a_values: Vec<f64>,
    pub inflection_i_c: f64,
    pub inflection_b_k: f64,
}

/// Represents the recovered parameter results from the diagnostics sweep.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct DiagnosticsResult {
    pub recovered_e_m: f64,
    pub recovered_v0: f64,
    pub alpha: f64,
    pub beta: f64,
    pub sweeps: Vec<UaSweep>,
}

/// Runs the diagnostics sweep for a given configuration, list of U_a values,
/// and number of I_c sweep points.
pub fn run_diagnostics_sweep(
    base_config: &MagnetronConfig,
    u_a_vals: &[f64],
    num_ic_points: usize,
    num_particles: usize,
    max_steps: usize,
) -> DiagnosticsResult {
    let mut sweeps = Vec::new();
    let e_over_m = E / M_E;
    let r_a = base_config.anode_radius;
    let r_c = base_config.cathode_radius;
    let geom_factor = (1.0 - (r_c / r_a).powi(2)).powi(2);

    for &u_a in u_a_vals {
        // Estimate critical B-field and current
        let b_k_est = (8.0 * u_a / (e_over_m * r_a * r_a * geom_factor)).sqrt();
        let i_c_est = compute_i_c(
            b_k_est,
            base_config.solenoid_turn_count,
            base_config.solenoid_length,
            base_config.solenoid_diameter,
        );

        // Sweep current from 0.6 * i_c_est to 1.6 * i_c_est
        let i_c_min = 0.6 * i_c_est;
        let i_c_max = 1.6 * i_c_est;

        let mut i_c_values = Vec::new();
        let mut i_a_values = Vec::new();

        for i in 0..num_ic_points {
            let i_c = i_c_min + (i_c_max - i_c_min) * (i as f64) / ((num_ic_points - 1) as f64);
            let mut config = base_config.clone();
            config.anode_voltage = u_a;
            config.solenoid_current = i_c;

            let fraction = run_coaxial_simulation(&config, num_particles, max_steps);
            let i_a = fraction * config.max_anode_current;

            i_c_values.push(i_c);
            i_a_values.push(i_a);
        }

        let inflection_i_c = find_inflection_point(&i_c_values, &i_a_values).unwrap_or(i_c_est);
        let inflection_b_k = compute_b_field(
            inflection_i_c,
            base_config.solenoid_turn_count,
            base_config.solenoid_length,
            base_config.solenoid_diameter,
        );

        sweeps.push(UaSweep {
            u_a,
            i_c_values,
            i_a_values,
            inflection_i_c,
            inflection_b_k,
        });
    }

    // Perform linear regression on B_k^2 = alpha * U_a + beta
    let x_reg: Vec<f64> = sweeps.iter().map(|s| s.u_a).collect();
    let y_reg: Vec<f64> = sweeps.iter().map(|s| s.inflection_b_k * s.inflection_b_k).collect();

    let (alpha, beta) = linear_regression(&x_reg, &y_reg);

    // Specific charge e/m = 8 / (alpha * R_a^2)
    let recovered_e_m = 8.0 / (alpha * r_a * r_a);

    // Initial velocity v0 = sqrt(8 * beta / (alpha^2 * R_a^2))
    let recovered_v0 = (8.0 * beta / (alpha * alpha * r_a * r_a)).sqrt();

    DiagnosticsResult {
        recovered_e_m,
        recovered_v0,
        alpha,
        beta,
        sweeps,
    }
}
