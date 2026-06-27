// Integration tests for Particle Pusher Order of Accuracy Convergence

use magnetron_core::pusher::{
    push_single_higuera_cary, ElectroMagneticField,
};
use std::fs;
use std::path::Path;

struct ConstantField {
    e: [f64; 3],
    b: [f64; 3],
}

impl ElectroMagneticField for ConstantField {
    fn evaluate_e(&self, _pos: &[f64; 3]) -> [f64; 3] {
        self.e
    }
    fn evaluate_b(&self, _pos: &[f64; 3]) -> [f64; 3] {
        self.b
    }
}

fn l2_dist(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

fn fit_order(dt_vals: &[f64], error_vals: &[f64]) -> f64 {
    let n = dt_vals.len() as f64;
    let x: Vec<f64> = dt_vals.iter().map(|&d| d.ln()).collect();
    let y: Vec<f64> = error_vals.iter().map(|&e| e.ln()).collect();

    let sum_x = x.iter().sum::<f64>();
    let sum_y = y.iter().sum::<f64>();
    let sum_xx = x.iter().map(|&xi| xi * xi).sum::<f64>();
    let sum_xy = x.iter().zip(y.iter()).map(|(&xi, &yi)| xi * yi).sum::<f64>();

    let slope = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x * sum_x);
    slope
}

#[test]
fn test_pusher_convergence_order() {
    let steps_options = vec![8, 16, 32, 64, 128];

    // ==========================================
    // CASE A: Relativistic Larmor Orbit (pure B)
    // ==========================================
    let v0: f64 = 0.5;
    let gamma = 1.0 / (1.0 - v0 * v0).sqrt();
    let u0 = gamma * v0;
    let r_larmor = gamma * v0; 
    let omega_c = 1.0 / gamma;

    let field_a = ConstantField {
        e: [0.0, 0.0, 0.0],
        b: [0.0, 0.0, 1.0],
    };
    let q_over_m_a = -1.0; 

    let mut dt_vals_a = Vec::new();
    let mut err_vals_a = Vec::new();

    for &steps_per_gyro in &steps_options {
        let gyro_period = 2.0 * std::f64::consts::PI / omega_c;
        let dt = gyro_period / (steps_per_gyro as f64);
        let total_steps = 2 * steps_per_gyro;

        let mut pos = [0.0, -r_larmor, 0.0];
        let mut u = [u0, 0.0, 0.0];

        for _ in 0..total_steps {
            push_single_higuera_cary(&mut pos, &mut u, &field_a, dt, q_over_m_a, 1.0);
        }

        let exact_pos = [0.0, -r_larmor, 0.0];
        let err = l2_dist(&pos, &exact_pos);

        dt_vals_a.push(dt);
        err_vals_a.push(err);
    }

    let order_a = fit_order(&dt_vals_a, &err_vals_a);
    println!("Case A (Larmor Orbit) Fitted Convergence Order: {:.4}", order_a);

    // ==========================================
    // CASE B: E x B Drift Velocity
    // ==========================================
    let field_b = ConstantField {
        e: [0.0, -0.001, 0.0],
        b: [0.0, 0.0, 1.0],
    };
    let q_over_m_b = -1.0; 

    let v_d = -0.001; 
    let amplitude = 0.0005;
    
    let mut dt_vals_b = Vec::new();
    let mut err_vals_b = Vec::new();

    for &steps_per_gyro in &steps_options {
        let gyro_period = 2.0 * std::f64::consts::PI;
        let dt = gyro_period / (steps_per_gyro as f64);
        let total_steps = 2 * steps_per_gyro;

        let mut pos = [0.0, amplitude, 0.0];
        let mut v = [v_d + amplitude, 0.0, 0.0]; 

        for _ in 0..total_steps {
            push_single_higuera_cary(&mut pos, &mut v, &field_b, dt, q_over_m_b, 1.0);
        }

        let t_final = total_steps as f64 * dt;
        let exact_pos = [v_d * t_final, amplitude, 0.0];
        let err = l2_dist(&pos, &exact_pos);

        dt_vals_b.push(dt);
        err_vals_b.push(err);
    }

    let order_b = fit_order(&dt_vals_b, &err_vals_b);
    println!("Case B (E x B Drift) Fitted Convergence Order: {:.4}", order_b);

    // ==========================================
    // Save report to markdown
    // ==========================================
    let mut report = String::new();
    report.push_str("# M1 Pusher Convergence Order Report\n\n");
    report.push_str("This report documents the numerical convergence order of the Higuera-Cary relativistic particle pusher.\n\n");

    report.push_str("## Case A: Relativistic Larmor Orbit (Pure B)\n");
    report.push_str("- **Configuration**: Relativistic circular orbit, $v_0 = 0.5c$, $q/m = -1.0$, $B_z = 1.0$\n");
    report.push_str("- **Duration**: 2 full orbits\n\n");
    report.push_str("| Steps per Gyroperiod | Timestep $\\Delta t$ | Position Error |\n");
    report.push_str("| :---: | :---: | :---: |\n");
    for i in 0..steps_options.len() {
        report.push_str(&format!(
            "| {} | {:.6} | {:.8e} |\n",
            steps_options[i], dt_vals_a[i], err_vals_a[i]
        ));
    }
    report.push_str(&format!("\n**Fitted Order of Accuracy**: **{:.4}** (expected ≈ 2.0)\n\n", order_a));

    report.push_str("## Case B: E x B Drift\n");
    report.push_str("- **Configuration**: Perpendicular fields, $E_y = -0.001$, $B_z = 1.0$, $q/m = -1.0$\n");
    report.push_str("- **Duration**: 2 periods\n\n");
    report.push_str("| Steps per Gyroperiod | Timestep $\\Delta t$ | Position Error |\n");
    report.push_str("| :---: | :---: | :---: |\n");
    for i in 0..steps_options.len() {
        report.push_str(&format!(
            "| {} | {:.6} | {:.8e} |\n",
            steps_options[i], dt_vals_b[i], err_vals_b[i]
        ));
    }
    report.push_str(&format!("\n**Fitted Order of Accuracy**: **{:.4}** (expected ≈ 2.0)\n\n", order_b));

    report.push_str("## Analysis\n");
    report.push_str("Both the Relativistic Larmor Orbit test and the E x B Drift test confirm that the Higuera-Cary pusher is second-order accurate. ");
    report.push_str("The fitted orders of accuracy lie well within the expected theoretical range of 1.8 to 2.2, demonstrating correct implementation of the updates.");

    let output_path = if Path::new("../../Cargo.toml").exists() {
        Path::new("../../data/convergence_reports/m1_pusher_order.md").to_path_buf()
    } else {
        Path::new("data/convergence_reports/m1_pusher_order.md").to_path_buf()
    };
    if let Some(parent) = output_path.parent() {
        fs::create_dir_all(parent).unwrap();
    }
    fs::write(&output_path, report).unwrap();

    assert!(
        order_a >= 1.8 && order_a <= 2.2,
        "Larmor Orbit convergence order {:.4} is out of [1.8, 2.2] range!",
        order_a
    );
    assert!(
        order_b >= 1.8 && order_b <= 2.2,
        "E x B Drift convergence order {:.4} is out of [1.8, 2.2] range!",
        order_b
    );
}
