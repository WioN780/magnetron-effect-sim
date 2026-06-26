use magnetron_core::config::MagnetronConfig;
use magnetron_core::diagnostics::run_diagnostics_sweep;

#[test]
fn test_sweep_resolution_convergence() {
    let mut config = MagnetronConfig::default();
    config.steps_per_gyroperiod = 32;

    let u_a_vals = vec![40.0, 50.0, 60.0, 70.0, 80.0];
    let n_particles = 25000;
    let max_steps = 1000;
    
    println!("=== Sweep Resolution Convergence (N=25000) ===");
    for &points in &[10, 20, 40, 80] {
        let result = run_diagnostics_sweep(&config, &u_a_vals, points, n_particles, max_steps);
        println!(
            "Points: {:2}, Recovered e/m: {:.4e} C/kg, v0: {:.4e} m/s (alpha: {:.4e}, beta: {:.4e})",
            points, result.recovered_e_m, result.recovered_v0, result.alpha, result.beta
        );
    }
}
