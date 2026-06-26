use magnetron_core::config::MagnetronConfig;
use magnetron_core::diagnostics::run_diagnostics_sweep;

#[test]
fn test_sweep_step_convergence() {
    let mut config = MagnetronConfig::default();
    config.steps_per_gyroperiod = 32;

    let u_a_vals = vec![40.0, 50.0, 60.0, 70.0, 80.0];
    let n_particles = 5000;
    let max_steps = 1000;
    
    // Print the details of the sweep for Ua = 40.0V to inspect the curve
    let result = run_diagnostics_sweep(&config, &u_a_vals, 20, n_particles, max_steps);
    
    println!("=== Ua = 40.0V Sweep Curve ===");
    for sweep in &result.sweeps {
        if (sweep.u_a - 40.0).abs() < 1e-6 {
            for i in 0..sweep.i_c_values.len() {
                println!("Ic = {:.6} A, Ia = {:.6} A (fraction = {:.4})", 
                    sweep.i_c_values[i], sweep.i_a_values[i], sweep.i_a_values[i] / config.max_anode_current
                );
            }
            println!("Detected inflection Ic = {:.6} A, Bk = {:.6} T", sweep.inflection_i_c, sweep.inflection_b_k);
        }
    }

    println!("=== Diagnostics Results ===");
    println!("Recovered e/m: {:.6e} C/kg", result.recovered_e_m);
    println!("Recovered v0: {:.6e} m/s", result.recovered_v0);
    println!("alpha: {:.6e}, beta: {:.6e}", result.alpha, result.beta);
}
