use clap::{Parser, Subcommand};
use magnetron_core::config::MagnetronConfig;
use magnetron_core::diagnostics::{run_diagnostics_sweep, run_coaxial_simulation_trajectories};
use std::fs;
use std::path::Path;
use serde::Serialize;

#[derive(Parser)]
#[command(name = "reference-cli")]
#[command(author, version, about = "Magnetron Reference CLI", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// Run parameter sweep
    Sweep {
        /// Mode of sweep: idealized or scl
        #[arg(long, default_value = "idealized")]
        mode: String,

        /// Output JSON file path
        #[arg(short, long, default_value = "data/runs/idealized_m4.json")]
        output: String,

        /// Number of particles
        #[arg(short, long, default_value_t = 25000)]
        particles: usize,

        /// Steps per gyroperiod
        #[arg(short, long, default_value_t = 32)]
        steps: u32,

        /// Number of sweep points for solenoid current Ic
        #[arg(short, long, default_value_t = 40)]
        points: usize,
    },
    /// Run convergence analysis (stub)
    Converge,
    /// Export replay data (stub)
    #[command(name = "replay-export")]
    ReplayExport,
}

#[derive(Serialize)]
struct TrajectoryRun {
    i_c: f64,
    trajectories: Vec<Vec<[f64; 3]>>,
}

#[derive(Serialize)]
struct SampleTrajectories {
    u_a: f64,
    passing: TrajectoryRun,
    cutoff: TrajectoryRun,
}

#[derive(Serialize)]
struct OutputDataset {
    recovered_e_m: f64,
    recovered_v0: f64,
    alpha: f64,
    beta: f64,
    u_a_sweeps: serde_json::Value,
    sample_trajectories: SampleTrajectories,
}

fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Some(Commands::Sweep { mode, output, particles, steps, points }) => {
            println!("Starting parameter sweep in mode: {}", mode);
            println!("Parameters: particles={}, steps={}, points={}", particles, steps, points);

            let mut config = MagnetronConfig::default();
            config.steps_per_gyroperiod = *steps;

            let u_a_vals = vec![40.0, 50.0, 60.0, 70.0, 80.0];
            let max_steps = 1000;

            // 1. Run the diagnostics sweep to recover parameters and get I_a(I_c) curves
            let diag_result = run_diagnostics_sweep(&config, &u_a_vals, *points, *particles, max_steps);

            // 2. Generate sample trajectories at Ua = 40 V
            // We'll capture a passing run (Ic = 0.0 A) and a cutoff run (Ic = 0.8 A)
            println!("Generating sample trajectories at Ua = 40 V...");
            
            let mut passing_config = config.clone();
            passing_config.anode_voltage = 40.0;
            passing_config.solenoid_current = 0.0;
            let passing_trajectories = run_coaxial_simulation_trajectories(&passing_config, *particles, max_steps, 5);

            let mut cutoff_config = config.clone();
            cutoff_config.anode_voltage = 40.0;
            cutoff_config.solenoid_current = 0.8;
            let cutoff_trajectories = run_coaxial_simulation_trajectories(&cutoff_config, *particles, max_steps, 5);

            let sample_trajectories = SampleTrajectories {
                u_a: 40.0,
                passing: TrajectoryRun {
                    i_c: 0.0,
                    trajectories: passing_trajectories,
                },
                cutoff: TrajectoryRun {
                    i_c: 0.8,
                    trajectories: cutoff_trajectories,
                },
            };

            // 3. Assemble the output dataset
            let u_a_sweeps_json = serde_json::to_value(&diag_result.sweeps).unwrap();
            let dataset = OutputDataset {
                recovered_e_m: diag_result.recovered_e_m,
                recovered_v0: diag_result.recovered_v0,
                alpha: diag_result.alpha,
                beta: diag_result.beta,
                u_a_sweeps: u_a_sweeps_json,
                sample_trajectories,
            };

            // 4. Write to JSON
            let out_path = Path::new(output);
            if let Some(parent) = out_path.parent() {
                fs::create_dir_all(parent).unwrap();
            }

            let file = fs::File::create(out_path).unwrap();
            serde_json::to_writer_pretty(file, &dataset).unwrap();

            println!("Sweep completed! Output written to {}", output);
            println!("Recovered e/m: {:.4e} C/kg", diag_result.recovered_e_m);
            println!("Recovered v0: {:.4e} m/s", diag_result.recovered_v0);
        }
        Some(Commands::Converge) => {
            println!("Converge subcommand (stub)");
        }
        Some(Commands::ReplayExport) => {
            println!("Replay-export subcommand (stub)");
        }
        None => {
            println!("Use --help to see available commands.");
        }
    }
}
