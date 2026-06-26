use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "reference-cli")]
#[command(author, version, about = "Magnetron Reference CLI", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// Run parameter sweep (stub)
    Sweep,
    /// Run convergence analysis (stub)
    Converge,
    /// Export replay data (stub)
    #[command(name = "replay-export")]
    ReplayExport,
}

fn main() {
    let cli = Cli::parse();
    
    // Call a dummy function in core to prove it depends on core
    println!("Core engine info: {}", magnetron_core::get_physics_info());

    match &cli.command {
        Some(Commands::Sweep) => {
            println!("Sweep subcommand (stub)");
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
