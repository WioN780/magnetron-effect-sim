// Magnetron Core Physics Engine Library

pub mod constants;
pub mod config;
pub mod pusher;
pub mod particles;
pub mod diagnostics;
pub mod poisson;

pub use config::{MagnetronConfig, Normalization};

pub fn get_physics_info() -> &'static str {
    "Magnetron Core Physics Engine"
}

/// A pure Rust simulation runner that encapsulates the state of the magnetron simulation.
/// This struct has no JS/WASM bindings and produces pure physics outputs.
pub struct Simulation {
    pub config: MagnetronConfig,
    pub norm: Normalization,
    pub field: particles::CoaxialDiodeField,
    pub phase_space: particles::PhaseSpace,
    pub active: Vec<bool>,
    pub hit_anode: Vec<bool>,
    pub step_count: usize,
    pub dt: f64,
    pub pos_scale: f64,
    pub charge_m_ratio: f64,
    pub r_c_norm: f64,
}

impl Simulation {
    /// Creates a new simulation with a given configuration and number of particles.
    pub fn new(config: MagnetronConfig, num_particles: usize) -> Self {
        let norm = config.normalization();
        let ln_ratio = (config.anode_radius / config.cathode_radius).ln();
        
        let e_coeff = -config.anode_voltage / (norm.e_0 * norm.l_0 * ln_ratio);
        let b_normalized_z = 1.0; 
        
        let field = particles::CoaxialDiodeField {
            e_coeff,
            b_z: b_normalized_z,
        };
        
        let mut phase_space = particles::PhaseSpace::new(num_particles);
        
        // Seed with a specific value for stability in reproduction
        diagnostics::initialize_cathode_particles_seeded(&mut phase_space, &config, &norm, 123456);
        
        let active = vec![true; num_particles];
        let hit_anode = vec![false; num_particles];
        
        let r_c_norm = norm.normalize_length(config.cathode_radius);
        let charge_m_ratio = -1.0;
        let dt = 2.0 * std::f64::consts::PI / (config.steps_per_gyroperiod as f64);
        let pos_scale = norm.v_0 * norm.t_0 / norm.l_0;
        
        Self {
            config,
            norm,
            field,
            phase_space,
            active,
            hit_anode,
            step_count: 0,
            dt,
            pos_scale,
            charge_m_ratio,
            r_c_norm,
        }
    }

    /// Advances the simulation by one step. Returns true if there are still active particles.
    pub fn step(&mut self) -> bool {
        particles::push_batch_higuera_cary(
            &mut self.phase_space,
            &self.active,
            &self.field,
            self.dt,
            self.charge_m_ratio,
            self.pos_scale,
        );
        
        self.step_count += 1;
        let num_particles = self.phase_space.num_particles();
        let pos_view = self.phase_space.positions.view();
        let mut active_count = 0;
        
        for i in 0..num_particles {
            if self.active[i] {
                let x = pos_view[[0, i]];
                let y = pos_view[[1, i]];
                let r = (x*x + y*y).sqrt();
                
                if r >= 1.0 {
                    self.hit_anode[i] = true;
                    self.active[i] = false;
                } else if r <= self.r_c_norm && self.step_count > 2 {
                    self.active[i] = false;
                } else {
                    active_count += 1;
                }
            }
        }
        
        active_count > 0
    }

    /// Returns a flat vector of physical coordinates (x, y, z) for all particles in meters.
    pub fn get_positions(&self) -> Vec<f64> {
        let num_particles = self.phase_space.num_particles();
        let mut result = Vec::with_capacity(num_particles * 3);
        let pos_view = self.phase_space.positions.view();
        let l_0 = self.norm.l_0;
        
        for i in 0..num_particles {
            result.push(pos_view[[0, i]] * l_0);
            result.push(pos_view[[1, i]] * l_0);
            result.push(pos_view[[2, i]] * l_0);
        }
        
        result
    }

    /// Returns a flat vector of physical velocities (vx, vy, vz) for all particles in m/s.
    pub fn get_velocities(&self) -> Vec<f64> {
        let num_particles = self.phase_space.num_particles();
        let mut result = Vec::with_capacity(num_particles * 3);
        let mom_view = self.phase_space.momenta.view();
        let c = constants::C;
        
        for i in 0..num_particles {
            let ux = mom_view[[0, i]];
            let uy = mom_view[[1, i]];
            let uz = mom_view[[2, i]];
            let u_sq = ux * ux + uy * uy + uz * uz;
            let gamma = (1.0 + u_sq).sqrt();
            result.push((ux / gamma) * c);
            result.push((uy / gamma) * c);
            result.push((uz / gamma) * c);
        }
        
        result
    }

    /// Returns the active state for all particles.
    pub fn get_active_states(&self) -> Vec<bool> {
        self.active.clone()
    }

    /// Computes the current collected at the anode.
    pub fn get_anode_current(&self) -> f64 {
        let num_particles = self.phase_space.num_particles();
        let hit_count = self.hit_anode.iter().filter(|&&h| h).count();
        let fraction = (hit_count as f64) / (num_particles as f64);
        fraction * self.config.max_anode_current
    }
}
