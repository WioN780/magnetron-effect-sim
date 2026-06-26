use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub fn hello() -> f64 {
    let _info = magnetron_core::get_physics_info();
    42.0
}

#[wasm_bindgen]
#[derive(Clone, Copy)]
pub struct WasmMagnetronConfig {
    pub anode_radius: f64,
    pub cathode_radius: f64,
    pub solenoid_length: f64,
    pub solenoid_diameter: f64,
    pub solenoid_turn_count: f64,
    pub anode_voltage: f64,
    pub filament_heating_voltage: f64,
    pub max_anode_current: f64,
    pub nominal_initial_velocity: f64,
    pub solenoid_current: f64,
    pub steps_per_gyroperiod: u32,
}

#[wasm_bindgen]
impl WasmMagnetronConfig {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        let default = magnetron_core::config::MagnetronConfig::default();
        Self {
            anode_radius: default.anode_radius,
            cathode_radius: default.cathode_radius,
            solenoid_length: default.solenoid_length,
            solenoid_diameter: default.solenoid_diameter,
            solenoid_turn_count: default.solenoid_turn_count,
            anode_voltage: default.anode_voltage,
            filament_heating_voltage: default.filament_heating_voltage,
            max_anode_current: default.max_anode_current,
            nominal_initial_velocity: default.nominal_initial_velocity,
            solenoid_current: default.solenoid_current,
            steps_per_gyroperiod: default.steps_per_gyroperiod,
        }
    }
}

impl From<WasmMagnetronConfig> for magnetron_core::config::MagnetronConfig {
    fn from(w: WasmMagnetronConfig) -> Self {
        Self {
            anode_radius: w.anode_radius,
            cathode_radius: w.cathode_radius,
            solenoid_length: w.solenoid_length,
            solenoid_diameter: w.solenoid_diameter,
            solenoid_turn_count: w.solenoid_turn_count,
            anode_voltage: w.anode_voltage,
            filament_heating_voltage: w.filament_heating_voltage,
            max_anode_current: w.max_anode_current,
            nominal_initial_velocity: w.nominal_initial_velocity,
            solenoid_current: w.solenoid_current,
            steps_per_gyroperiod: w.steps_per_gyroperiod,
        }
    }
}

#[wasm_bindgen]
pub struct WasmSimulation {
    config: magnetron_core::config::MagnetronConfig,
    norm: magnetron_core::config::Normalization,
    field: magnetron_core::particles::CoaxialDiodeField,
    phase_space: magnetron_core::particles::PhaseSpace,
    active: Vec<bool>,
    hit_anode: Vec<bool>,
    step_count: usize,
    dt: f64,
    pos_scale: f64,
    charge_m_ratio: f64,
    r_c_norm: f64,
}

#[wasm_bindgen]
impl WasmSimulation {
    #[wasm_bindgen(constructor)]
    pub fn new(config: WasmMagnetronConfig, num_particles: usize) -> Self {
        let config: magnetron_core::config::MagnetronConfig = config.into();
        let norm = config.normalization();
        let ln_ratio = (config.anode_radius / config.cathode_radius).ln();
        
        let e_coeff = -config.anode_voltage / (norm.e_0 * norm.l_0 * ln_ratio);
        let b_normalized_z = 1.0; 
        
        let field = magnetron_core::particles::CoaxialDiodeField {
            e_coeff,
            b_z: b_normalized_z,
        };
        
        let mut phase_space = magnetron_core::particles::PhaseSpace::new(num_particles);
        
        // Seed with a specific value for stability in reproduction
        magnetron_core::diagnostics::initialize_cathode_particles_seeded(&mut phase_space, &config, &norm, 123456);
        
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

    pub fn step(&mut self) -> bool {
        magnetron_core::particles::push_batch_higuera_cary(
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

    pub fn get_active_states(&self) -> Vec<u8> {
        self.active.iter().map(|&a| if a { 1 } else { 0 }).collect()
    }

    pub fn get_anode_current(&self) -> f64 {
        let num_particles = self.phase_space.num_particles();
        let hit_count = self.hit_anode.iter().filter(|&&h| h).count();
        let fraction = (hit_count as f64) / (num_particles as f64);
        fraction * self.config.max_anode_current
    }
    
    pub fn get_anode_radius(&self) -> f64 {
        self.config.anode_radius
    }

    pub fn get_cathode_radius(&self) -> f64 {
        self.config.cathode_radius
    }
}

#[wasm_bindgen]
pub fn run_sweep_js(
    config: WasmMagnetronConfig,
    u_a_vals: Vec<f64>,
    num_ic_points: usize,
    num_particles: usize,
    max_steps: usize,
) -> Result<String, String> {
    let core_config: magnetron_core::config::MagnetronConfig = config.into();
    let result = magnetron_core::diagnostics::run_diagnostics_sweep(
        &core_config,
        &u_a_vals,
        num_ic_points,
        num_particles,
        max_steps,
    );
    serde_json::to_string(&result).map_err(|e| e.to_string())
}
