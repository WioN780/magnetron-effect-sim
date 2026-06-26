// Particle Swarm Phase Space and Initialization

use ndarray::{Array2, Zip};
use rand::Rng;
use rand_distr::{Normal, Distribution};
use crate::config::{MagnetronConfig, Normalization};
use crate::pusher::{ElectroMagneticField, push_single_higuera_cary};

/// Particle Phase Space container in Structure of Arrays (SoA) layout.
#[derive(Debug, Clone)]
pub struct PhaseSpace {
    /// Positions: shape (3, N) containing x, y, z coordinates of N particles
    pub positions: Array2<f64>,
    /// Momenta: shape (3, N) containing normalized momenta (u = gamma * v / c) of N particles
    pub momenta: Array2<f64>,
}

impl PhaseSpace {
    /// Creates a new PhaseSpace container for N particles initialized to zero
    pub fn new(n: usize) -> Self {
        Self {
            positions: Array2::zeros((3, n)),
            momenta: Array2::zeros((3, n)),
        }
    }

    /// Returns the number of particles in the phase space
    pub fn num_particles(&self) -> usize {
        self.positions.shape()[1]
    }
}

/// Initialize N particles uniformly distributed along the cathode surface with Maxwell-Boltzmann velocities
pub fn initialize_cathode_particles(
    phase_space: &mut PhaseSpace,
    config: &MagnetronConfig,
    norm: &Normalization,
) {
    let mut rng = rand::thread_rng();
    
    // Normalized cathode radius
    let r_c_norm = norm.normalize_length(config.cathode_radius);
    // Normalized solenoid length (along z)
    let z_max_norm = norm.normalize_length(config.solenoid_length / 2.0);
    let z_min_norm = -z_max_norm;
    
    // Thermal speed normalized by c
    let v_th_norm = norm.normalize_velocity(config.nominal_initial_velocity);
    
    let normal_dist = Normal::new(0.0, v_th_norm).unwrap();
    
    let num_particles = phase_space.num_particles();
    let mut pos_view = phase_space.positions.view_mut();
    let mut mom_view = phase_space.momenta.view_mut();
    
    for i in 0..num_particles {
        // Sample angle theta uniformly in [0, 2*pi)
        let theta = rng.gen_range(0.0..2.0 * std::f64::consts::PI);
        // Sample position z uniformly along the cathode length
        let z = rng.gen_range(z_min_norm..z_max_norm);
        
        let x = r_c_norm * theta.cos();
        let y = r_c_norm * theta.sin();
        
        // Sample local velocities
        // For thermionic emission, radial velocity vr >= 0 is Rayleigh distributed.
        // We sample Rayleigh by taking the norm of two independent 1D Gaussians.
        let vr1 = normal_dist.sample(&mut rng);
        let vr2 = normal_dist.sample(&mut rng);
        let vr = (vr1 * vr1 + vr2 * vr2).sqrt();
        
        let vtheta = normal_dist.sample(&mut rng);
        let vz = normal_dist.sample(&mut rng);
        
        // Transform local velocities (radial, azimuthal, axial) to Cartesian coordinates
        let vx = vr * theta.cos() - vtheta * theta.sin();
        let vy = vr * theta.sin() + vtheta * theta.cos();
        
        // Proper momentum: u = gamma * v
        let v_sq: f64 = vx * vx + vy * vy + vz * vz;
        let gamma = if v_sq >= 1.0 {
            // Safety fallback to prevent division by zero or imaginary numbers
            100.0
        } else {
            1.0 / (1.0 - v_sq).sqrt()
        };
        
        pos_view[[0, i]] = x;
        pos_view[[1, i]] = y;
        pos_view[[2, i]] = z;
        
        mom_view[[0, i]] = gamma * vx;
        mom_view[[1, i]] = gamma * vy;
        mom_view[[2, i]] = gamma * vz;
    }
}

/// Vectorized Higuera-Cary step via ndarray Zip
pub fn push_batch_higuera_cary(
    phase_space: &mut PhaseSpace,
    active: &[bool],
    field: &dyn ElectroMagneticField,
    dt: f64,
    charge_m_ratio: f64,
    pos_scale: f64,
) {
    let mut pos_view = phase_space.positions.view_mut();
    let mut mom_view = phase_space.momenta.view_mut();

    Zip::from(pos_view.columns_mut())
        .and(mom_view.columns_mut())
        .and(ndarray::aview1(active))
        .for_each(|mut pos, mut mom, &is_active| {
            if is_active {
                let mut p_arr = [pos[0], pos[1], pos[2]];
                let mut u_arr = [mom[0], mom[1], mom[2]];
                push_single_higuera_cary(&mut p_arr, &mut u_arr, field, dt, charge_m_ratio, pos_scale);
                pos[0] = p_arr[0];
                pos[1] = p_arr[1];
                pos[2] = p_arr[2];
                mom[0] = u_arr[0];
                mom[1] = u_arr[1];
                mom[2] = u_arr[2];
            }
        });
}

/// Coaxial diode field with radial electric field (1/r) and uniform magnetic field (along z).
#[derive(Debug, Clone)]
pub struct CoaxialDiodeField {
    pub e_coeff: f64,
    pub b_z: f64,
}

impl ElectroMagneticField for CoaxialDiodeField {
    fn evaluate_e(&self, pos: &[f64; 3]) -> [f64; 3] {
        let r = (pos[0] * pos[0] + pos[1] * pos[1]).sqrt();
        if r < 1e-12 {
            [0.0, 0.0, 0.0]
        } else {
            let e_r = self.e_coeff / r;
            [e_r * pos[0] / r, e_r * pos[1] / r, 0.0]
        }
    }

    fn evaluate_b(&self, _pos: &[f64; 3]) -> [f64; 3] {
        [0.0, 0.0, self.b_z]
    }
}

/// Runs a full coaxial diode simulation for a set of particles.
/// Returns the anode fraction (number of particles reaching the anode / total particles).
pub fn run_coaxial_simulation(
    config: &MagnetronConfig,
    num_particles: usize,
    max_steps: usize,
) -> f64 {
    let norm = config.normalization();
    let ln_ratio = (config.anode_radius / config.cathode_radius).ln();
    
    // Calculate normalized fields
    let e_coeff = -config.anode_voltage / (norm.e_0 * norm.l_0 * ln_ratio);
    let b_normalized_z = 1.0; 
    
    let field = CoaxialDiodeField {
        e_coeff,
        b_z: b_normalized_z,
    };
    
    let mut ps = PhaseSpace::new(num_particles);
    initialize_cathode_particles(&mut ps, config, &norm);
    
    let mut active = vec![true; num_particles];
    let mut hit_anode = vec![false; num_particles];
    
    let r_c_norm = norm.normalize_length(config.cathode_radius);
    let charge_m_ratio = -1.0;
    
    // Compute normalized dt: gyroperiod is 2pi.
    let dt = 2.0 * std::f64::consts::PI / (config.steps_per_gyroperiod as f64);
    let pos_scale = norm.v_0 * norm.t_0 / norm.l_0;
    
    for step in 0..max_steps {
        // Push active particles
        push_batch_higuera_cary(&mut ps, &active, &field, dt, charge_m_ratio, pos_scale);
        
        let mut active_count = 0;
        let pos_view = ps.positions.view();
        
        for i in 0..num_particles {
            if active[i] {
                let x = pos_view[[0, i]];
                let y = pos_view[[1, i]];
                let r = (x*x + y*y).sqrt();
                
                if r >= 1.0 {
                    // Hit anode
                    hit_anode[i] = true;
                    active[i] = false;
                } else if r <= r_c_norm && step > 2 {
                    // Returned to cathode
                    active[i] = false;
                } else {
                    active_count += 1;
                }
            }
        }
        
        if active_count == 0 {
            break;
        }
    }
    
    // Anode fraction
    let hit_count = hit_anode.iter().filter(|&&h| h).count();
    (hit_count as f64) / (num_particles as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cathode_initialization() {
        let config = MagnetronConfig::default();
        let norm = config.normalization();
        
        let n = 1000;
        let mut ps = PhaseSpace::new(n);
        initialize_cathode_particles(&mut ps, &config, &norm);
        
        let r_c_norm = norm.normalize_length(config.cathode_radius);
        
        // Assert all initialized particles are on the cathode surface
        for i in 0..n {
            let x = ps.positions[[0, i]];
            let y = ps.positions[[1, i]];
            let r = (x*x + y*y).sqrt();
            assert!((r - r_c_norm).abs() < 1e-12, "Particle {} not on cathode surface! r={}, r_c={}", i, r, r_c_norm);
        }
    }
}
