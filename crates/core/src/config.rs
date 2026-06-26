// Magnetron Configuration and Normalization Definitions

use crate::constants::{C, E, M_E, MU_0};

/// Configuration for the Magnetron apparatus (based on 2D2S vacuum diode).
#[derive(Debug, Clone)]
pub struct MagnetronConfig {
    /// Anode radius (m)
    pub anode_radius: f64,
    /// Cathode (filament) radius (m)
    pub cathode_radius: f64,
    /// Solenoid length (m)
    pub solenoid_length: f64,
    /// Solenoid diameter (m)
    pub solenoid_diameter: f64,
    /// Solenoid turns count
    pub solenoid_turn_count: f64,
    /// Anode operating voltage (V)
    pub anode_voltage: f64,
    /// Filament heating voltage (V)
    pub filament_heating_voltage: f64,
    /// Maximum anode current (A)
    pub max_anode_current: f64,
    /// Nominal initial thermal velocity (m/s)
    pub nominal_initial_velocity: f64,
    /// Solenoid current (A)
    pub solenoid_current: f64,
    /// Number of integration steps per local cyclotron gyroperiod
    pub steps_per_gyroperiod: u32,
}

impl Default for MagnetronConfig {
    fn default() -> Self {
        Self {
            anode_radius: 5.1e-3,            // 5.1 mm
            cathode_radius: 6.25e-5,         // 62.5 um
            solenoid_length: 0.167,          // 16.7 cm
            solenoid_diameter: 0.062,        // 6.2 cm
            solenoid_turn_count: 2300.0,
            anode_voltage: 40.0,
            filament_heating_voltage: 1.5,   // ~1.2 - 1.7 V
            max_anode_current: 0.040,        // 40 mA
            nominal_initial_velocity: 1.1e6, // 1.1 * 10^6 m/s
            solenoid_current: 1.0,           // 1.0 A nominal
            steps_per_gyroperiod: 32,
        }
    }
}

impl MagnetronConfig {
    /// Computes the central homogeneous magnetic field of the solenoid (T)
    pub fn central_b_field(&self) -> f64 {
        let l = self.solenoid_length;
        let d = self.solenoid_diameter;
        (MU_0 * self.solenoid_turn_count * self.solenoid_current) / (l * l + d * d).sqrt()
    }

    /// Computes the characteristic cyclotron frequency at the center (rad/s)
    pub fn characteristic_cyclotron_frequency(&self) -> f64 {
        let b = self.central_b_field();
        (E * b) / M_E
    }

    /// Returns a new Normalization helper based on this configuration
    pub fn normalization(&self) -> Normalization {
        Normalization::new(self.anode_radius, self.central_b_field())
    }
}

/// Normalization helper to map physical (SI) quantities to dimensionless values.
///
/// Length is normalized by anode radius R_a
/// Time is normalized by characteristic cyclotron period (1 / omega_c0)
/// Velocity is normalized by the speed of light c
#[derive(Debug, Clone)]
pub struct Normalization {
    /// Length scale factor (R_a)
    pub l_0: f64,
    /// Time scale factor (1 / omega_c)
    pub t_0: f64,
    /// Velocity scale factor (c)
    pub v_0: f64,
    /// Magnetic field scale factor (B_0)
    pub b_0: f64,
    /// Electric field scale factor (c * B_0)
    pub e_0: f64,
}

impl Normalization {
    /// Creates a new normalization system given length scale R_a and magnetic field scale B_0
    pub fn new(r_a: f64, b_0: f64) -> Self {
        // Handle potential zero magnetic field safely by falling back to 1.0 T scale
        let b_scale = if b_0.abs() < 1e-12 { 1.0 } else { b_0 };
        let omega_c0 = (E * b_scale) / M_E;
        Self {
            l_0: r_a,
            t_0: 1.0 / omega_c0,
            v_0: C,
            b_0: b_scale,
            e_0: C * b_scale,
        }
    }

    pub fn normalize_length(&self, x: f64) -> f64 { x / self.l_0 }
    pub fn denormalize_length(&self, x: f64) -> f64 { x * self.l_0 }

    pub fn normalize_pos(&self, pos: [f64; 3]) -> [f64; 3] {
        [pos[0] / self.l_0, pos[1] / self.l_0, pos[2] / self.l_0]
    }
    pub fn denormalize_pos(&self, pos: [f64; 3]) -> [f64; 3] {
        [pos[0] * self.l_0, pos[1] * self.l_0, pos[2] * self.l_0]
    }

    pub fn normalize_time(&self, t: f64) -> f64 { t / self.t_0 }
    pub fn denormalize_time(&self, t: f64) -> f64 { t * self.t_0 }

    pub fn normalize_velocity(&self, v: f64) -> f64 { v / self.v_0 }
    pub fn denormalize_velocity(&self, v: f64) -> f64 { v * self.v_0 }

    pub fn normalize_vel(&self, vel: [f64; 3]) -> [f64; 3] {
        [vel[0] / self.v_0, vel[1] / self.v_0, vel[2] / self.v_0]
    }
    pub fn denormalize_vel(&self, vel: [f64; 3]) -> [f64; 3] {
        [vel[0] * self.v_0, vel[1] * self.v_0, vel[2] * self.v_0]
    }

    pub fn normalize_b(&self, b: [f64; 3]) -> [f64; 3] {
        [b[0] / self.b_0, b[1] / self.b_0, b[2] / self.b_0]
    }
    pub fn denormalize_b(&self, b: [f64; 3]) -> [f64; 3] {
        [b[0] * self.b_0, b[1] * self.b_0, b[2] * self.b_0]
    }

    pub fn normalize_e(&self, e: [f64; 3]) -> [f64; 3] {
        [e[0] / self.e_0, e[1] / self.e_0, e[2] / self.e_0]
    }
    pub fn denormalize_e(&self, e: [f64; 3]) -> [f64; 3] {
        [e[0] * self.e_0, e[1] * self.e_0, e[2] * self.e_0]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normalization_round_trip() {
        let config = MagnetronConfig::default();
        let norm = config.normalization();

        let length = 0.00345;
        let norm_len = norm.normalize_length(length);
        let back_len = norm.denormalize_length(norm_len);
        assert!((length - back_len).abs() < 1e-14);

        let pos = [0.001, -0.002, 0.003];
        let norm_pos = norm.normalize_pos(pos);
        let back_pos = norm.denormalize_pos(norm_pos);
        for i in 0..3 {
            assert!((pos[i] - back_pos[i]).abs() < 1e-14);
        }

        let time = 1.234e-8;
        let norm_time = norm.normalize_time(time);
        let back_time = norm.denormalize_time(norm_time);
        assert!((time - back_time).abs() < 1e-14);

        let vel = [1e6, -2e6, 3e6];
        let norm_vel = norm.normalize_vel(vel);
        let back_vel = norm.denormalize_vel(norm_vel);
        for i in 0..3 {
            assert!((vel[i] - back_vel[i]).abs() < 1e-14);
        }

        let b = [0.01, -0.02, 0.03];
        let norm_b = norm.normalize_b(b);
        let back_b = norm.denormalize_b(norm_b);
        for i in 0..3 {
            assert!((b[i] - back_b[i]).abs() < 1e-14);
        }

        let e = [100.0, -200.0, 300.0];
        let norm_e = norm.normalize_e(e);
        let back_e = norm.denormalize_e(norm_e);
        for i in 0..3 {
            assert!((e[i] - back_e[i]).abs() < 1e-14);
        }
    }
}
