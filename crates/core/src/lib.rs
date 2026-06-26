// Magnetron Core Physics Engine Library

pub mod constants;
pub mod config;
pub mod pusher;
pub mod particles;
pub mod diagnostics;
pub mod poisson;

pub fn get_physics_info() -> &'static str {
    "Magnetron Core Physics Engine"
}
