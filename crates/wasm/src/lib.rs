use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub fn hello() -> f64 {
    // Ensure we reference core to prove the dependency works
    let _info = magnetron_core::get_physics_info();
    42.0
}
