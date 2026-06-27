#![allow(non_snake_case)]

use dioxus::prelude::*;
use dioxus_web::WebEventExt;
use wasm_bindgen::prelude::*;
use wasm_bindgen::JsCast;
use web_sys::{CanvasRenderingContext2d, HtmlCanvasElement};
use magnetron_core::{Simulation, MagnetronConfig};
use std::rc::Rc;
use std::cell::RefCell;
use std::future::Future;
use std::pin::Pin;
use std::task::{Context, Poll};

// Include the stylesheet using the Dioxus 0.7 asset system
const MY_STYLES: Asset = asset!("/assets/main.css");

/// Future that resolves on the next browser animation frame (requestAnimationFrame).
pub struct NextFrame {
    state: Rc<RefCell<NextFrameState>>,
}

struct NextFrameState {
    waker: Option<std::task::Waker>,
    fired: bool,
    timestamp: f64,
}

impl Future for NextFrame {
    type Output = f64;
    fn poll(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Self::Output> {
        let mut state = self.state.borrow_mut();
        if state.fired {
            Poll::Ready(state.timestamp)
        } else {
            state.waker = Some(cx.waker().clone());
            Poll::Pending
        }
    }
}

pub fn next_frame() -> NextFrame {
    let state = Rc::new(RefCell::new(NextFrameState {
        waker: None,
        fired: false,
        timestamp: 0.0,
    }));
    
    let state_clone = state.clone();
    let closure = Closure::once(move |timestamp: f64| {
        let mut state = state_clone.borrow_mut();
        state.fired = true;
        state.timestamp = timestamp;
        if let Some(waker) = state.waker.take() {
            waker.wake();
        }
    });
    
    web_sys::window()
        .unwrap()
        .request_animation_frame(closure.as_ref().unchecked_ref())
        .unwrap();
        
    closure.forget();
    
    NextFrame { state }
}

fn main() {
    dioxus::launch(App);
}

fn App() -> Element {
    // Dynamic Geometry Config (Autoreset on change)
    let mut anode_radius = use_signal(|| 5.1); // in mm (SI: 5.1e-3 m)
    let mut cathode_radius = use_signal(|| 0.0625); // in mm (SI: 62.5e-6 m)
    let mut solenoid_turns = use_signal(|| 2300.0);
    let mut initial_velocity = use_signal(|| 1.1); // in million m/s (SI: 1.1e6 m/s)
    let mut num_particles = use_signal(|| 150);

    // Live Config (Dynamic update without reset)
    let mut anode_voltage = use_signal(|| 40.0); // in V
    let mut solenoid_current = use_signal(|| 1.0); // in A
    let mut filament_heating_voltage = use_signal(|| 1.5); // in V

    // Simulation Visual Factors
    let mut time_scale = use_signal(|| 1.0); // 0.1 to 3.0 (time dilation)
    let mut trace_length = use_signal(|| 15); // 0 to 40 step histories

    // State Playback Controls
    let mut is_paused = use_signal(|| false);
    let mut reset_trigger = use_signal(|| 0);

    // Simulation Stats Signals
    let mut step_count = use_signal(|| 0);
    let mut active_count = use_signal(|| 0);
    let mut anode_current = use_signal(|| 0.0);

    // Canvas Reference Signal
    let mut canvas_ref = use_signal(|| None::<HtmlCanvasElement>);

    // Compute live physics parameters for display and critical calculations
    let current_config = MagnetronConfig {
        anode_radius: anode_radius() / 1000.0,
        cathode_radius: cathode_radius() / 1000.0,
        solenoid_turn_count: solenoid_turns(),
        nominal_initial_velocity: initial_velocity() * 1e6,
        anode_voltage: anode_voltage(),
        solenoid_current: solenoid_current(),
        filament_heating_voltage: filament_heating_voltage(),
        ..MagnetronConfig::default()
    };
    
    // Homogeneous central field B_0
    let b_field = current_config.central_b_field();

    // Hull critical cutoff magnetic field:
    // B_c = sqrt(8 * m_e * V_a / (e * R_a^2 * (1 - R_c^2 / R_a^2)^2))
    let ra_val = current_config.anode_radius;
    let rc_val = current_config.cathode_radius;
    let denom = ra_val * ra_val * (1.0 - (rc_val * rc_val) / (ra_val * ra_val)).powi(2);
    let b_cutoff = if denom > 0.0 && current_config.anode_voltage >= 0.0 {
        ((8.0 * magnetron_core::constants::M_E * current_config.anode_voltage) 
         / (magnetron_core::constants::E * denom)).sqrt()
    } else {
        0.0
    };

    let is_cutoff = b_field > b_cutoff;

    // Persistent simulation loop running on requestAnimationFrame
    let _sim_loop = use_future(move || async move {
        // Track the previous states of structural parameters to auto-reset simulation
        let mut last_anode_radius = *anode_radius.peek();
        let mut last_cathode_radius = *cathode_radius.peek();
        let mut last_solenoid_turns = *solenoid_turns.peek();
        let mut last_initial_velocity = *initial_velocity.peek();
        let mut last_num_particles = *num_particles.peek();
        let mut last_reset = *reset_trigger.peek();

        // Instantiate simulation with initial parameters
        let mut sim = Simulation::new(
            MagnetronConfig {
                anode_radius: last_anode_radius / 1000.0,
                cathode_radius: last_cathode_radius / 1000.0,
                solenoid_turn_count: last_solenoid_turns,
                nominal_initial_velocity: last_initial_velocity * 1e6,
                anode_voltage: anode_voltage.peek().clone(),
                solenoid_current: solenoid_current.peek().clone(),
                filament_heating_voltage: filament_heating_voltage.peek().clone(),
                ..MagnetronConfig::default()
            },
            last_num_particles
        );
        let mut particle_history = vec![std::collections::VecDeque::new(); last_num_particles];

        let mut last_timestamp = 0.0;
        let mut accumulator = 0.0;
        let mut last_physics_positions = sim.get_positions();
        let mut current_physics_positions = sim.get_positions();

        loop {
            // Wait until next frame paint and read timestamp
            let timestamp = next_frame().await;

            // Check if any structural parameters changed or if reset button was clicked
            let cur_anode_radius = *anode_radius.read();
            let cur_cathode_radius = *cathode_radius.read();
            let cur_solenoid_turns = *solenoid_turns.read();
            let cur_initial_velocity = *initial_velocity.read();
            let cur_num_particles = *num_particles.read();
            let cur_reset = *reset_trigger.read();

            let should_reset = cur_reset != last_reset
                || cur_anode_radius != last_anode_radius
                || cur_cathode_radius != last_cathode_radius
                || cur_solenoid_turns != last_solenoid_turns
                || cur_initial_velocity != last_initial_velocity
                || cur_num_particles != last_num_particles;

            if should_reset {
                last_reset = cur_reset;
                last_anode_radius = cur_anode_radius;
                last_cathode_radius = cur_cathode_radius;
                last_solenoid_turns = cur_solenoid_turns;
                last_initial_velocity = cur_initial_velocity;
                last_num_particles = cur_num_particles;

                sim = Simulation::new(
                    MagnetronConfig {
                        anode_radius: cur_anode_radius / 1000.0,
                        cathode_radius: cur_cathode_radius / 1000.0,
                        solenoid_turn_count: cur_solenoid_turns,
                        nominal_initial_velocity: cur_initial_velocity * 1e6,
                        anode_voltage: *anode_voltage.read(),
                        solenoid_current: *solenoid_current.read(),
                        filament_heating_voltage: *filament_heating_voltage.read(),
                        ..MagnetronConfig::default()
                    },
                    cur_num_particles
                );
                particle_history = vec![std::collections::VecDeque::new(); cur_num_particles];
                last_physics_positions = sim.get_positions();
                current_physics_positions = sim.get_positions();
                accumulator = 0.0;
            }

            // Sync dynamic slider updates to running simulation physics parameters
            sim.config.anode_voltage = *anode_voltage.read();
            sim.config.solenoid_current = *solenoid_current.read();
            sim.config.filament_heating_voltage = *filament_heating_voltage.read();
            
            let norm = sim.config.normalization();
            let ln_ratio = (sim.config.anode_radius / sim.config.cathode_radius).ln();
            let e_coeff = -sim.config.anode_voltage / (norm.e_0 * norm.l_0 * ln_ratio);
            sim.norm = norm;
            sim.field.e_coeff = e_coeff;
            
            // FIXED_DT stays exactly the value validated (never scaled by the speed control)
            let fixed_dt = 2.0 * std::f64::consts::PI / (sim.config.steps_per_gyroperiod as f64);
            sim.dt = fixed_dt;
            sim.pos_scale = sim.norm.v_0 * sim.norm.t_0 / sim.norm.l_0;

            // Calculate elapsed frame time in seconds
            let real_dt_seconds = if last_timestamp > 0.0 {
                (timestamp - last_timestamp) / 1000.0
            } else {
                1.0 / 60.0
            };
            last_timestamp = timestamp;

            // Cap elapsed time to prevent physics explosions on background suspension
            let real_dt_seconds = real_dt_seconds.min(0.1);

            let time_scale_val = *time_scale.read();

            // Run simulation steps if not paused
            if !*is_paused.read() {
                // Accumulate normalized simulation time: 1.0 timeScale at 60fps = 1 fixed_dt per frame
                accumulator += real_dt_seconds * 60.0 * fixed_dt * time_scale_val;

                let mut steps_run = 0;
                let n_parts = sim.phase_space.num_particles();

                // Ensure history structure is properly aligned
                if particle_history.len() != n_parts {
                    particle_history = vec![std::collections::VecDeque::new(); n_parts];
                }

                while accumulator >= fixed_dt && steps_run < 240 {
                    last_physics_positions = current_physics_positions.clone();
                    sim.step();
                    current_physics_positions = sim.get_positions();
                    accumulator -= fixed_dt;
                    steps_run += 1;

                    // Push intermediate step to history for smooth traces
                    let active = sim.get_active_states();
                    for i in 0..n_parts {
                        if active[i] {
                            let x = current_physics_positions[i * 3];
                            let y = current_physics_positions[i * 3 + 1];
                            let history = &mut particle_history[i];
                            history.push_back([x, y]);
                            
                            let limit = *trace_length.read();
                            while history.len() > limit {
                                history.pop_front();
                            }
                        } else {
                            let history = &mut particle_history[i];
                            if !history.is_empty() {
                                history.pop_front();
                            }
                        }
                    }
                }

                // If loop cap hit, discard remaining accumulator
                if steps_run >= 240 {
                    accumulator = 0.0;
                }

                // Update UI Stats Signals
                let active = sim.get_active_states();
                step_count.set(sim.step_count);
                active_count.set(active.iter().filter(|&&a| a).count());
                anode_current.set(sim.get_anode_current());
            } else {
                // Clear accumulator during pause to prevent visual jumps on resume
                accumulator = 0.0;
            }

            // Interpolate draw positions for rendering
            let blend = (accumulator / fixed_dt).min(1.0).max(0.0);
            let mut draw_positions = vec![0.0; current_physics_positions.len()];
            for i in 0..current_physics_positions.len() {
                if i < last_physics_positions.len() && i < current_physics_positions.len() {
                    draw_positions[i] = last_physics_positions[i] * (1.0 - blend) + current_physics_positions[i] * blend;
                }
            }

            // Draw to Canvas
            if let Some(canvas) = canvas_ref.read().as_ref() {
                if let Ok(Some(ctx_val)) = canvas.get_context("2d") {
                    if let Ok(ctx) = ctx_val.dyn_into::<CanvasRenderingContext2d>() {
                        let width = canvas.width() as f64;
                        let height = canvas.height() as f64;
                        let cx = width / 2.0;
                        let cy = height / 2.0;
                        let r_draw = cx.min(cy) * 0.94;
                        
                        // Clear canvas viewport
                        ctx.clear_rect(0.0, 0.0, width, height);

                        // 1. Draw Anode Ring (Neutral Gray)
                        ctx.begin_path();
                        ctx.arc(cx, cy, r_draw, 0.0, std::f64::consts::PI * 2.0).unwrap();
                        ctx.set_fill_style_str("#FFFFFF");
                        ctx.fill();
                        ctx.set_stroke_style_str("#CBD5E1");
                        ctx.set_line_width(4.0);
                        ctx.stroke();

                        // 2. Draw Cathode Filament (Soft Red, heated)
                        let r_c_draw = (sim.config.cathode_radius / sim.config.anode_radius) * r_draw;
                        let r_c_draw = r_c_draw.max(5.0); // Minimum visible radius
                        ctx.begin_path();
                        ctx.arc(cx, cy, r_c_draw, 0.0, std::f64::consts::PI * 2.0).unwrap();
                        ctx.set_fill_style_str("#FCA5A5");
                        ctx.fill();
                        ctx.set_stroke_style_str("#E2E8F0");
                        ctx.set_line_width(1.0);
                        ctx.stroke();

                        // 3. Draw Electron Traces (Fading pastel blue: #93C5FD)
                        let history_limit = *trace_length.read();
                        if history_limit > 0 {
                            ctx.begin_path();
                            for i in 0..sim.phase_space.num_particles() {
                                if i >= particle_history.len() {
                                    continue;
                                }
                                let history = &particle_history[i];
                                if history.len() < 2 {
                                    continue;
                                }
                                let p0 = history[0];
                                let px = cx + (p0[0] / sim.config.anode_radius) * r_draw;
                                let py = cy + (p0[1] / sim.config.anode_radius) * r_draw;
                                ctx.move_to(px, py);
                                for j in 1..history.len() {
                                    let pj = history[j];
                                    let pxj = cx + (pj[0] / sim.config.anode_radius) * r_draw;
                                    let pyj = cy + (pj[1] / sim.config.anode_radius) * r_draw;
                                    ctx.line_to(pxj, pyj);
                                }
                            }
                            ctx.set_stroke_style_str("rgba(147, 197, 253, 0.55)");
                            ctx.set_line_width(1.5);
                            ctx.stroke();
                        }

                        // 4. Draw Electron Heads (Deep royal blue: #2563EB)
                        ctx.begin_path();
                        let active = sim.get_active_states();
                        for i in 0..sim.phase_space.num_particles() {
                            if i < active.len() && active[i] && (i * 3 + 1) < draw_positions.len() {
                                let x = draw_positions[i * 3];
                                let y = draw_positions[i * 3 + 1];
                                let px = cx + (x / sim.config.anode_radius) * r_draw;
                                let py = cy + (y / sim.config.anode_radius) * r_draw;
                                ctx.move_to(px + 2.5, py);
                                ctx.arc(px, py, 2.5, 0.0, std::f64::consts::PI * 2.0).unwrap();
                            }
                        }
                        ctx.set_fill_style_str("#2563EB");
                        ctx.fill();
                    }
                }
            }
        }
    });

    rsx! {
        // Enforce CSS injection into Document Head
        document::Link { rel: "stylesheet", href: MY_STYLES }

        div { class: "app-container",
            // Controls sidebar panel
            div { class: "sidebar",
                h1 { class: "sidebar-title", "Magnetron Sim" }
                
                div { class: "sidebar-section",
                    // 1. DIODE GEOMETRY CONFIG
                    div { class: "sidebar-group",
                        h3 { class: "sidebar-group-title", "Diode Geometry" }
                        
                        // Anode Radius Slider
                        div { class: "control-group",
                            label { class: "control-label",
                                span { "Anode Radius (R_a)" }
                                span { class: "control-value", "{anode_radius:.1} mm" }
                            }
                            input {
                                class: "control-input",
                                r#type: "range",
                                min: "2.0",
                                max: "12.0",
                                step: "0.1",
                                value: "{anode_radius}",
                                oninput: move |e| {
                                    if let Ok(val) = e.value().parse::<f64>() {
                                        anode_radius.set(val);
                                    }
                                }
                            }
                        }

                        // Cathode Radius Slider
                        div { class: "control-group",
                            label { class: "control-label",
                                span { "Cathode Radius (R_c)" }
                                span { class: "control-value", "{cathode_radius * 1000.0:.0} um" }
                            }
                            input {
                                class: "control-input",
                                r#type: "range",
                                min: "0.02",
                                max: "0.50",
                                step: "0.01",
                                value: "{cathode_radius}",
                                oninput: move |e| {
                                    if let Ok(val) = e.value().parse::<f64>() {
                                        cathode_radius.set(val);
                                    }
                                }
                            }
                        }

                        // Solenoid Turns
                        div { class: "control-group",
                            label { class: "control-label",
                                span { "Solenoid Turns (N)" }
                                span { class: "control-value", "{solenoid_turns}" }
                            }
                            input {
                                class: "control-input",
                                r#type: "range",
                                min: "500",
                                max: "5000",
                                step: "100",
                                value: "{solenoid_turns}",
                                oninput: move |e| {
                                    if let Ok(val) = e.value().parse::<f64>() {
                                        solenoid_turns.set(val);
                                    }
                                }
                            }
                        }
                    }

                    // 2. ELECTROMAGNETIC CONFIG
                    div { class: "sidebar-group",
                        h3 { class: "sidebar-group-title", "Electromagnetics" }
                        
                        // Anode Voltage Control
                        div { class: "control-group",
                            label { class: "control-label",
                                span { "Anode Voltage (V_a)" }
                                span { class: "control-value", "{anode_voltage} V" }
                            }
                            input {
                                class: "control-input",
                                r#type: "range",
                                min: "0",
                                max: "150",
                                step: "1",
                                value: "{anode_voltage}",
                                oninput: move |e| {
                                    if let Ok(val) = e.value().parse::<f64>() {
                                        anode_voltage.set(val);
                                    }
                                }
                            }
                        }

                        // Solenoid Current Control
                        div { class: "control-group",
                            label { class: "control-label",
                                span { "Solenoid Current (I_s)" }
                                span { class: "control-value", "{solenoid_current:.2} A" }
                            }
                            input {
                                class: "control-input",
                                r#type: "range",
                                min: "0.0",
                                max: "4.0",
                                step: "0.05",
                                value: "{solenoid_current}",
                                oninput: move |e| {
                                    if let Ok(val) = e.value().parse::<f64>() {
                                        solenoid_current.set(val);
                                    }
                                }
                            }
                        }

                        // Filament Heating Voltage Control
                        div { class: "control-group",
                            label { class: "control-label",
                                span { "Filament Voltage" }
                                span { class: "control-value", "{filament_heating_voltage:.2} V" }
                            }
                            input {
                                class: "control-input",
                                r#type: "range",
                                min: "0.1",
                                max: "3.0",
                                step: "0.05",
                                value: "{filament_heating_voltage}",
                                oninput: move |e| {
                                    if let Ok(val) = e.value().parse::<f64>() {
                                        filament_heating_voltage.set(val);
                                    }
                                }
                            }
                        }
                    }

                    // 3. EMISSION SWARM CONFIG
                    div { class: "sidebar-group",
                        h3 { class: "sidebar-group-title", "Emission Swarm" }
                        
                        // Particles / Electrons Count Control
                        div { class: "control-group",
                            label { class: "control-label",
                                span { "Electrons Count" }
                                span { class: "control-value", "{num_particles}" }
                            }
                            input {
                                class: "control-input",
                                r#type: "range",
                                min: "10",
                                max: "500",
                                step: "10",
                                value: "{num_particles}",
                                oninput: move |e| {
                                    if let Ok(val) = e.value().parse::<usize>() {
                                        num_particles.set(val);
                                    }
                                }
                            }
                        }

                        // Initial Thermal Velocity Slider
                        div { class: "control-group",
                            label { class: "control-label",
                                span { "Initial Vel (v_th)" }
                                span { class: "control-value", "{initial_velocity:.1}M m/s" }
                            }
                            input {
                                class: "control-input",
                                r#type: "range",
                                min: "0.1",
                                max: "3.0",
                                step: "0.1",
                                value: "{initial_velocity}",
                                oninput: move |e| {
                                    if let Ok(val) = e.value().parse::<f64>() {
                                        initial_velocity.set(val);
                                    }
                                }
                            }
                        }
                    }

                    // 4. SIMULATION VISUAL CONFIG
                    div { class: "sidebar-group",
                        h3 { class: "sidebar-group-title", "Simulation & Visuals" }

                        // Time Scale Slider
                        div { class: "control-group",
                            label { class: "control-label",
                                span { "Time Flow Rate" }
                                span { class: "control-value", "{time_scale:.2}x" }
                            }
                            input {
                                class: "control-input",
                                r#type: "range",
                                min: "0.05",
                                max: "3.00",
                                step: "0.05",
                                value: "{time_scale}",
                                oninput: move |e| {
                                    if let Ok(val) = e.value().parse::<f64>() {
                                        time_scale.set(val);
                                    }
                                }
                            }
                        }

                        // Electron Traces Length
                        div { class: "control-group",
                            label { class: "control-label",
                                span { "Trace Tail Length" }
                                span { class: "control-value", "{trace_length}" }
                            }
                            input {
                                class: "control-input",
                                r#type: "range",
                                min: "0",
                                max: "40",
                                step: "1",
                                value: "{trace_length}",
                                oninput: move |e| {
                                    if let Ok(val) = e.value().parse::<usize>() {
                                        trace_length.set(val);
                                    }
                                }
                            }
                        }
                    }

                    // Simulation State Buttons
                    div { class: "btn-container",
                        button {
                            class: "btn",
                            onclick: move |_| {
                                is_paused.toggle();
                            },
                            if is_paused() { "Resume" } else { "Pause" }
                        }
                        button {
                            class: "btn",
                            onclick: move |_| {
                                reset_trigger.set(reset_trigger() + 1);
                            },
                            "Reset"
                        }
                    }
                }
            }

            // Main rendering view
            div { class: "main-view",
                div { class: "canvas-container",
                    canvas {
                        id: "magnetron-canvas",
                        width: "600",
                        height: "600",
                        onmounted: move |evt| {
                            let web_sys_element = evt.as_web_event();
                            if let Ok(canvas) = web_sys_element.dyn_into::<HtmlCanvasElement>() {
                                canvas_ref.set(Some(canvas));
                            }
                        }
                    }
                }

                // Statistics panels
                div { class: "stats-panel",
                    div { class: "stat-item",
                        span { class: "stat-label", "Step Count" }
                        span { class: "stat-value", "{step_count}" }
                    }
                    div { class: "stat-item",
                        span { class: "stat-label", "Active Electrons" }
                        span { class: "stat-value", "{active_count} / {num_particles}" }
                    }
                    div { class: "stat-item",
                        span { class: "stat-label", "Anode Current" }
                        span { class: "stat-value", "{anode_current * 1000.0:.2} mA" }
                    }
                    div { class: "stat-item",
                        span { class: "stat-label", "Magnetic Field" }
                        span { class: "stat-value", "{b_field * 1000.0:.1} mT" }
                        if is_cutoff {
                            span { class: "badge badge-cutoff", "Hull Cutoff" }
                        } else {
                            span { class: "badge badge-conducting", "Conducting" }
                        }
                    }
                    div { class: "stat-item",
                        span { class: "stat-label", "Cutoff Limit (B_c)" }
                        span { class: "stat-value", "{b_cutoff * 1000.0:.1} mT" }
                    }
                }
            }
        }
    }
}
