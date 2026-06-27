use ndarray::Array2;
use std::f64::consts::PI;

/// Cylindrical coordinate (r, z) grid definition.
/// Radial spacing is exponential (clustered near cathode), axial spacing is uniform.
#[derive(Debug, Clone)]
pub struct PoissonGrid {
    pub r: Vec<f64>,
    pub z: Vec<f64>,
    pub nr: usize,
    pub nz: usize,
    pub r_c: f64,
    pub r_a: f64,
    pub z_min: f64,
    pub z_max: f64,
    pub delta_s: f64, // Logical spacing in r: ln(r_a / r_c) / nr
    pub delta_z: f64, // Spacing in z
}

impl PoissonGrid {
    /// Creates a new cylindrical grid with exponential radial spacing and uniform axial spacing.
    pub fn new(r_c: f64, r_a: f64, z_min: f64, z_max: f64, nr: usize, nz: usize) -> Self {
        let mut r = vec![0.0; nr + 1];
        let mut z = vec![0.0; nz + 1];

        let delta_s = (r_a / r_c).ln() / (nr as f64);
        for i in 0..=nr {
            r[i] = r_c * (delta_s * i as f64).exp();
        }

        let delta_z = (z_max - z_min) / (nz as f64);
        for j in 0..=nz {
            z[j] = z_min + delta_z * j as f64;
        }

        Self {
            r,
            z,
            nr,
            nz,
            r_c,
            r_a,
            z_min,
            z_max,
            delta_s,
            delta_z,
        }
    }

    /// Maps physical coordinates (r, z) to logical coordinate space (u, v).
    /// Clamps coordinates to grid boundaries.
    pub fn map_to_logical(&self, r_p: f64, z_p: f64) -> (f64, f64) {
        let r_clamped = r_p.clamp(self.r_c, self.r_a);
        let z_clamped = z_p.clamp(self.z_min, self.z_max);

        let u = (r_clamped / self.r_c).ln() / self.delta_s;
        let v = (z_clamped - self.z_min) / self.delta_z;

        (u, v)
    }

    /// Computes the control cell volume associated with node (i, j).
    pub fn cell_volume(&self, i: usize, j: usize) -> f64 {
        // Cell width in r: dr_i = (r_{i+1} - r_{i-1}) / 2
        let dr = if i == 0 {
            (self.r[1] - self.r[0]) / 2.0
        } else if i == self.nr {
            (self.r[self.nr] - self.r[self.nr - 1]) / 2.0
        } else {
            (self.r[i + 1] - self.r[i - 1]) / 2.0
        };

        // Cell width in z: dz_j = h_z for interior, h_z / 2 for boundaries
        let dz = if j == 0 || j == self.nz {
            self.delta_z / 2.0
        } else {
            self.delta_z
        };

        // Control volume V = 2 * pi * r_i * dr_i * dz_j
        2.0 * PI * self.r[i] * dr * dz
    }

    /// Deposit charge from particle positions using Cloud-In-Cell (CIC) bilinear interpolation.
    pub fn deposit_cic(&self, positions: &Array2<f64>, q_p: &[f64]) -> Array2<f64> {
        let mut q_grid = Array2::<f64>::zeros((self.nr + 1, self.nz + 1));
        let num_particles = positions.shape()[1];

        for p in 0..num_particles {
            let x = positions[[0, p]];
            let y = positions[[1, p]];
            let z_p = positions[[2, p]];
            let r_p = (x * x + y * y).sqrt();

            let (u, v) = self.map_to_logical(r_p, z_p);

            let i = (u.floor() as usize).min(self.nr - 1);
            let j = (v.floor() as usize).min(self.nz - 1);

            let du = u - i as f64;
            let dv = v - j as f64;

            let w00 = (1.0 - du) * (1.0 - dv);
            let w10 = du * (1.0 - dv);
            let w01 = (1.0 - du) * dv;
            let w11 = du * dv;

            q_grid[[i, j]] += q_p[p] * w00;
            q_grid[[i + 1, j]] += q_p[p] * w10;
            q_grid[[i, j + 1]] += q_p[p] * w01;
            q_grid[[i + 1, j + 1]] += q_p[p] * w11;
        }

        // Convert charge to density: rho_ij = q_grid_ij / V_ij
        let mut rho = Array2::<f64>::zeros((self.nr + 1, self.nz + 1));
        for i in 0..=self.nr {
            for j in 0..=self.nz {
                rho[[i, j]] = q_grid[[i, j]] / self.cell_volume(i, j);
            }
        }

        rho
    }

    /// Deposit charge from particle positions using Triangular-Shaped-Cloud (TSC) quadratic spline.
    pub fn deposit_tsc(&self, positions: &Array2<f64>, q_p: &[f64]) -> Array2<f64> {
        let mut q_grid = Array2::<f64>::zeros((self.nr + 1, self.nz + 1));
        let num_particles = positions.shape()[1];

        for p in 0..num_particles {
            let x = positions[[0, p]];
            let y = positions[[1, p]];
            let z_p = positions[[2, p]];
            let r_p = (x * x + y * y).sqrt();

            let (u, v) = self.map_to_logical(r_p, z_p);

            let i = u.round() as i32;
            let j = v.round() as i32;

            let du = u - i as f64;
            let dv = v - j as f64;

            // TSC 1D weights
            let w_u = [
                0.5 * (0.5 - du).powi(2),
                0.75 - du.powi(2),
                0.5 * (0.5 + du).powi(2),
            ];

            let w_v = [
                0.5 * (0.5 - dv).powi(2),
                0.75 - dv.powi(2),
                0.5 * (0.5 + dv).powi(2),
            ];

            for di in -1..=1 {
                let target_i = i + di;
                let final_i = target_i.clamp(0, self.nr as i32) as usize;

                for dj in -1..=1 {
                    let target_j = j + dj;
                    // Mirror reflection for Neumann boundaries in z
                    let final_j = if target_j < 0 {
                        (-target_j) as usize
                    } else if target_j > self.nz as i32 {
                        (2 * self.nz as i32 - target_j) as usize
                    } else {
                        target_j as usize
                    };

                    let weight = w_u[(di + 1) as usize] * w_v[(dj + 1) as usize];
                    q_grid[[final_i, final_j]] += q_p[p] * weight;
                }
            }
        }

        // Convert charge to density: rho_ij = q_grid_ij / V_ij
        let mut rho = Array2::<f64>::zeros((self.nr + 1, self.nz + 1));
        for i in 0..=self.nr {
            for j in 0..=self.nz {
                rho[[i, j]] = q_grid[[i, j]] / self.cell_volume(i, j);
            }
        }

        rho
    }

    /// Interpolates the electric field (Er, Ez) to logical coordinates (u, v) using bilinear interpolation.
    pub fn interpolate_field_cic(
        &self,
        e_r: &Array2<f64>,
        e_z: &Array2<f64>,
        u: f64,
        v: f64,
    ) -> (f64, f64) {
        let i = (u.floor() as usize).min(self.nr - 1);
        let j = (v.floor() as usize).min(self.nz - 1);

        let du = u - i as f64;
        let dv = v - j as f64;

        let interp = |f: &Array2<f64>| {
            (1.0 - du) * (1.0 - dv) * f[[i, j]]
                + du * (1.0 - dv) * f[[i + 1, j]]
                + (1.0 - du) * dv * f[[i, j + 1]]
                + du * dv * f[[i + 1, j + 1]]
        };

        (interp(e_r), interp(e_z))
    }

    /// Interpolates the electric field (Er, Ez) to logical coordinates (u, v) using quadratic spline (TSC).
    pub fn interpolate_field_tsc(
        &self,
        e_r: &Array2<f64>,
        e_z: &Array2<f64>,
        u: f64,
        v: f64,
    ) -> (f64, f64) {
        let i = u.round() as i32;
        let j = v.round() as i32;

        let du = u - i as f64;
        let dv = v - j as f64;

        let w_u = [
            0.5 * (0.5 - du).powi(2),
            0.75 - du.powi(2),
            0.5 * (0.5 + du).powi(2),
        ];

        let w_v = [
            0.5 * (0.5 - dv).powi(2),
            0.75 - dv.powi(2),
            0.5 * (0.5 + dv).powi(2),
        ];

        let mut er_interp = 0.0;
        let mut ez_interp = 0.0;

        for di in -1..=1 {
            let target_i = i + di;
            let final_i = target_i.clamp(0, self.nr as i32) as usize;

            for dj in -1..=1 {
                let target_j = j + dj;
                let final_j = if target_j < 0 {
                    (-target_j) as usize
                } else if target_j > self.nz as i32 {
                    (2 * self.nz as i32 - target_j) as usize
                } else {
                    target_j as usize
                };

                let w = w_u[(di + 1) as usize] * w_v[(dj + 1) as usize];
                er_interp += e_r[[final_i, final_j]] * w;
                ez_interp += e_z[[final_i, final_j]] * w;
            }
        }

        (er_interp, ez_interp)
    }

    /// Checks the Debye length resolution criterion: grid spacing delta_r <= lambda_D / 2
    /// at all nodes where the space charge density exceeds 1% of the maximum density in the domain.
    pub fn check_debye_length(
        &self,
        rho: &Array2<f64>,
        v0: f64,       // physical initial thermal velocity
        e_0: f64,      // physical electric field scale
        l_0: f64,      // physical length scale (R_a)
        e_over_m: f64, // physical specific charge
    ) -> Result<(), String> {
        let mut max_rho = 0.0;
        for &val in rho.iter() {
            if val.abs() > max_rho {
                max_rho = val.abs();
            }
        }

        let threshold = 0.01 * max_rho;
        if threshold < 1e-12 {
            return Ok(()); // Space charge is negligible everywhere
        }

        for i in 0..self.nr {
            let dr = self.r[i + 1] - self.r[i];
            for j in 0..=self.nz {
                let rho_val = rho[[i, j]];
                if rho_val.abs() > threshold {
                    // Compute normalized Debye length
                    // lambda_D = sqrt( v0^2 / ( |rho| * l_0 * e_0 * (e/m) ) )
                    let den = rho_val.abs() * l_0 * e_0 * e_over_m;
                    if den > 1e-20 {
                        let lambda_d = (v0 * v0 / den).sqrt();
                        if dr > lambda_d / 2.0 {
                            return Err(format!(
                                "Debye length resolution check failed: at node i={}, j={}, radial grid spacing dr={:.6e} exceeds lambda_D/2={:.6e} (local density rho={:.6e}, threshold={:.6e})",
                                i, j, dr, lambda_d / 2.0, rho_val, threshold
                            ));
                        }
                    }
                }
            }
        }

        Ok(())
    }
}

/// Preconditioned Conjugate Gradient Poisson solver for cylindrical coordinates.
pub struct PoissonSolver {
    pub grid: PoissonGrid,
    pub potential: Array2<f64>,
    pub e_r: Array2<f64>,
    pub e_z: Array2<f64>,
    pub diag_precond: Array2<f64>, // Jacobi preconditioner diagonal values (1 / M_ii)
}

impl PoissonSolver {
    pub fn new(grid: PoissonGrid) -> Self {
        let nr = grid.nr;
        let nz = grid.nz;
        let potential = Array2::<f64>::zeros((nr + 1, nz + 1));
        let e_r = Array2::<f64>::zeros((nr + 1, nz + 1));
        let e_z = Array2::<f64>::zeros((nr + 1, nz + 1));

        // Precompute Jacobi preconditioner diagonal elements.
        // We'll build the SPD operator A_SPD = -M, and get its diagonal.
        let mut diag_precond = Array2::<f64>::zeros((nr + 1, nz + 1));
        for i in 1..nr {
            let h_r_curr = grid.r[i + 1] - grid.r[i];
            let h_r_prev = grid.r[i] - grid.r[i - 1];
            let r_half_plus = (grid.r[i] + grid.r[i + 1]) / 2.0;
            let r_half_minus = (grid.r[i - 1] + grid.r[i]) / 2.0;

            for j in 0..=nz {
                let dz = if j == 0 || j == nz {
                    grid.delta_z / 2.0
                } else {
                    grid.delta_z
                };

                // M_r contribution: -2 * pi * dz_j * [r_{i+1/2}/h_r + r_{i-1/2}/h_{r-1}]
                let diag_r = 2.0 * PI * dz * (r_half_plus / h_r_curr + r_half_minus / h_r_prev);

                // M_z contribution: -2 * pi * r_i * dr_i * [2 / h_z] (for interior) or [1 / h_z] (for boundary)
                let dr = (grid.r[i + 1] - grid.r[i - 1]) / 2.0;
                let diag_z = if j == 0 || j == nz {
                    2.0 * PI * grid.r[i] * dr / grid.delta_z
                } else {
                    4.0 * PI * grid.r[i] * dr / grid.delta_z
                };

                diag_precond[[i, j]] = 1.0 / (diag_r + diag_z);
            }
        }

        Self {
            grid,
            potential,
            e_r,
            e_z,
            diag_precond,
        }
    }

    /// Computes the matrix-vector product A_SPD * x = -M * x.
    /// Assumes homogeneous boundary conditions (x = 0 at i=0 and i=nr).
    pub fn multiply_m_spd(&self, x: &Array2<f64>, out: &mut Array2<f64>) {
        let nr = self.grid.nr;
        let nz = self.grid.nz;
        out.fill(0.0);

        for i in 1..nr {
            let h_r_curr = self.grid.r[i + 1] - self.grid.r[i];
            let h_r_prev = self.grid.r[i] - self.grid.r[i - 1];
            let r_half_plus = (self.grid.r[i] + self.grid.r[i + 1]) / 2.0;
            let r_half_minus = (self.grid.r[i - 1] + self.grid.r[i]) / 2.0;
            let dr = (self.grid.r[i + 1] - self.grid.r[i - 1]) / 2.0;

            for j in 0..=nz {
                let dz = if j == 0 || j == nz {
                    self.grid.delta_z / 2.0
                } else {
                    self.grid.delta_z
                };

                // -T_r term:
                let val_r = 2.0
                    * PI
                    * dz
                    * (r_half_plus * (x[[i, j]] - if i + 1 == nr { 0.0 } else { x[[i + 1, j]] })
                        / h_r_curr
                        + r_half_minus
                            * (x[[i, j]] - if i - 1 == 0 { 0.0 } else { x[[i - 1, j]] })
                            / h_r_prev);

                // -T_z term:
                let val_z = if j == 0 {
                    2.0 * PI * self.grid.r[i] * dr * (x[[i, 0]] - x[[i, 1]]) / self.grid.delta_z
                } else if j == nz {
                    2.0 * PI * self.grid.r[i] * dr * (x[[i, nz]] - x[[i, nz - 1]])
                        / self.grid.delta_z
                } else {
                    2.0 * PI
                        * self.grid.r[i]
                        * dr
                        * (2.0 * x[[i, j]] - x[[i, j + 1]] - x[[i, j - 1]])
                        / self.grid.delta_z
                };

                out[[i, j]] = val_r + val_z;
            }
        }
    }

    /// Solves the Poisson equation M * phi = -V * rho using Jacobi-preconditioned conjugate gradient.
    /// Dirichlet boundary conditions: phi = V(r_c, z) = z * (U_f / l_c) at cathode (i=0) and phi = U_a at anode (i=nr).
    /// Neumann boundary conditions: dphi/dz = 0 at z_min (j=0) and z_max (j=nz).
    pub fn solve(
        &mut self,
        rho: &Array2<f64>,
        u_a: f64,
        u_f: f64,
        tol: f64,
        max_iter: usize,
    ) -> usize {
        let nr = self.grid.nr;
        let nz = self.grid.nz;
        let l_c_norm = self.grid.z_max - self.grid.z_min;

        // 1. Build the right hand side b_SPD = V * rho - M * phi_boundary
        let mut b = Array2::<f64>::zeros((nr + 1, nz + 1));
        for i in 1..nr {
            let dr = (self.grid.r[i + 1] - self.grid.r[i - 1]) / 2.0;
            for j in 0..=nz {
                let dz = if j == 0 || j == nz {
                    self.grid.delta_z / 2.0
                } else {
                    self.grid.delta_z
                };
                let vol = 2.0 * PI * self.grid.r[i] * dr * dz;
                b[[i, j]] = vol * rho[[i, j]];
            }
        }

        // Boundary contribution to the right-hand side (since we solved M_SPD * phi = b_SPD)
        for j in 0..=nz {
            let dz = if j == 0 || j == nz {
                self.grid.delta_z / 2.0
            } else {
                self.grid.delta_z
            };

            // Anode contribution (i = nr) to cell i = nr - 1
            let h_r_anode = self.grid.r[nr] - self.grid.r[nr - 1];
            let r_half_anode = (self.grid.r[nr - 1] + self.grid.r[nr]) / 2.0;
            b[[nr - 1, j]] += 2.0 * PI * dz * r_half_anode * u_a / h_r_anode;

            // Cathode contribution (i = 0) to cell i = 1
            let h_r_cathode = self.grid.r[1] - self.grid.r[0];
            let r_half_cathode = (self.grid.r[0] + self.grid.r[1]) / 2.0;
            let phi_cathode_j = self.grid.z[j] * (u_f / l_c_norm);
            b[[1, j]] += 2.0 * PI * dz * r_half_cathode * phi_cathode_j / h_r_cathode;
        }

        // 2. Initialize PCG variables.
        // We use the current potential as the initial guess to leverage temporal coherence.
        // Clamp boundaries of potential first
        for j in 0..=nz {
            self.potential[[0, j]] = self.grid.z[j] * (u_f / l_c_norm);
            self.potential[[nr, j]] = u_a;
        }

        // Compute initial residual r0 = b - A * phi0
        let mut temp = Array2::<f64>::zeros((nr + 1, nz + 1));
        self.multiply_m_spd(&self.potential, &mut temp);

        let mut r = Array2::<f64>::zeros((nr + 1, nz + 1));
        for i in 1..nr {
            for j in 0..=nz {
                r[[i, j]] = b[[i, j]] - temp[[i, j]];
            }
        }

        // Check initial convergence
        let mut r_norm = 0.0;
        let mut b_norm = 0.0;
        for i in 1..nr {
            for j in 0..=nz {
                r_norm += r[[i, j]].powi(2);
                b_norm += b[[i, j]].powi(2);
            }
        }
        r_norm = r_norm.sqrt();
        b_norm = b_norm.sqrt();

        if b_norm < 1e-12 {
            b_norm = 1.0;
        }

        if r_norm / b_norm < tol {
            self.compute_electric_fields();
            return 0;
        }

        // z0 = P^-1 * r0
        let mut z = Array2::<f64>::zeros((nr + 1, nz + 1));
        for i in 1..nr {
            for j in 0..=nz {
                z[[i, j]] = r[[i, j]] * self.diag_precond[[i, j]];
            }
        }

        // p0 = z0
        let mut p = z.clone();

        // gamma0 = dot(r0, z0)
        let mut gamma = 0.0;
        for i in 1..nr {
            for j in 0..=nz {
                gamma += r[[i, j]] * z[[i, j]];
            }
        }

        let mut iter = 0;
        let mut w = Array2::<f64>::zeros((nr + 1, nz + 1));

        while iter < max_iter {
            // w_k = A * p_k
            self.multiply_m_spd(&p, &mut w);

            // alpha_k = gamma_k / dot(p_k, w_k)
            let mut p_dot_w = 0.0;
            for i in 1..nr {
                for j in 0..=nz {
                    p_dot_w += p[[i, j]] * w[[i, j]];
                }
            }
            if p_dot_w.abs() < 1e-20 {
                break;
            }
            let alpha = gamma / p_dot_w;

            // phi_k+1 = phi_k + alpha * p_k
            // r_k+1 = r_k - alpha * w_k
            let mut r_norm_sq = 0.0;
            for i in 1..nr {
                for j in 0..=nz {
                    self.potential[[i, j]] += alpha * p[[i, j]];
                    r[[i, j]] -= alpha * w[[i, j]];
                    r_norm_sq += r[[i, j]].powi(2);
                }
            }
            let r_norm_curr = r_norm_sq.sqrt();

            if r_norm_curr / b_norm < tol {
                iter += 1;
                break;
            }

            // z_k+1 = P^-1 * r_k+1
            for i in 1..nr {
                for j in 0..=nz {
                    z[[i, j]] = r[[i, j]] * self.diag_precond[[i, j]];
                }
            }

            // gamma_k+1 = dot(r_k+1, z_k+1)
            let mut gamma_new = 0.0;
            for i in 1..nr {
                for j in 0..=nz {
                    gamma_new += r[[i, j]] * z[[i, j]];
                }
            }

            // beta_k = gamma_k+1 / gamma_k
            let beta = gamma_new / gamma;
            gamma = gamma_new;

            // p_k+1 = z_k+1 + beta * p_k
            for i in 1..nr {
                for j in 0..=nz {
                    p[[i, j]] = z[[i, j]] + beta * p[[i, j]];
                }
            }

            iter += 1;
        }

        self.compute_electric_fields();
        iter
    }

    /// Computes the electric fields Er and Ez from the solved potential.
    /// E_r = -dphi/dr, E_z = -dphi/dz.
    pub fn compute_electric_fields(&mut self) {
        let nr = self.grid.nr;
        let nz = self.grid.nz;

        // 1. Compute Er (radial electric field)
        // For interior radial nodes: central difference
        for i in 1..nr {
            let dr_curr = self.grid.r[i + 1] - self.grid.r[i];
            let dr_prev = self.grid.r[i] - self.grid.r[i - 1];

            for j in 0..=nz {
                // Non-uniform grid central difference formula:
                // df/dx = [ (f_{i+1} - f_i)*dx_{i-1}/dx_i + (f_i - f_{i-1})*dx_i/dx_{i-1} ] / (dx_i + dx_{i-1})
                self.e_r[[i, j]] = -((self.potential[[i + 1, j]] - self.potential[[i, j]])
                    * dr_prev
                    / dr_curr
                    + (self.potential[[i, j]] - self.potential[[i - 1, j]]) * dr_curr / dr_prev)
                    / (dr_curr + dr_prev);
            }
        }

        // Boundaries: one-sided differences
        for j in 0..=nz {
            let dr0 = self.grid.r[1] - self.grid.r[0];
            self.e_r[[0, j]] = -(self.potential[[1, j]] - self.potential[[0, j]]) / dr0;

            let dr_last = self.grid.r[nr] - self.grid.r[nr - 1];
            self.e_r[[nr, j]] = -(self.potential[[nr, j]] - self.potential[[nr - 1, j]]) / dr_last;
        }

        // 2. Compute Ez (axial electric field)
        // For interior axial nodes: central difference
        for i in 0..=nr {
            for j in 1..nz {
                self.e_z[[i, j]] = -(self.potential[[i, j + 1]] - self.potential[[i, j - 1]])
                    / (2.0 * self.grid.delta_z);
            }
            // Boundaries (dphi/dz = 0 at z_min and z_max)
            self.e_z[[i, 0]] = 0.0;
            self.e_z[[i, nz]] = 0.0;
        }
    }
}
