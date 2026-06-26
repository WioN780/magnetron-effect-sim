# M1 Pusher Convergence Order Report

This report documents the numerical convergence order of the Higuera-Cary relativistic particle pusher.

## Case A: Relativistic Larmor Orbit (Pure B)
- **Configuration**: Relativistic circular orbit, $v_0 = 0.5c$, $q/m = -1.0$, $B_z = 1.0$
- **Duration**: 2 full orbits

| Steps per Gyroperiod | Timestep $\Delta t$ | Position Error |
| :---: | :---: | :---: |
| 8 | 0.906900 | 2.42831081e-1 |
| 16 | 0.453450 | 5.90191769e-2 |
| 32 | 0.226725 | 1.46170451e-2 |
| 64 | 0.113362 | 3.64514820e-3 |
| 128 | 0.056681 | 9.10709637e-4 |

**Fitted Order of Accuracy**: **2.0135** (expected ≈ 2.0)

## Case B: E x B Drift
- **Configuration**: Perpendicular fields, $E_y = -0.001$, $B_z = 1.0$, $q/m = -1.0$
- **Duration**: 2 periods

| Steps per Gyroperiod | Timestep $\Delta t$ | Position Error |
| :---: | :---: | :---: |
| 8 | 0.785398 | 3.13427642e-4 |
| 16 | 0.392699 | 8.03581027e-5 |
| 32 | 0.196350 | 2.01725959e-5 |
| 64 | 0.098175 | 5.05241411e-6 |
| 128 | 0.049087 | 1.26863797e-6 |

**Fitted Order of Accuracy**: **1.9889** (expected ≈ 2.0)

## Analysis
Both the Relativistic Larmor Orbit test and the E x B Drift test confirm that the Higuera-Cary pusher is second-order accurate. The fitted orders of accuracy lie well within the expected theoretical range of 1.8 to 2.2, demonstrating correct implementation of the updates.