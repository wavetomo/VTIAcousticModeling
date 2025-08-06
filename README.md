# 2D VTI Acoustic Wave Simulation using Staggered Grid FD Method
Tests for symmetric and non-symmetric VTI acoustic wave equation

This project implements 2D acoustic wave forward modeling in vertical transverse isotropic (VTI) media using staggered grid finite difference (FD) schemes. It includes symmetric and non-symmetric formulations.

## 📁 File Structure

### 🔧 Main Forward Modeling Scripts

| File Name                                                          | Description                                                  |
|--------------------------------------------------------------------|--------------------------------------------------------------|
| `staggered_grid_FD_symmetric_VTIacoustic.m`                        | Symmetric VTI acoustic wave forward modeling main script     |
| `staggered_grid_FD_stable_adjoint_non_symmetric_VTIacoustic.m`     | Non-symmetric, stable formulation of VTI adjoint modeling    |
| `staggered_grid_FD_unstable_adjoint_non_symmetric_VTIacoustic.m`   | Non-symmetric, *unstable* formulation (for testing purpose)  |
| `staggered_grid_FD_nonsymmetric_VTIacoustic.m`                     | Non-symmetric, VTI acoustic wave forward modeling main script|
### 📐 Model Parameter Files (`.vel` and `.bin`)

These files store 2D spatial models in **float32 binary** format, saved **column-wise (Fortran-style)**.

#### Marmousi Model Files

| File Name                                  | Description                        | Size (nz × nx) |
|-------------------------------------------|------------------------------------|----------------|
| `MODEL_P-WAVE_VELOCITY_351x1301_10m.vel`  | P-wave velocity model              | 351 × 1301     |
| `MODEL_DENSITY_351x1301_10m.vel`          | Density model                      | 351 × 1301     |
| `MODEL_DELTA_351x1301_10m.vel`            | Thomsen anisotropy parameter δ     | 351 × 1301     |
| `MODEL_EPSLION_351x1301_10m.vel`          | Thomsen anisotropy parameter ε     | 351 × 1301     |

#### Hess Model Files

| File Name                       | Description              | Size (nz × nx) |
|--------------------------------|--------------------------|----------------|
| `vp_hess_366_882_25m.bin`      | Hess model for Vp        | 366 × 882      |
| `rho_hess_366_882_25m.bin`     | Hess model for density   | 366 × 882      |
| `delta_hess_366_882_25m.bin`   | Hess model for delta     | 366 × 882      |
| `epsilon_hess_366_882_25m.bin` | Hess model for epsilon   | 366 × 882      |

### 📂 Utility Scripts

| File Name         | Purpose                                                      |
|-------------------|--------------------------------------------------------------|
| `read_matrix.m`   | Read `.bin` and `.vel` binary matrix files into MATLAB       |
| `Ricker.m`        | Generate Ricker wavelet                                      |
| `boundary2dcut.m` | Apply absorbing boundary or sponge layers                    |
| `modpad2d.m`      | Pad 2D models for boundary conditions                        |
| `meal2d.m`        | Possibly for wavefield visualization or post-processing      |

### 📄 Documentation

| File Name                                       | Description                           |
|------------------------------------------------|---------------------------------------|
| `README.md`                                    | Project description and file guide    |
| `Tests for symmetric and non-symmetric VTI acoustic wave.pptx` | Test case illustrations |

---

## 🚀 Quick Start

```matlab
% Execute the script of symmetric VTI acoustic forward modeling in matlab
staggered_grid_FD_nonsymmetric_VTIacoustic.m
staggered_grid_FD_stable_adjoint_non_symmetric_VTIacoustic.m
staggered_grid_FD_symmetric_VTIacoustic.m
staggered_grid_FD_unstable_adjoint_non_symmetric_VTIacoustic.m



