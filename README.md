# Seismic Properties of Pore Fluids (MATLAB)

## Table of Contents
- [Seismic Properties of Pore Fluids (MATLAB)](#seismic-properties-of-pore-fluids-matlab)
  - [Contents](#contents)
    - [Functions](#functions)
    - [Demo scripts](#demo-scripts)
  - [Requirements](#requirements)
  - [Quick start](#quick-start)
  - [Reference](#reference)
  - [Maintainer](#maintainer)

This repository contains **MATLAB** implementations of the Batzle & Wang (1992) correlations for **brine**, **dead oil**, and **hydrocarbon gas**.  
Each fluid has a function that computes **density (kg/m³)**, **acoustic velocity (m/s)**, and **bulk modulus (GPa)** and a companion **demo script** that sweeps one variable while keeping the others fixed and plots **velocity and density**.

All functions accept **arrays (1D/2D)** for the input variables (matching sizes) and include **explicit sanity checks** and `arguments` blocks for validation.

---

## Contents

### Functions
| File | Purpose | Inputs (units) | Outputs (units) | Notes |
|---|---|---|---|---|
| `brineBatzleWang.m` | Brine properties from Batzle–Wang | `T` (°C), `P` (MPa), `Salinity` (ppm) | `Rho_Brine` (kg/m³), `V_Brine` (m/s), `K_Brine` (GPa) | Equations **(27a)**, **(27b)**, **(28)**, **(29)** |
| `deadOilBatzleWang.m` | Dead oil (black oil) properties | `T` (°C), `P` (MPa), `API` (°API) | `Rho_Oil` (kg/m³), `V_Oil` (m/s), `K_Oil` (GPa) | Equations **(14)**, **(18)**, **(19)**, **(20b)** |
| `gasBatzleWang.m` | Hydrocarbon gas properties | `T` (°C), `P` (MPa), `G` (gas gravity) | `Rho_Gas` (kg/m³), `V_Gas` (m/s), `K_Gas` (GPa) | Equations **(9a–9b)**, **(10a–10c)**, **(11a–11b)** |

### Demo scripts
| Script | What it plots (dual y-axes) | Fixed inputs (defaults in script) |
|---|---|---|
| `brine_demo_script.m` | Velocity & density vs **salinity**, **temperature**, **pressure** | Holds the other two variables constant in each sweep |
| `dead_oil_demo_script.m` | Velocity & density vs **API**, **temperature**, **pressure** | Same pattern as above |
| `gas_demo_script.m` | Velocity & density vs **gas gravity G**, **temperature**, **pressure** | Same pattern as above |

> All scripts are **standalone**: just run them in MATLAB to regenerate the figures.

---

## Requirements

- MATLAB **R2019b or newer** (uses `arguments` blocks and modern language features).
- No toolboxes required.

---

## Quick start

### 1) Clone and open in MATLAB
```bash
git clone https://github.com/borgesf/BatzleWang.git
cd BatzleWang
```

### 2) Run any demo script
```matlab
brine_demo_script
dead_oil_demo_script
gas_demo_script
```

### 3) Call functions directly (arrays allowed)
```matlab
% Brine
T = [40 80]'; P = [20 20]'; S = [35000 60000]';
[rho_b, v_b, K_b] = brineBatzleWang(T, P, S);

% Dead oil
API = 35*ones(2,1);
[rho_o, v_o, K_o] = deadOilBatzleWang(T, P, API);

% Gas
G = 0.7*ones(2,1);
[rho_g, v_g, K_g] = gasBatzleWang(T, P, G);
```

---

## Reference

Batzle, M., & Wang, Z. (1992). Seismic properties of pore fluids. Geophysics, 57(11), 1396–1408. https://doi.org/10.1190/1.1443207

---

## Maintainer

FIAD (August 2025)

