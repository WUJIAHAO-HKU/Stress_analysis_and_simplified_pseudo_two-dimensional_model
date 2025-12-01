# 🔋 JuBat - Battery Modeling Framework

<div align="center">

**A high-performance Julia-based framework for advanced battery modeling and stress analysis**

[![Julia](https://img.shields.io/badge/Julia-1.6+-9558B2?style=flat&logo=julia)](https://julialang.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.softx.2024.101760-blue)](https://doi.org/10.1016/j.softx.2024.101760)

</div>

---

## 📋 Table of Contents

- [Overview](#-overview)
- [Key Features](#-key-features)
- [Models Supported](#-models-supported)
- [Installation](#-installation)
- [Quick Start](#-quick-start)
- [Examples](#-examples)
- [Project Structure](#-project-structure)
- [Research Applications](#-research-applications)
- [Contributing](#-contributing)
- [References](#-references)
- [License](#-license)

---

## 🎯 Overview

**JuBat** is a cutting-edge Julia-based framework for battery modeling, built on Newman's electrochemical models. It combines the computational efficiency of Julia with high-accuracy second-order Finite Element Method (FEM) to provide fast, accurate, and reliable battery simulations.

### Why JuBat?

- 🚀 **High Performance**: Leverages Julia's speed for rapid computations
- 🎯 **High Accuracy**: 2nd order FEM for superior precision
- 📐 **Comprehensive Models**: P2D, SPM, SPMe, and simplified sP2D models
- 🔧 **Extensible**: Easy to add custom models and features
- 📊 **Mechanical Stress Analysis**: Advanced stress analysis capabilities for battery safety research

---

## ✨ Key Features

### Core Capabilities

- **Multi-Physics Modeling**: Electrochemical, thermal, and mechanical coupling
- **Advanced Models**: Full P2D, SPM, SPMe, and innovative sP2D models
- **Mechanical Stress Analysis**: Surface tangential stress prediction for degradation studies
- **Flexible Current Profiles**: Support for dynamic loading conditions (UDDS, US06, pulse, etc.)
- **High-Order FEM**: 2nd order finite element implementation for better convergence
- **Parameter Optimization**: Built-in tools for model parameter identification

### Battery Cells Supported

- ✅ Enertech cell parameters
- ✅ LG M50 cell parameters
- ✅ Northrop reference parameters
- 🔧 Easy to add custom cell parameters

---

## 🔬 Models Supported

### 1. **P2D (Pseudo-Two-Dimensional) Model**
The full electrochemical model based on Newman's porous electrode theory, solving coupled PDEs for solid and electrolyte phases.

### 2. **SPM (Single Particle Model)**
Simplified model assuming uniform current distribution, suitable for fast simulations.

### 3. **SPMe (Single Particle Model with Electrolyte)**
Enhanced SPM including electrolyte dynamics for improved accuracy.

### 4. **sP2D (Simplified Pseudo-2D) Model**
Innovative simplified P2D model using piecewise polynomial approximations, achieving near-P2D accuracy with significantly reduced computational cost.

### 5. **Thermal Model**
Coupled thermal analysis for temperature distribution prediction.

### 6. **Mechanical Stress Model**
Advanced mechanical stress analysis for particle surface stress prediction and degradation studies.

---

## 🛠 Installation

### Prerequisites

- Julia 1.6 or higher
- Required Julia packages: `LinearAlgebra`, `SparseArrays`, `Plots`, `Parameters`, `CSV`, `DataFrames`

### Setup

1. Clone the repository:
```bash
git clone https://github.com/WUJIAHAO-HKU/Stress_analysis_and_simplified_pseudo_two-dimensional_model.git
cd Stress_analysis_and_simplified_pseudo_two-dimensional_model
```

2. Install dependencies:
```julia
using Pkg
include("src/install.jl")
```

---

## 🚀 Quick Start

### Basic P2D Simulation

```julia
using Plots
include("src/JuBat.jl")

# Select battery cell parameters
param_dim = JuBat.ChooseCell("Enertech")

# Configure simulation options
opt = JuBat.Option()
opt.mechanicalmodel = "full"  # Enable stress analysis
Crate = 1  # C-rate
opt.Current = x -> 2.28 * Crate  # Constant current
opt.time = [0, 3600]  # Simulation time [s]
opt.model = "P2D"  # Choose model: "P2D", "SPM", "SPMe", or "sP2D"

# Run simulation
case1 = JuBat.SetCase(param_dim, opt)
result = JuBat.Solve(case1)

# Plot results
time = result["time [s]"]
voltage = result["cell voltage [V]"]
plot(time, voltage, xlabel="Time [s]", ylabel="Voltage [V]")
```

### Mechanical Stress Analysis

```julia
# Extract stress data
stress = result["negative particle surface tangential stress[Pa]"][1, :]
concentration = result["negative particle surface lithium concentration [mol/m^3]"][1, :]

# Visualize stress evolution
plot(time, stress, xlabel="Time [s]", ylabel="Stress [Pa]", 
     title="Surface Tangential Stress During Discharge")
```

---

## 📚 Examples

The `example/` directory contains comprehensive examples:

- **`minimal_example.jl`**: Basic P2D simulation
- **`thermal_example.jl`**: Thermal coupling simulation
- **`mechanical_example.jl`**: Full mechanical stress analysis
- **`sP2D_example.jl`**: Simplified P2D model demonstration
- **`change_current.jl`**: Dynamic current profile examples
- **`change_model.jl`**: Switching between different models

For advanced applications, see the `最新的进展/` directory for:
- Dynamic loading conditions (UDDS, US06 drive cycles)
- High-frequency pulse current analysis
- Fast charge/discharge scenarios
- Multi-rate discharge comparisons
- Spatial distribution analysis

---

## 📂 Project Structure

```
.
├── src/                          # Core source code
│   ├── JuBat.jl                 # Main module
│   ├── P2D.jl                   # P2D model implementation
│   ├── sP2D.jl                  # Simplified P2D model
│   ├── SPM.jl / SPMe.jl         # Single particle models
│   ├── Mechanical.jl            # Stress analysis
│   ├── Thermal.jl               # Thermal model
│   ├── SetParams.jl             # Parameter management
│   ├── Solve.jl                 # Solver routines
│   └── parameters/              # Battery cell parameters
│       ├── Enertech.jl
│       ├── LGM50.jl
│       └── Northrop.jl
├── example/                      # Example scripts
├── 最新的进展/                   # Advanced research applications
├── .vscode/                      # VS Code configuration
├── README.md                     # This file
└── LICENSE                       # License information
```

---

## 🔬 Research Applications

This framework has been used for cutting-edge battery research:

### Stress Analysis & Safety
- Surface tangential stress prediction during charge/discharge cycles
- Mechanical degradation analysis
- Particle fracture risk assessment

### Model Comparison & Validation
- P2D vs. sP2D accuracy-efficiency trade-offs
- Multi-rate discharge performance prediction
- Validation against experimental data

### Dynamic Operating Conditions
- Electric vehicle drive cycles (UDDS, US06, WLTC)
- High-frequency pulse current response
- Fast charging protocols
- Temperature-dependent behavior

### Parameter Optimization
- Model calibration using experimental data
- Uncertainty quantification
- Sensitivity analysis

---

## 🤝 Contributing

We welcome contributions from the community! Here's how you can help:

1. **Add New Features**: Implement new models or analysis tools
2. **Add Battery Parameters**: Contribute validated parameter sets
3. **Improve Documentation**: Enhance examples and tutorials
4. **Report Issues**: Help us identify and fix bugs
5. **Share Research**: Add your publications to `Citation.jl`

### Contribution Guidelines

- Ensure code is well-documented and includes examples
- Add references to `Citation.jl` for proper attribution
- Follow Julia best practices and coding style
- Test your contributions thoroughly

---

## 📖 References

### Primary Publications

[1] **W. Ai, Y. Liu**, "JuBat: A Julia-based framework for battery modelling using finite element method", *SoftwareX*, Vol. 27, 2024, 101760.  
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.softx.2024.101760-blue)](https://doi.org/10.1016/j.softx.2024.101760)

[2] **W. Ai, Y. Liu**, "Improving the convergence rate of Newman's battery model using 2nd order finite element method", *Journal of Energy Storage*, Vol. 67, 2023, 107512.  
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.est.2023.107512-blue)](https://doi.org/10.1016/j.est.2023.107512)

[3] **W. Ai, Y. Liu**, "sP2D: Simplified pseudo 2D battery model by piecewise sinusoidal/quadratic functions of potential curves", *Journal of Energy Storage*, Vol. 86, 2024, 111386.  
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.est.2024.111386-blue)](https://doi.org/10.1016/j.est.2024.111386)

### Related Work

For more research applications and validation studies, please refer to the publications listed in `src/Citation.jl`.

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- Thanks to all contributors and researchers using JuBat
- Built on the foundation of Newman's pioneering battery models
- Supported by the Julia community

---

## 📧 Contact & Support

For questions, suggestions, or collaborations:

- 🐛 **Issues**: [GitHub Issues](https://github.com/WUJIAHAO-HKU/Stress_analysis_and_simplified_pseudo_two-dimensional_model/issues)
- 💬 **Discussions**: [GitHub Discussions](https://github.com/WUJIAHAO-HKU/Stress_analysis_and_simplified_pseudo_two-dimensional_model/discussions)
- 📧 **Email**: Contact the repository owner

---

<div align="center">

**⭐ Star this repository if you find it useful! ⭐**

Made with ❤️ for the battery modeling community

</div>

