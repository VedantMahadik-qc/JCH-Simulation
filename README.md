# Jaynes-Cummings-Hubbard (JCH) Model Simulation

This repository contains simulations of the Jaynes-Cummings-Hubbard model, investigating photon dynamics, phase transitions, and lattice effects. The project implements simulations using three different frameworks to benchmark performance and verify results:

1.  **Python**: Using the [QuTiP](http://qutip.org/) library.
2.  **Julia**: Using the [QuantumOptics.jl](https://qojulia.org/) framework.
3.  **Julia**: Using [QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl) (QuTiP-like syntax in Julia).

## Project Structure

```text
jch-simulation/
├── Qutip/                  # Python implementations using QuTiP
│   ├── Cavity-Cavity.../   # Cavity interaction simulations
│   └── Jaynes Cummings.../ # Single-site JC model benchmarks
├── QuantumOptics.jl/       # Julia implementations using QuantumOptics.jl
│   ├── Cavity-Cavity/      # Sparse matrix operations & scaling tests
│   └── Jaynes-Cummings/    # Fundamental JC dynamics
├── QuantumToolbox.jl/      # Julia implementations using QuantumToolbox.jl
└── 2D_Lattice/             # 2D Lattice simulations (Photon spreading, etc.)