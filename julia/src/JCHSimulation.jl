"""
    JCHSimulation

Jaynes-Cummings-Hubbard model on 1-D chains and 2-D lattices, built on
QuantumOptics.jl.

Conventions (ħ = 1), matching the Python `jch` package:

    H = Σ_i [ ω_c a_i† a_i + (ω_a/2) σ_z^i + g (a_i† σ_-^i + a_i σ_+^i) ]
        - J Σ_<ij> (a_i† a_j + h.c.)

Note the minus sign on the hopping term (Bose-Hubbard convention).
"""
module JCHSimulation

using QuantumOptics
using LinearAlgebra

include("hamiltonians.jl")
include("meanfield.jl")

export JCHParams,
       build_jch_chain,
       build_2d_cavity_lattice,
       site_basis,
       neighbour_table,
       meanfield_lattice,
       photon_numbers

end # module
