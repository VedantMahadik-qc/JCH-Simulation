using JCH_Simulation
using Test
using QuantumOptics

@testset "JCH Hamiltonian Tests" begin
    # 1. Define simple parameters
    N_sites = 2
    wc = 1.0
    J = 0.1
    g = 0.05
    dim = 3

    # 2. Build Hamiltonian using YOUR package
    H, b_total = build_jch_hamiltonian(N_sites, wc, J, g, dim)

    # 3. Verify it is a valid Operator
    @test isa(H, Operator)

    # 4. Verify it is Hermitian (Physical Hamiltonians must be Hermitian)
    @test ishermitian(H)

    # 5. Verify the dimension size (dim^N_sites)
    # For 2 sites with dimension 3, total size should be 3^2 = 9
    @test length(b_total) == 9
end