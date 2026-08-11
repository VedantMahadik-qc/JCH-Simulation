using JCHSimulation
using QuantumOptics
using LinearAlgebra
using Test

@testset "JCHSimulation" begin

    @testset "parameter validation" begin
        @test_throws ArgumentError JCHParams(n_cutoff = 1)
        @test_throws ArgumentError JCHParams(g = -1.0)
        @test JCHParams().n_cutoff == 3
    end

    @testset "chain Hamiltonian" begin
        p = JCHParams(n_cutoff = 3, ω_c = 1.0, ω_a = 1.0, g = 0.05, J = 0.1)
        H, b, a, n = build_jch_chain(p, 2)

        @test H isa AbstractOperator
        @test ishermitian(dense(H))
        @test length(b) == (3 * 2)^2          # (n_cutoff · 2)^n_sites
        @test length(a) == 2 && length(n) == 2
    end

    @testset "single site reduces to Jaynes-Cummings" begin
        g = 0.05
        p = JCHParams(n_cutoff = 12, ω_c = 1.0, ω_a = 1.0, g = g, J = 0.0)
        H, _, _, _ = build_jch_chain(p, 1)
        e = sort(real.(eigenenergies(dense(H))))
        # first excitation manifold splits by exactly 2g on resonance
        @test isapprox(e[3] - e[2], 2g; atol = 1e-9)
    end

    @testset "total excitation number is conserved" begin
        p = JCHParams(n_cutoff = 4, J = 0.1)
        H, b, a, n = build_jch_chain(p, 2)
        _, loc = site_basis(p)
        Nexc = sum(n) + sum(embed(b, i, 0.5 * (loc.sz + one(loc.sz))) for i in 1:2)
        @test norm((dense(H) * dense(Nexc) - dense(Nexc) * dense(H)).data) < 1e-10
    end

    @testset "neighbour table" begin
        nb = neighbour_table(3, 3)
        @test length(nb) == 9
        @test sort(nb[5]) == [2, 4, 6, 8]     # centre of a 3x3 grid
        @test sort(nb[1]) == [2, 4]           # corner
        @test sum(length, nb) == 2 * 12       # 12 bonds, counted twice
    end

    @testset "2D bare-cavity lattice" begin
        H, b, n_ops = build_2d_cavity_lattice(2, 2, 1.0, 0.1, 2)
        @test length(b) == 2^4
        @test length(n_ops) == 4
        @test ishermitian(dense(H))
    end

    @testset "mean field spreads the photon outward" begin
        p = JCHParams(n_cutoff = 3, J = 0.1)
        sol, b_site, ops, nb = meanfield_lattice(p, 3, 3; tspan = (0.0, 5.0), saveat = 0.5)
        n_t = photon_numbers(sol, b_site, ops, 9)

        @test size(n_t) == (length(sol.t), 9)
        @test all(n_t .>= -1e-10)              # occupations cannot be negative
        @test n_t[1, 5] > n_t[1, 1]            # starts localised at the centre
        @test n_t[end, 1] > n_t[1, 1]          # and leaks to the corner
    end
end
