using OrdinaryDiffEq

"""
    meanfield_lattice(p::JCHParams, Lx, Ly; α0, centre, tspan, saveat)

Gutzwiller / product-state mean-field dynamics on an `Lx × Ly` lattice.

Instead of one wavefunction in a Hilbert space of dimension `(2·n_cutoff)^N`,
we carry `N` *independent* site wavefunctions and couple them only through the
expectation value of their neighbours' field:

    |Ψ⟩ ≈ ⊗_i |ψ_i⟩,        α_i = ⟨ψ_i| a |ψ_i⟩

    H_eff^i = H_local  −  J ( a† Σ_{j∈nn(i)} α_j  +  a Σ_{j∈nn(i)} α_j* )

This turns an exponential problem into a linear one: cost scales as
`N · d_site` rather than `d_site^N`. The price is that inter-site entanglement
is discarded by construction, so the method is quantitative deep in the
superfluid phase and only qualitative near the Mott transition.

Returns `(sol, b_site, ops, neighbours)`.
"""
function meanfield_lattice(p::JCHParams, Lx::Int, Ly::Int;
                           α0::ComplexF64 = ComplexF64(1.0),
                           centre::Union{Nothing,Int} = nothing,
                           tspan::Tuple{Float64,Float64} = (0.0, 20.0),
                           saveat::Float64 = 0.1)
    N = Lx * Ly
    b_site, loc = site_basis(p)
    d_site = length(b_site)
    neighbours = neighbour_table(Lx, Ly)

    H_local = p.ω_c * loc.n + 0.5 * p.ω_a * loc.sz +
              p.g * (loc.a' * loc.sm + loc.a * loc.sm')

    function rhs!(dψ_flat, ψ_flat, _, _)
        kets = [Ket(b_site, ψ_flat[(i-1)*d_site+1 : i*d_site]) for i in 1:N]
        α = [expect(loc.a, k) for k in kets]

        for i in 1:N
            field = isempty(neighbours[i]) ? zero(ComplexF64) :
                    sum(α[j] for j in neighbours[i])
            H_eff = H_local - p.J * (loc.a' * field + loc.a * conj(field))
            dψ_flat[(i-1)*d_site+1 : i*d_site] = (-1.0im * (H_eff * kets[i])).data
        end
    end

    # A coherent seed is essential: a Fock state has ⟨a⟩ = 0, so under mean field
    # the neighbours would feel nothing and the photon would never spread.
    c = something(centre, cld(N, 2))
    kets0 = [fockstate(FockBasis(p.n_cutoff - 1), 0) ⊗ spindown(SpinBasis(1//2)) for _ in 1:N]
    kets0[c] = coherentstate(FockBasis(p.n_cutoff - 1), α0) ⊗ spindown(SpinBasis(1//2))
    u0 = reduce(vcat, [k.data for k in kets0])

    prob = ODEProblem(rhs!, u0, tspan)
    sol = solve(prob, Tsit5(); saveat = saveat)

    return sol, b_site, loc, neighbours
end

"""
    photon_numbers(sol, b_site, ops, N)

Reconstruct ⟨n_i⟩(t) from the flattened mean-field solution.
Returns a `length(sol.t) × N` matrix.
"""
function photon_numbers(sol, b_site, ops, N::Int)
    d_site = length(b_site)
    n_t = zeros(length(sol.t), N)
    for (k, u) in enumerate(sol.u), i in 1:N
        ψ = Ket(b_site, u[(i-1)*d_site+1 : i*d_site])
        n_t[k, i] = real(expect(ops.n, ψ))
    end
    return n_t
end
