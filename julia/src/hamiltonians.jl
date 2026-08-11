"""
    JCHParams

Parameters of a Jaynes-Cummings-Hubbard lattice.

- `ω_c`      cavity frequency
- `ω_a`      atomic transition frequency
- `g`        atom-cavity coupling
- `J`        photon hopping amplitude between neighbouring cavities
- `n_cutoff` Fock truncation: each cavity keeps |0⟩ … |n_cutoff-1⟩
"""
Base.@kwdef struct JCHParams
    ω_c::Float64 = 1.0
    ω_a::Float64 = 1.0
    g::Float64 = 0.05
    J::Float64 = 0.1
    n_cutoff::Int = 3

    function JCHParams(ω_c, ω_a, g, J, n_cutoff)
        n_cutoff >= 2 || throw(ArgumentError("n_cutoff must be ≥ 2, got $n_cutoff"))
        g >= 0 || throw(ArgumentError("g must be non-negative, got $g"))
        new(ω_c, ω_a, g, J, n_cutoff)
    end
end

"""
    site_basis(p::JCHParams)

Local Hilbert space of one site: Fock(n_cutoff) ⊗ Spin(1/2),
together with the local operators `(a, σ_-, σ_z, n)`.
"""
function site_basis(p::JCHParams)
    b_cav = FockBasis(p.n_cutoff - 1)
    b_atom = SpinBasis(1//2)
    b_site = b_cav ⊗ b_atom

    a = destroy(b_cav) ⊗ one(b_atom)
    sm = one(b_cav) ⊗ sigmam(b_atom)
    sz = one(b_cav) ⊗ sigmaz(b_atom)
    n = number(b_cav) ⊗ one(b_atom)

    return b_site, (a = a, sm = sm, sz = sz, n = n)
end

"""
    build_jch_chain(p::JCHParams, n_sites::Int)

Full (exact) JCH Hamiltonian on a 1-D open chain.

Returns `(H, b_total, a_ops, n_ops)`.

!!! warning "Exponential wall"
    `dim = (n_cutoff * 2)^n_sites`. At `n_cutoff = 3`, six sites is already
    46 656 states and eight sites is 1.7 million. Use `meanfield_lattice`
    beyond that.
"""
function build_jch_chain(p::JCHParams, n_sites::Int)
    n_sites >= 1 || throw(ArgumentError("n_sites must be ≥ 1"))
    b_site, loc = site_basis(p)
    b_total = CompositeBasis([b_site for _ in 1:n_sites]...)

    a  = [embed(b_total, i, loc.a)  for i in 1:n_sites]
    sm = [embed(b_total, i, loc.sm) for i in 1:n_sites]
    sz = [embed(b_total, i, loc.sz) for i in 1:n_sites]
    n  = [embed(b_total, i, loc.n)  for i in 1:n_sites]

    H = SparseOperator(b_total)
    for i in 1:n_sites
        H += p.ω_c * a[i]' * a[i] + 0.5 * p.ω_a * sz[i]
        H += p.g * (a[i]' * sm[i] + a[i] * sm[i]')
        if i < n_sites
            H -= p.J * (a[i]' * a[i+1] + a[i] * a[i+1]')
        end
    end

    return H, b_total, a, n
end

"""
    neighbour_table(Lx, Ly)

Nearest-neighbour indices for an `Lx × Ly` grid with open boundaries.
Site `(x, y)` maps to linear index `i = x + (y-1)*Lx`.
"""
function neighbour_table(Lx::Int, Ly::Int)
    N = Lx * Ly
    neighbours = [Int[] for _ in 1:N]
    for y in 1:Ly, x in 1:Lx
        i = x + (y - 1) * Lx
        if x < Lx
            j = (x + 1) + (y - 1) * Lx
            push!(neighbours[i], j); push!(neighbours[j], i)
        end
        if y < Ly
            j = x + y * Lx
            push!(neighbours[i], j); push!(neighbours[j], i)
        end
    end
    return neighbours
end

"""
    build_2d_cavity_lattice(Lx, Ly, ω_c, J, n_cutoff)

Bare-cavity (no atoms) 2-D lattice — the photonic tight-binding model used for
the wavepacket-spreading animations.

Returns `(H, b_total, n_ops)`.
"""
function build_2d_cavity_lattice(Lx::Int, Ly::Int, ω_c::Float64, J::Float64, n_cutoff::Int)
    N = Lx * Ly
    b_mode = FockBasis(n_cutoff - 1)
    b_total = tensor([b_mode for _ in 1:N]...)

    a = destroy(b_mode); at = create(b_mode); n = number(b_mode)
    n_ops = [embed(b_total, i, n) for i in 1:N]

    H = SparseOperator(b_total)
    neighbours = neighbour_table(Lx, Ly)
    for i in 1:N
        H += ω_c * n_ops[i]
        for j in neighbours[i]
            j > i || continue                      # count each bond once
            H -= J * (embed(b_total, i, at) * embed(b_total, j, a) +
                      embed(b_total, i, a)  * embed(b_total, j, at))
        end
    end

    return H, b_total, n_ops
end
