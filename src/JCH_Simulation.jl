module JCHModels

using QuantumOptics

export build_jch_hamiltonian

"""
build_jch_hamiltonian(N_sites, wc, J, g, dim)
Constructs the Hamiltonian for a standard JCH chain.
"""
function build_jch_hamiltonian(N_sites::Int, wc::Float64, J::Float64, g::Float64, dim::Int)
    # Define basis
    b_site = FockBasis(dim)
    b_total = CompositeBasis([b_site for _ in 1:N_sites]...)

    # Operators
    a = [embed(b_total, i, destroy(b_site)) for i in 1:N_sites]
    sm = [embed(b_total, i, sigmam(b_site)) for i in 1:N_sites]

    # Build H (Hopping + Interaction)
    H = 0 * create(b_total) * destroy(b_total) # Zero operator

    for i in 1:N_sites
        # Cavity energy + Atom energy
        H += wc * a[i]' * a[i] + wc * sm[i]' * sm[i]
        # Interaction
        H += g * (a[i]' * sm[i] + a[i] * sm[i]')

        # Hopping (Periodic Boundary Conditions)
        if i < N_sites
            H -= J * (a[i]' * a[i+1] + a[i] * a[i+1]')
        end
    end
    return H, b_total
end

end