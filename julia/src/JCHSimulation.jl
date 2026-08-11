module JCH_Simulation

using QuantumOptics

export build_jch_hamiltonian

"""
    build_jch_hamiltonian(N_sites, wc, J, g, dim_fock)

Constructs the Jaynes-Cummings-Hubbard Hamiltonian.
Each site is a Tensor Product: Cavity(Fock) ⊗ Atom(Spin).
"""
function build_jch_hamiltonian(N_sites::Int, wc::Float64, J::Float64, g::Float64, dim_fock::Int)
    # 1. Define Local Basis: Cavity ⊗ Atom
    b_cav = FockBasis(dim_fock)
    b_atom = SpinBasis(1//2)
    b_site = b_cav ⊗ b_atom
    
    # 2. Define Total Basis
    b_total = CompositeBasis([b_site for _ in 1:N_sites]...)
    
    # 3. Create Operators for Hamiltonian construction
    # a[i] = Destroy cavity photon at site i
    # sm[i] = Lower atom spin at site i
    # sz[i] = Sigma Z atom spin at site i
    a  = [embed(b_total, i, destroy(b_cav) ⊗ one(b_atom)) for i in 1:N_sites]
    sm = [embed(b_total, i, one(b_cav) ⊗ sigmam(b_atom)) for i in 1:N_sites]
    sz = [embed(b_total, i, one(b_cav) ⊗ sigmaz(b_atom)) for i in 1:N_sites]
    
    # 4. Build Hamiltonian
    H = DenseOperator(b_total) # Start with zeros
    
    for i in 1:N_sites
        # Site Energy: Cavity (ω*a'a) + Atom (0.5*ω*sz)
        H += wc * a[i]' * a[i] + 0.5 * wc * sz[i]
        
        # Interaction: g(a'sm + asm')
        H += g * (a[i]' * sm[i] + a[i] * sm[i]')
        
        # Hopping Interaction: -J(a_i' a_{i+1} + h.c.)
        if i < N_sites
            H -= J * (a[i]' * a[i+1] + a[i] * a[i+1]')
        end
    end
    
    return H, b_total
end

end