# File: jch_solver.jl
using QuantumOptics

# This is the function Python will call
function solve_jch(N, ω, g, J, t_max)
    # 1. Define Bases
    b_cav = FockBasis(N)
    b_atom = SpinBasis(1//2)
    b_site = b_cav ⊗ b_atom
    
    # 2. Define Operators
    a_local = destroy(b_cav) ⊗ one(b_atom)
    sm_local = one(b_cav) ⊗ sigmam(b_atom)
    sz_local = one(b_cav) ⊗ sigmaz(b_atom)
    
    # Expand to 2 sites
    a1 = a_local ⊗ one(b_site)
    sm1 = sm_local ⊗ one(b_site)
    sz1 = sz_local ⊗ one(b_site)
    
    a2 = one(b_site) ⊗ a_local
    sm2 = one(b_site) ⊗ sm_local
    sz2 = one(b_site) ⊗ sz_local
    
    # 3. Hamiltonian
    H = ω*a1'*a1 + 0.5*ω*sz1 + g*(a1'*sm1 + a1*sm1') +
        ω*a2'*a2 + 0.5*ω*sz2 + g*(a2'*sm2 + a2*sm2') - 
        J*(a1'*a2 + a1*a2')
        
    # 4. Initial State: Atom A Excited, Atom B Ground
    psi0 = fockstate(b_cav, 0) ⊗ spinup(b_atom) ⊗ fockstate(b_cav, 0) ⊗ spindown(b_atom)
    
    # 5. Time Evolution
    T = [0:0.5:t_max;]
    tout, psi_t = timeevolution.schroedinger(T, psi0, H)
    
    # 6. Calculate Observables (Atom B Excitation Probability)
    P_exc_op = (sz2 + one(sz2)) / 2
    prob_B = real(expect(P_exc_op, psi_t))
    
    # Return raw data arrays to Python
    return tout, prob_B
end