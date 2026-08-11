import numpy as np
from qutip import basis, coherent, destroy, tensor, thermal_dm

from jch.dynamics import evolve, rabi_period
from jch.hamiltonians import JCParams, jaynes_cummings
from jch.observables import entanglement_entropy, g2_zero


def test_g2_classifies_the_three_light_sources():
    N, a = 20, destroy(20)
    assert np.isclose(g2_zero(basis(N, 1), a), 0.0, atol=1e-9)      # Fock: antibunched
    assert np.isclose(g2_zero(coherent(N, 1.0), a), 1.0, atol=1e-2)  # laser: Poissonian
    assert g2_zero(thermal_dm(N, 1.0), a) > 1.7                      # thermal: bunched


def test_fock_state_g2_follows_1_minus_1_over_n():
    N, a = 25, destroy(25)
    for n in (2, 3, 4):
        assert np.isclose(g2_zero(basis(N, n), a), 1 - 1 / n, atol=1e-9)


def test_entanglement_entropy_zero_for_product_state():
    psi = tensor(basis(10, 0), basis(2, 0))
    assert np.isclose(entanglement_entropy(psi, subsystem=1), 0.0, atol=1e-10)


def test_entropy_peaks_near_ln2_at_quarter_rabi_period():
    """A quarter period makes the maximally entangled (|1,g> + |0,e>)/sqrt(2)."""
    p = JCParams(n_cutoff=10, g=0.05)
    H, _ = jaynes_cummings(p)
    psi0 = tensor(basis(10, 0), basis(2, 0))
    tlist = np.linspace(0, rabi_period(p.g), 200)
    out = evolve(H, psi0, tlist)
    entropies = [entanglement_entropy(s, 1) for s in out.states]
    assert max(entropies) > 0.68           # ln 2 = 0.6931
    assert max(entropies) <= np.log(2) + 1e-9
