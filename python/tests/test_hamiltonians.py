"""Physics-level tests. Each one asserts something a viva examiner could ask."""

import numpy as np
import pytest
from qutip import basis, tensor

from jch.hamiltonians import (
    JCParams,
    cavity_chain,
    coupled_cavities,
    hilbert_dim,
    jaynes_cummings,
    jch_two_site,
)


def test_jc_is_hermitian():
    H, _ = jaynes_cummings(JCParams(n_cutoff=10))
    assert H.isherm


def test_jc_dimension_is_2N():
    H, _ = jaynes_cummings(JCParams(n_cutoff=7))
    assert H.shape == (14, 14)


def test_jc_vacuum_energy_is_minus_half_wa():
    """|0,g> is an exact eigenstate: no photon to absorb, no excitation to emit."""
    p = JCParams(n_cutoff=10, w_c=1.0, w_a=1.0, g=0.05)
    H, _ = jaynes_cummings(p)
    ground = tensor(basis(10, 0), basis(2, 1))
    assert np.isclose(H.matrix_element(ground, ground), -0.5 * p.w_a)


def test_vacuum_rabi_splitting_equals_2g():
    """On resonance the first excitation manifold splits by exactly 2g."""
    g = 0.05
    H, _ = jaynes_cummings(JCParams(n_cutoff=20, w_c=1.0, w_a=1.0, g=g))
    e = np.sort(H.eigenenergies())
    assert np.isclose(e[2] - e[1], 2 * g, atol=1e-9)


def test_jc_conserves_excitation_number():
    """[H, N_exc] = 0 with N_exc = a^dag a + sigma^+ sigma^-. This is why the
    JC model is analytically solvable - it block-diagonalises."""
    p = JCParams(n_cutoff=12)
    H, ops = jaynes_cummings(p)
    n_exc = ops["n"] + 0.5 * (ops["sz"] + 1)
    assert np.allclose((H * n_exc - n_exc * H).full(), 0, atol=1e-12)


def test_coupled_cavities_normal_modes_split_by_2J():
    """Two degenerate cavities hybridise into bonding / anti-bonding at w -/+ J."""
    J = 0.05
    H, _ = coupled_cavities(6, 1.0, 1.0, J)
    e = np.sort(H.eigenenergies())
    assert np.isclose(e[2] - e[1], 2 * J, atol=1e-9)


def test_jch_dimer_dimension_and_hermiticity():
    H, ops = jch_two_site(JCParams(n_cutoff=3), J=0.05)
    assert H.shape == (36, 36)          # (3 * 2) ** 2
    assert H.isherm


def test_chain_dimension_matches_formula():
    H, psi0, n_ops = cavity_chain(4, w_c=1.0, J=0.1, n_cutoff=3)
    assert H.shape == (81, 81)          # 3 ** 4
    assert len(n_ops) == 4
    assert psi0.shape == (81, 1)


def test_chain_starts_with_exactly_one_photon_on_site_zero():
    from qutip import expect
    _, psi0, n_ops = cavity_chain(4, 1.0, 0.1, 3)
    occupancies = [float(np.real(expect(op, psi0))) for op in n_ops]
    assert np.allclose(occupancies, [1, 0, 0, 0])


@pytest.mark.parametrize("bad", [0, 1, -3])
def test_cutoff_below_two_is_rejected(bad):
    with pytest.raises(ValueError):
        JCParams(n_cutoff=bad)


def test_hilbert_dim_matches_the_number_in_the_writeup():
    """3x3 lattice, 5 photons max: (5+1)^9 * 2^9 ~ 5e9 - the exponential wall."""
    assert hilbert_dim(9, 6, with_atoms=True) == 6**9 * 2**9
