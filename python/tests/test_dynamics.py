import numpy as np
from qutip import basis, tensor

from jch.dynamics import evolve, rabi_period
from jch.hamiltonians import JCParams, coupled_cavities, jaynes_cummings


def test_vacuum_rabi_full_swap_at_half_period():
    """Atom excited, cavity empty. After t = pi/(2g) the photon is in the cavity."""
    p = JCParams(n_cutoff=10, g=0.05)
    H, ops = jaynes_cummings(p)
    psi0 = tensor(basis(10, 0), basis(2, 0))     # |0, e>
    t_half = rabi_period(p.g) / 2
    out = evolve(H, psi0, np.linspace(0, t_half, 200), e_ops=[ops["n"]])
    assert out.expect[0][-1] > 0.99               # one photon now in the cavity


def test_photon_swap_between_bare_cavities():
    """A single photon hops A -> B in time pi/(2J)."""
    J = 0.05
    H, ops = coupled_cavities(6, 1.0, 1.0, J)
    psi0 = tensor(basis(6, 1), basis(6, 0))
    tlist = np.linspace(0, np.pi / (2 * J), 200)
    out = evolve(H, psi0, tlist, e_ops=[ops["n_a"], ops["n_b"]])
    assert out.expect[0][-1] < 0.01
    assert out.expect[1][-1] > 0.99


def test_closed_evolution_conserves_total_excitation():
    p = JCParams(n_cutoff=12, g=0.05)
    H, ops = jaynes_cummings(p)
    psi0 = tensor(basis(12, 0), basis(2, 0))
    n_exc = ops["n"] + 0.5 * (ops["sz"] + 1)
    out = evolve(H, psi0, np.linspace(0, 100, 50), e_ops=[n_exc])
    assert np.allclose(out.expect[0], 1.0, atol=1e-6)


def test_cavity_decay_drains_the_photon():
    """Adding a collapse operator sqrt(kappa) a must make <n> decay to zero."""
    p = JCParams(n_cutoff=12, g=0.05)
    H, ops = jaynes_cummings(p)
    psi0 = tensor(basis(12, 0), basis(2, 0))
    c_ops = [np.sqrt(0.05) * ops["a"]]
    out = evolve(H, psi0, np.linspace(0, 500, 200), c_ops=c_ops, e_ops=[ops["n"]])
    assert out.expect[0][-1] < 1e-3
