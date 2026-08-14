"""jch - Jaynes-Cummings-Hubbard simulation toolkit (QuTiP backend)."""

from jch.dynamics import (
    evolve,
    rabi_period,
    rabi_splitting_sweep,
    transmission_spectrum,
)
from jch.hamiltonians import (
    JCParams,
    cavity_chain,
    coupled_cavities,
    driven_coupled_cavities,
    hilbert_dim,
    jaynes_cummings,
    jch_two_site,
)
from jch.observables import (
    entanglement_entropy,
    g2_zero,
    photon_numbers,
    wigner_grid,
)

__version__ = "0.1.0"

__all__ = [
    # parameters
    "JCParams",
    # Hamiltonians
    "jaynes_cummings",
    "coupled_cavities",
    "driven_coupled_cavities",
    "jch_two_site",
    "cavity_chain",
    "hilbert_dim",
    # dynamics
    "evolve",
    "rabi_period",
    "rabi_splitting_sweep",
    "transmission_spectrum",
    # observables
    "g2_zero",
    "entanglement_entropy",
    "photon_numbers",
    "wigner_grid",
]
