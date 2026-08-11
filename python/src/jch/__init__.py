"""jch - Jaynes-Cummings-Hubbard simulation toolkit (QuTiP backend)."""

from jch.dynamics import evolve, rabi_splitting_sweep, transmission_spectrum
from jch.hamiltonians import (
    cavity_chain,
    coupled_cavities,
    driven_coupled_cavities,
    jaynes_cummings,
    jch_two_site,
)
from jch.observables import entanglement_entropy, g2_zero, photon_numbers

__version__ = "0.1.0"

__all__ = [
    "jaynes_cummings",
    "coupled_cavities",
    "driven_coupled_cavities",
    "jch_two_site",
    "cavity_chain",
    "g2_zero",
    "entanglement_entropy",
    "photon_numbers",
    "evolve",
    "transmission_spectrum",
    "rabi_splitting_sweep",
]
