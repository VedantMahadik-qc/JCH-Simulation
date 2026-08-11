"""Observables: photon statistics, entanglement, phase-space quasi-probability."""

from __future__ import annotations

import numpy as np
from qutip import Qobj, entropy_vn, expect, wigner

__all__ = ["g2_zero", "entanglement_entropy", "photon_numbers", "wigner_grid"]


def g2_zero(state: Qobj, a: Qobj) -> float:
    """Second-order coherence at zero delay.

        g2(0) = <a^dag a^dag a a> / <a^dag a>^2

    Interpretation:
        g2 < 1  antibunched, non-classical (a Fock state |n> gives 1 - 1/n)
        g2 = 1  Poissonian - coherent state (a laser)
        g2 = 2  bunched, thermal / chaotic light
    """
    n_avg = expect(a.dag() * a, state)
    if np.isclose(n_avg, 0.0):
        return 0.0
    return float(np.real(expect(a.dag() * a.dag() * a * a, state) / n_avg**2))


def entanglement_entropy(state: Qobj, subsystem: int = 1, base: str = "e") -> float:
    """Von Neumann entropy of one subsystem of a bipartite pure state.

        S = -Tr(rho_A ln rho_A)

    For a pure global state this measures atom-field entanglement. The maximum
    for a two-level subsystem is ln(2) ~ 0.693.
    """
    rho = state.ptrace(subsystem)
    return float(entropy_vn(rho, base=np.e if base == "e" else 2))


def photon_numbers(states: list[Qobj], n_ops: list[Qobj]) -> np.ndarray:
    """Photon number at every site, for every time. Shape (n_times, n_sites)."""
    return np.array([[float(np.real(expect(op, s))) for op in n_ops] for s in states])


def wigner_grid(
    state: Qobj, subsystem: int | None = None, span: float = 5.0, points: int = 100
):
    """Wigner quasi-probability of a (possibly reduced) cavity state.

    Negative regions are the signature of non-classicality; they are what makes
    the Fock and Schrodinger-cat plots in this project interesting.

    Returns (xvec, W).
    """
    rho = state if subsystem is None else state.ptrace(subsystem)
    xvec = np.linspace(-span, span, points)
    return xvec, wigner(rho, xvec, xvec)
