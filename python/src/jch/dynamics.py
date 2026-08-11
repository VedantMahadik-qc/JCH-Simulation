"""Time evolution and steady-state drivers - one place that calls the solver."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from qutip import Qobj, expect, mesolve, steadystate

__all__ = ["evolve", "transmission_spectrum", "rabi_splitting_sweep", "rabi_period"]


def evolve(
    H: Qobj,
    psi0: Qobj,
    tlist: np.ndarray,
    c_ops: Sequence[Qobj] | None = None,
    e_ops: Sequence[Qobj] | None = None,
):
    """Lindblad master-equation evolution.

    With c_ops empty this reduces to the Schrodinger equation; with c_ops
    supplied (e.g. sqrt(kappa) * a for cavity leakage) it is the full open-system
    Lindblad equation

        d(rho)/dt = -i[H, rho] + sum_k ( C_k rho C_k^dag - 1/2 {C_k^dag C_k, rho} )

    Pass e_ops=None to get the full states back (needed for Wigner / entropy).
    """
    return mesolve(H, psi0, tlist, c_ops or [], e_ops=e_ops or [])


def transmission_spectrum(
    build_H,
    w_laser_list: np.ndarray,
    c_ops_fn,
    probe_op_fn,
) -> np.ndarray:
    """Sweep the drive frequency and record the steady-state probe expectation.

    This is numerical spectroscopy: for each laser frequency we solve
    L(rho_ss) = 0 rather than integrating in time, which is far cheaper and
    gives the t -> infinity answer exactly.
    """
    out = []
    for w_L in w_laser_list:
        H, ops = build_H(w_L)
        rho_ss = steadystate(H, c_ops_fn(ops))
        out.append(float(np.real(expect(probe_op_fn(ops), rho_ss))))
    return np.array(out)


def rabi_splitting_sweep(build_H, w_a_list: np.ndarray, n_levels: int = 5) -> np.ndarray:
    """Eigen-energies vs atomic frequency - the avoided-crossing plot.

    At resonance the |1,g> and |0,e> states hybridise into polaritons split by
    2g; away from resonance the branches revert to the bare states.

    Returns an array of shape (len(w_a_list), n_levels).
    """
    rows = []
    for w_a in w_a_list:
        H, _ = build_H(w_a)
        rows.append(H.eigenenergies()[:n_levels])
    return np.array(rows)


def rabi_period(g: float) -> float:
    """Vacuum Rabi period: full excitation swap takes t = pi / g."""
    if g <= 0:
        raise ValueError("g must be positive")
    return float(np.pi / g)
