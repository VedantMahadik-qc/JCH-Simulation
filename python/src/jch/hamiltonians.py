"""Hamiltonian builders for cavity-QED and Jaynes-Cummings-Hubbard systems.

Conventions used throughout (hbar = 1):

    Single site (Jaynes-Cummings), Hilbert space  Fock(N) (x) Spin(1/2):

        H = w_c a^dag a  +  (w_a/2) sigma_z  +  g (a^dag sigma_- + a sigma_+)

    Hopping between neighbouring cavities carries a MINUS sign, matching the
    Bose-Hubbard convention used in the dissertation:

        H_hop = -J (a_i^dag a_j + a_i a_j^dag)
"""

from __future__ import annotations

from dataclasses import dataclass

from qutip import Qobj, basis, destroy, num, qeye, sigmam, sigmaz, tensor

__all__ = [
    "JCParams",
    "jaynes_cummings",
    "coupled_cavities",
    "driven_coupled_cavities",
    "jch_two_site",
    "cavity_chain",
]


@dataclass(frozen=True)
class JCParams:
    """Parameters of a Jaynes-Cummings site.

    Attributes
    ----------
    n_cutoff : Fock-space truncation. The cavity keeps states |0> .. |n_cutoff-1>.
    w_c      : cavity frequency.
    w_a      : atomic transition frequency.
    g        : atom-cavity (dipole) coupling.
    """

    n_cutoff: int = 15
    w_c: float = 1.0
    w_a: float = 1.0
    g: float = 0.05

    def __post_init__(self) -> None:
        if self.n_cutoff < 2:
            raise ValueError(f"n_cutoff must be >= 2, got {self.n_cutoff}")
        if self.g < 0:
            raise ValueError(f"g must be non-negative, got {self.g}")

    @property
    def detuning(self) -> float:
        """Atom-cavity detuning  Delta = w_a - w_c."""
        return self.w_a - self.w_c


def jaynes_cummings(params: JCParams) -> tuple[Qobj, dict[str, Qobj]]:
    """Build the single-site Jaynes-Cummings Hamiltonian.

    Returns
    -------
    H   : the Hamiltonian on Fock(n_cutoff) (x) Spin(1/2).
    ops : dict of the operators used, so callers never rebuild them by hand.
          Keys: "a", "sm", "sz", "n".
    """
    n = params.n_cutoff
    a = tensor(destroy(n), qeye(2))
    sm = tensor(qeye(n), sigmam())
    sz = tensor(qeye(n), sigmaz())

    H = (
        params.w_c * a.dag() * a
        + 0.5 * params.w_a * sz
        + params.g * (a.dag() * sm + a * sm.dag())
    )
    return H, {"a": a, "sm": sm, "sz": sz, "n": a.dag() * a}


def coupled_cavities(
    n_cutoff: int, w_a: float, w_b: float, J: float
) -> tuple[Qobj, dict[str, Qobj]]:
    """Two bare cavities (no atoms) exchanging photons at rate J.

        H = w_a a^dag a + w_b b^dag b + J (a^dag b + a b^dag)

    Note the PLUS sign here: this is the beam-splitter convention, which is why
    the normal modes sit at w_c -/+ J and a single photon swaps with period
    T = pi / J.
    """
    if n_cutoff < 2:
        raise ValueError(f"n_cutoff must be >= 2, got {n_cutoff}")
    a = tensor(destroy(n_cutoff), qeye(n_cutoff))
    b = tensor(qeye(n_cutoff), destroy(n_cutoff))

    H = w_a * a.dag() * a + w_b * b.dag() * b + J * (a.dag() * b + a * b.dag())
    ops = {
        "a": a,
        "b": b,
        "n_a": tensor(num(n_cutoff), qeye(n_cutoff)),
        "n_b": tensor(qeye(n_cutoff), num(n_cutoff)),
    }
    return H, ops


def driven_coupled_cavities(
    n_cutoff: int, w_c: float, J: float, F: float, w_laser: float
) -> tuple[Qobj, dict[str, Qobj]]:
    """Two coupled cavities, cavity A driven by a laser, in the rotating frame.

    Moving to the frame rotating at the laser frequency replaces w_c by the
    detuning  delta = w_c - w_laser  and makes the drive time-independent:

        H = delta (a^dag a + b^dag b) + J (a^dag b + a b^dag) + F (a^dag + a)
    """
    delta = w_c - w_laser
    a = tensor(destroy(n_cutoff), qeye(n_cutoff))
    b = tensor(qeye(n_cutoff), destroy(n_cutoff))

    H = (
        delta * a.dag() * a
        + delta * b.dag() * b
        + J * (a.dag() * b + a * b.dag())
        + F * (a.dag() + a)
    )
    return H, {"a": a, "b": b}


def jch_two_site(
    params: JCParams, J: float
) -> tuple[Qobj, dict[str, Qobj]]:
    """Two-site Jaynes-Cummings-Hubbard dimer.

    Tensor ordering is (cavity_1, atom_1, cavity_2, atom_2), so the total
    Hilbert-space dimension is (n_cutoff * 2)^2.
    """
    n = params.n_cutoff
    Ic, Ia = qeye(n), qeye(2)

    a1 = tensor(destroy(n), Ia, Ic, Ia)
    sm1 = tensor(Ic, sigmam(), Ic, Ia)
    sz1 = tensor(Ic, sigmaz(), Ic, Ia)
    a2 = tensor(Ic, Ia, destroy(n), Ia)
    sm2 = tensor(Ic, Ia, Ic, sigmam())
    sz2 = tensor(Ic, Ia, Ic, sigmaz())

    def site(a, sm, sz):
        return (
            params.w_c * a.dag() * a
            + 0.5 * params.w_a * sz
            + params.g * (a.dag() * sm + a * sm.dag())
        )

    H = site(a1, sm1, sz1) + site(a2, sm2, sz2) - J * (a1.dag() * a2 + a1 * a2.dag())

    ops = {
        "a1": a1, "a2": a2, "sm1": sm1, "sm2": sm2, "sz1": sz1, "sz2": sz2,
        "n1": a1.dag() * a1, "n2": a2.dag() * a2,
        "p_exc_1": (sz1 + 1) / 2, "p_exc_2": (sz2 + 1) / 2,
    }
    return H, ops


def cavity_chain(
    n_sites: int, w_c: float, J: float, n_cutoff: int
) -> tuple[Qobj, Qobj, list[Qobj]]:
    """1-D chain of `n_sites` bare cavities with open boundary conditions.

        H = w_c sum_i n_i  -  J sum_<ij> (a_i^dag a_j + h.c.)

    Returns (H, psi0, n_ops) where psi0 has one photon on site 0 and n_ops[i]
    is the photon-number operator on site i.

    Warning: the Hilbert-space dimension is n_cutoff ** n_sites. This is the
    exponential wall - n_sites=14 at n_cutoff=3 is already 4.8 million states.
    """
    if n_sites < 1:
        raise ValueError(f"n_sites must be >= 1, got {n_sites}")

    ident = qeye(n_cutoff)
    a = destroy(n_cutoff)

    def embed(op, site):
        ops = [ident] * n_sites
        ops[site] = op
        return tensor(ops)

    H = sum(w_c * embed(a.dag() * a, i) for i in range(n_sites))
    for i in range(n_sites - 1):
        H -= J * (
            embed(a.dag(), i) * embed(a, i + 1) + embed(a, i) * embed(a.dag(), i + 1)
        )

    psi = [basis(n_cutoff, 0)] * n_sites
    psi[0] = basis(n_cutoff, 1)
    psi0 = tensor(psi)

    n_ops = [embed(num(n_cutoff), i) for i in range(n_sites)]
    return H, psi0, n_ops


def hilbert_dim(n_sites: int, n_cutoff: int, with_atoms: bool = True) -> int:
    """Dimension of the JCH Hilbert space - the number quoted in the write-up."""
    per_site = n_cutoff * (2 if with_atoms else 1)
    return per_site ** n_sites
