"""Reproducible scaling benchmarks - results go to CSV, never hard-coded."""

from __future__ import annotations

import csv
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from qutip import basis, tensor

from jch.dynamics import evolve
from jch.hamiltonians import JCParams, cavity_chain, jch_two_site

__all__ = ["BenchmarkResult", "scale_cutoff", "scale_sites", "write_csv"]


@dataclass
class BenchmarkResult:
    backend: str
    axis: str          # what was varied: "n_cutoff" or "n_sites"
    value: int
    dim: int           # Hilbert-space dimension actually solved
    seconds: float


def _timed(fn, repeats: int = 1) -> float:
    """Best-of-N wall time. Best-of, not mean: it suppresses OS scheduling noise."""
    best = float("inf")
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best


def scale_cutoff(
    cutoffs: list[int], t_max: float = 50.0, n_times: int = 200, repeats: int = 1
) -> list[BenchmarkResult]:
    """Fix 2 JCH sites, grow the photon cutoff. Dimension = (2N)^2."""
    tlist = np.linspace(0, t_max, n_times)
    results = []
    for n in cutoffs:
        params = JCParams(n_cutoff=n)
        H, ops = jch_two_site(params, J=0.05)
        psi0 = tensor(basis(n, 0), basis(2, 0), basis(n, 0), basis(2, 1))
        secs = _timed(lambda H=H, psi0=psi0: evolve(H, psi0, tlist), repeats)
        results.append(
            BenchmarkResult("qutip", "n_cutoff", n, H.shape[0], round(secs, 6))
        )
        print(f"  n_cutoff={n:>3}  dim={H.shape[0]:>7}  {secs:.4f}s")
    return results


def scale_sites(
    site_counts: list[int], n_cutoff: int = 3, t_max: float = 50.0, n_times: int = 100
) -> list[BenchmarkResult]:
    """Fix the cutoff, grow the chain. Dimension = n_cutoff ** n_sites.

    This is the run that exposes the exponential wall: each extra site
    multiplies the state space by n_cutoff, and the runtime by roughly the same.
    """
    tlist = np.linspace(0, t_max, n_times)
    results = []
    for n_sites in site_counts:
        H, psi0, n_ops = cavity_chain(n_sites, w_c=1.0, J=0.1, n_cutoff=n_cutoff)
        secs = _timed(lambda H=H, psi0=psi0, n_ops=n_ops: evolve(H, psi0, tlist, e_ops=n_ops))
        results.append(
            BenchmarkResult("qutip", "n_sites", n_sites, H.shape[0], round(secs, 6))
        )
        print(f"  n_sites={n_sites:>3}  dim={H.shape[0]:>9}  {secs:.4f}s")
    return results


def write_csv(results: list[BenchmarkResult], path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(asdict(results[0])))
        writer.writeheader()
        writer.writerows(asdict(r) for r in results)
    return path


def main() -> None:
    """Entry point for the  console script."""
    print("2-site JCH, growing photon cutoff:")
    write_csv(scale_cutoff([3, 5, 8, 10, 12, 15]), "results/data/scaling_cutoff.csv")
    print("Cavity chain, growing site count:")
    write_csv(scale_sites([2, 3, 4, 5, 6, 7, 8]), "results/data/scaling_sites.csv")
    print("Wrote results/data/*.csv")


if __name__ == "__main__":
    main()
