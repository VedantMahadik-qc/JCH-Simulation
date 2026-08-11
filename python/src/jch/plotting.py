"""Plot helpers. Every function takes an Axes so notebooks compose freely."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

__all__ = ["plot_populations", "plot_wigner", "plot_scaling", "plot_lattice_heatmap"]


def plot_populations(tlist, series: dict[str, np.ndarray], ax=None, title: str = ""):
    """Population / photon-number traces against time."""
    ax = ax or plt.subplots(figsize=(8, 5))[1]
    for label, y in series.items():
        ax.plot(tlist, y, label=label, lw=2)
    ax.set_xlabel("Time")
    ax.set_ylabel("Population")
    if title:
        ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.4)
    return ax


def plot_wigner(xvec, W, ax=None, title: str = ""):
    """Wigner function with a symmetric diverging colour map.

    The colour scale is forced symmetric about zero so negative (non-classical)
    lobes are visually unambiguous.
    """
    ax = ax or plt.subplots(figsize=(6, 6))[1]
    lim = float(np.max(np.abs(W)))
    cont = ax.contourf(xvec, xvec, W, 100, cmap="RdBu_r",
                       levels=np.linspace(-lim, lim, 100))
    ax.set_xlabel("Position q")
    ax.set_ylabel("Momentum p")
    ax.set_aspect("equal")
    if title:
        ax.set_title(title)
    return ax, cont


def plot_scaling(x, series: dict[str, np.ndarray], ax=None,
                 xlabel: str = "Hilbert space cutoff (N)", logy: bool = True):
    """Runtime-vs-size comparison across backends."""
    ax = ax or plt.subplots(figsize=(9, 6))[1]
    markers = ["o-", "s-", "^-", "d-"]
    for (label, y), m in zip(series.items(), markers, strict=False):
        y = np.asarray(y, dtype=float)
        y = np.where(y <= 0, 1e-4, y)   # log axis cannot take exact zeros
        ax.plot(x, y, m, label=label, lw=2)
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Execution time (s)")
    ax.legend()
    ax.grid(True, which="both", ls="--", alpha=0.5)
    return ax


def plot_lattice_heatmap(densities: np.ndarray, shape: tuple[int, int], ax=None,
                         title: str = "", vmax: float | None = None):
    """Photon density on an Lx x Ly lattice at a single instant."""
    ax = ax or plt.subplots(figsize=(5, 5))[1]
    grid = np.asarray(densities).reshape(shape)
    im = ax.imshow(grid, cmap="viridis", origin="lower", vmin=0, vmax=vmax)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    if title:
        ax.set_title(title)
    return ax, im
