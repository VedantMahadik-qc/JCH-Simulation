# Jaynes–Cummings–Hubbard Simulation

[![python](https://github.com/VedantMahadik-qc/JCH-Simulation/actions/workflows/python.yml/badge.svg)](https://github.com/VedantMahadik-qc/JCH-Simulation/actions/workflows/python.yml)
[![julia](https://github.com/VedantMahadik-qc/JCH-Simulation/actions/workflows/julia.yml/badge.svg)](https://github.com/VedantMahadik-qc/JCH-Simulation/actions/workflows/julia.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Numerical study of the **Jaynes–Cummings–Hubbard (JCH)** model — an array of coupled
cavity-QED sites where light–matter coupling competes with photon hopping — with a
cross-language performance benchmark of **QuTiP (Python)**, **QuantumOptics.jl** and
**QuantumToolbox.jl (Julia)**, plus a CUDA path for larger lattices.

MSc Quantum Technologies project, University College London.

---

## The model

With ħ = 1, on a lattice of sites *i* with nearest-neighbour bonds ⟨ij⟩:

$$
H \;=\; \sum_i \Big[\, \omega_c\, a_i^\dagger a_i \;+\; \tfrac{\omega_a}{2}\,\sigma_z^i \;+\; g\,(a_i^\dagger \sigma_-^i + a_i \sigma_+^i)\,\Big] \;-\; J \sum_{\langle ij\rangle} \big(a_i^\dagger a_j + \mathrm{h.c.}\big)
$$

| Symbol | Meaning |
|---|---|
| `ω_c` | cavity mode frequency |
| `ω_a` | atomic transition frequency |
| `g` | atom–cavity dipole coupling (Rabi frequency) |
| `J` | photon hopping amplitude between neighbouring cavities |
| `κ`, `γ` | cavity leakage / atomic spontaneous emission (Lindblad) |

Open-system dynamics use the Lindblad master equation

$$
\dot\rho = -i[H,\rho] + \sum_k \Big( C_k \rho C_k^\dagger - \tfrac12\{C_k^\dagger C_k, \rho\}\Big),
\qquad C_k \in \{\sqrt{\kappa}\,a_i,\ \sqrt{\gamma}\,\sigma_-^i\}.
$$

**The central difficulty.** A lattice of *M* sites truncated at *N* photons per cavity
lives in a Hilbert space of dimension `(N·2)^M`. A 3×3 grid with `N = 6` is already
≈ 5×10⁹ complex amplitudes — beyond exact diagonalisation. Everything in this
repository is organised around that wall: exact solvers below it, mean-field above it.

---

## Results

### Backend scaling — 1-D cavity chain, cutoff 3 photons/site

| Sites | Hilbert dim | QuTiP (s) | QuantumOptics.jl (s) | QuantumToolbox.jl (s) |
|---:|---:|---:|---:|---:|
| 10 | 59 049 | 3.18 | 0.82 | 1.46 |
| 12 | 531 441 | 30.30 | 10.08 | 15.00 |
| 14 | 4 782 969 | 255.90 | **85.41** | 142.10 |

**≈ 3.0× speed-up** for QuantumOptics.jl over QuTiP at the largest tractable size,
with the gap widening as the Hilbert space grows. QuantumToolbox.jl sits between the
two — it trades some raw speed for QuTiP-compatible syntax.

Raw numbers: [`results/data/backend_scaling_chain.csv`](results/data/backend_scaling_chain.csv).
The QuTiP column regenerates reproducibly with `python -m jch.benchmarks`.

**Read the small-N rows with care.** Below roughly 8 sites every backend finishes in
under a second, and the ordering there reflects timer resolution and JIT warm-up
rather than solver throughput — several QuantumOptics.jl entries read exactly
`0.0000` for that reason. The comparison only becomes meaningful from about N = 10,
where runtimes clear the noise floor by orders of magnitude. Full methodology and
limitations in [`docs/benchmarks.md`](docs/benchmarks.md).

### Physics reproduced

| Phenomenon | Where |
|---|---|
| Vacuum Rabi oscillations; splitting = 2g at resonance | `notebooks/01_jc_single_site/` |
| Collapse and revival under a coherent drive | `notebooks/01_jc_single_site/` |
| Schrödinger-cat formation (negative Wigner lobes) | `notebooks/01_jc_single_site/` |
| Atom–field entanglement entropy → ln 2 | `notebooks/01_jc_single_site/` |
| g⁽²⁾(0) for Fock / coherent / thermal light | `notebooks/01_jc_single_site/` |
| Photon hopping and normal-mode splitting at ω_c ± J | `notebooks/02_coupled_cavities/` |
| Two-site JCH state transfer | `notebooks/02_coupled_cavities/` |
| 2-D wavepacket spreading on a 3×3 lattice | `notebooks/03_lattice/` |
| Gutzwiller mean-field dynamics, 5×5 lattice | `notebooks/03_lattice/` |
| GPU (CUDA) vs CPU master-equation solve | `notebooks/04_gpu/` |

---

## Repository layout

```
JCH-Simulation/
├── python/                   pip-installable `jch` package (QuTiP backend)
│   ├── src/jch/
│   │   ├── hamiltonians.py     JC, coupled cavities, JCH dimer, cavity chain
│   │   ├── dynamics.py         mesolve wrapper, steady-state spectroscopy
│   │   ├── observables.py      g2(0), von Neumann entropy, Wigner grids
│   │   ├── benchmarks.py       timing runs → CSV (no hard-coded numbers)
│   │   └── plotting.py         reusable figure helpers
│   └── tests/                  21 tests — physics assertions, not smoke tests
├── julia/                    JCHSimulation.jl
│   ├── src/
│   │   ├── hamiltonians.jl     1-D chain, 2-D lattice, neighbour tables
│   │   └── meanfield.jl        Gutzwiller product-state solver
│   └── test/runtests.jl        20 tests across 7 test sets
├── notebooks/                thin demos — they import, they do not define
│   ├── 01_jc_single_site/      Rabi, collapse/revival, Wigner, entropy, g2
│   ├── 02_coupled_cavities/    hopping, spectroscopy, tomography, benchmarks
│   ├── 03_lattice/             2-D spreading, mean field
│   ├── 04_gpu/                 CUDA vs CPU
│   ├── julia_quantumoptics/    same physics, QuantumOptics.jl
│   └── julia_quantumtoolbox/   same physics, QuantumToolbox.jl
├── results/
│   ├── data/                   benchmark CSVs + provenance notes (tracked)
│   ├── figures/                showcase PNGs tracked, scratch output ignored
│   └── animations/             showcase GIFs tracked, scratch output ignored
└── docs/
    ├── theory.md               derivation from cavity QED first principles
    └── benchmarks.md           methodology and full results
```

---

## Quick start

### Python

```bash
cd python
python -m pip install -e ".[dev,notebooks]"
pytest -q                       # 21 tests
```

```python
import numpy as np
from jch import JCParams, jaynes_cummings, evolve
from jch.observables import entanglement_entropy

params = JCParams(n_cutoff=10, w_c=1.0, w_a=1.0, g=0.05)
H, ops = jaynes_cummings(params)

from qutip import basis, tensor
psi0 = tensor(basis(10, 0), basis(2, 0))          # |0 photons, atom excited⟩
result = evolve(H, psi0, np.linspace(0, 200, 500), e_ops=[ops["n"]])
```

### Julia

```bash
julia --project=julia -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

```julia
using JCHSimulation
p = JCHParams(n_cutoff=3, g=0.05, J=0.1)
H, basis, a, n = build_jch_chain(p, 4)

sol, b_site, ops, _ = meanfield_lattice(p, 5, 5)   # 5×5, mean-field
n_t = photon_numbers(sol, b_site, ops, 25)
```

---

## Methods

**Exact diagonalisation / master equation** for ≤ ~10⁵ states. Sparse operators
throughout; the JCH Hamiltonian is block-diagonal in total excitation number, which
QuTiP and QuantumOptics.jl both exploit.

**Gutzwiller mean-field** above that. The many-body state is approximated as a product
`⊗ᵢ |ψᵢ⟩` and sites couple only through `αⱼ = ⟨aⱼ⟩`, giving an effective single-site
Hamiltonian `H_eff = H_local − J(a†Σαⱼ + a Σαⱼ*)`. Cost becomes linear in the number of
sites. Inter-site entanglement is discarded by construction, so results are quantitative
in the superfluid regime and qualitative near the Mott transition. A coherent seed state
is required — a Fock state has `⟨a⟩ = 0` and would never couple.

**GPU offload** via `CUDA.jl` + `QuantumToolbox.jl`: matrices are moved with `cu(...)`
and the same `mesolve` call runs on device. Worth it only once the Hilbert space is
large enough to amortise the transfer.

---

## Citation

If this code is useful to you:

```bibtex
@software{mahadik_jch_simulation,
  author  = {Mahadik, Vedant},
  title   = {JCH-Simulation: benchmarking Python and Julia for
             Jaynes--Cummings--Hubbard dynamics},
  year    = {2026},
  url     = {https://github.com/VedantMahadik-qc/JCH-Simulation}
}
```

## Licence

MIT — see [LICENSE](LICENSE).
