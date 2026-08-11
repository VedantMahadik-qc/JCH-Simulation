# Benchmark data

Raw timings behind the figures in the project README. Recovered from the original
notebooks, where they had been pasted in as literal Python lists.

## `backend_scaling_chain.csv`

1-D chain of bare coupled cavities, photon cutoff fixed at 3 per site, site count
swept from 2 to 14. Hilbert-space dimension is `n_cutoff ** n_sites`, so the last
row solves 4 782 969 states.

Headline result: at 14 sites QuTiP takes **255.90 s**, QuantumOptics.jl takes
**85.41 s** — a **3.0x** speed-up. QuantumToolbox.jl lands in between at 142.10 s,
trading raw speed for QuTiP-compatible syntax.

Below about 8 sites all three backends finish in well under a second, and the
ordering there is dominated by timer resolution and JIT warm-up rather than by
solver throughput. Several QuantumOptics.jl entries read exactly `0.0000` for this
reason. Do not read anything into the small-N rows; the trend only becomes
meaningful from roughly N = 10 upward, where the runtimes exceed the noise floor
by orders of magnitude.

## `backend_scaling_cutoff.csv`

2-site Jaynes-Cummings-Hubbard dimer, site count fixed at 2, photon cutoff swept
from 3 to 20. Hilbert-space dimension is `(n_cutoff * 2) ** 2`.

These runs stay in the sub-10-millisecond range throughout, which is why the
cutoff sweep is the *weaker* of the two benchmarks: it never leaves the regime
where fixed overheads dominate. It is retained because it isolates the effect of
growing the local Fock space rather than the lattice.

## Caveats

Both datasets were collected on a single machine in a single session, without
repeat trials or a reported hardware specification, so they are indicative rather
than rigorous. `python -m jch.benchmarks` regenerates the QuTiP column
reproducibly, with best-of-N timing and Hilbert-space dimensions recorded per row;
the Julia columns still need an equivalent script before the comparison can be
called properly controlled.
