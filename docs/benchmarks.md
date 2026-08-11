# Benchmark methodology and results

## What is being compared

Three open-quantum-systems libraries solving identical physics:

| Backend | Language | Role |
|---|---|---|
| [QuTiP](https://qutip.org) 5.x | Python | community standard, C-accelerated inner loops via Cython |
| [QuantumOptics.jl](https://qojulia.org) 1.2 | Julia | JIT-compiled, native sparse types |
| [QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl) | Julia | QuTiP-compatible syntax on a Julia core |

The question is narrow and worth stating precisely: **for a fixed physical problem
and a fixed integration tolerance, how does wall-clock solver time grow with
Hilbert-space dimension in each backend?** Not "which library is better" — that
depends on ecosystem, documentation and the cost of rewriting existing code.

## What was measured

Two sweeps, both on a 1-D system with open boundary conditions.

**Chain sweep** (`results/data/backend_scaling_chain.csv`). Bare coupled cavities,
no atoms. Photon cutoff fixed at 3 per site; site count swept 2 → 14. Hilbert-space
dimension `n_cutoff ** n_sites`, reaching 4 782 969 states at N = 14.

**Cutoff sweep** (`results/data/backend_scaling_cutoff.csv`). Two-site JCH dimer
with atoms. Site count fixed at 2; photon cutoff swept 3 → 20. Dimension
`(n_cutoff * 2) ** 2`, reaching 1600 states.

Only the solver call is timed. Hamiltonian construction and operator setup happen
outside the timed region, because they are one-off costs that do not scale with the
integration and would otherwise flatter whichever library builds operators fastest.

Julia runs include a warm-up call before timing, to exclude JIT compilation. This
matters enormously: an un-warmed first call can take seconds regardless of problem
size, which would invert the result at small N.

## Results

At the largest tractable size (N = 14, 4.8 M states):

| Backend | Time | Relative |
|---|---:|---:|
| QuTiP | 255.90 s | 1.00× |
| QuantumToolbox.jl | 142.10 s | 1.80× |
| QuantumOptics.jl | 85.41 s | **3.00×** |

The ordering is stable from N ≈ 10 onward and the gap widens with dimension,
consistent with Python's per-call interpreter overhead being amortised poorly as the
number of solver steps grows.

QuantumToolbox.jl sitting between the two is the interesting result. It is written
in Julia but mirrors QuTiP's API, and it pays for that compatibility — roughly 1.7×
slower than QuantumOptics.jl. The trade is a real one: it is by far the cheapest
migration path for an existing QuTiP codebase.

## Limitations

Stated plainly, because a benchmark without its caveats is an anecdote.

**Small-N rows are noise, not signal.** Below about 8 sites all three backends
finish in under a second. Several QuantumOptics.jl entries read exactly `0.0000` —
the measurement fell below timer resolution. Nothing should be concluded from rows
where the runtime is comparable to the clock granularity.

**Single trial, single machine.** No repeat runs, no error bars, no reported CPU or
RAM. Wall-clock timings on a desktop under normal load can vary by tens of percent.
The 3× headline gap is far larger than that variation, so the conclusion survives —
but a 1.2× difference measured this way would not be meaningful.

**Asymmetric reproducibility.** The QuTiP column can be regenerated with
`python -m jch.benchmarks`, which records Hilbert-space dimension per row and uses
best-of-N timing. The Julia columns are recovered from single-shot runs in
notebooks and have no equivalent script. Until one exists, the comparison is not
properly controlled.

**Solver defaults not pinned.** Each library's default integrator and tolerances
were used rather than a matched (method, rtol, atol) triple. Different adaptive
steppers can take different numbers of steps for the same problem, so part of the
measured gap may reflect solver choice rather than language or implementation.

## Next steps

1. Write `julia/benchmarks/scaling.jl` mirroring `jch.benchmarks`, so both columns
   are reproducible from one command each.
2. Pin the integrator and tolerances explicitly across all three backends.
3. Report hardware, library versions, and best-of-5 timings with spread.
4. Extend the chain sweep past N = 14 for the Julia backends alone, where QuTiP is
   no longer tractable in reasonable time — that is where the argument for Julia
   actually gets made.
