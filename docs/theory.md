# Theory

> **Skeleton — write this in your own words.** The section headings and the
> prompts under each are a scaffold; the content should come from your literature
> review. Every prompt below is phrased as a question an examiner could plausibly
> ask, so filling them in is viva preparation as much as documentation.
>
> Delete this block when you've written the sections.

## 1. The Jaynes-Cummings model

Start from a two-level atom coupled to a single quantised cavity mode.

- Write down `H_field`, `H_atom`, `H_int` separately before combining them.
- The dipole interaction is `H_int = -d · E`. Show how this becomes
  `ℏg(σ₊ + σ₋)(a + a†)`.
- **The rotating-wave approximation.** Four terms appear; two are kept, two
  discarded. Which two, and what physical process does each discarded term
  describe? Under what condition on `g`, `ω_a`, `ω_c` is dropping them justified —
  and where does that condition fail? (Ultrastrong coupling is the answer to the
  last part, and it is the single likeliest place an examiner will push.)
- Final form: `H_JC = ω_c a†a + (ω_a/2)σ_z + g(a†σ₋ + aσ₊)`.

## 2. Why the JC model is exactly solvable

- Show `[H, N̂] = 0` for `N̂ = a†a + σ₊σ₋`.
- Explain what that conservation buys you: `H` block-diagonalises into 2×2 blocks
  spanned by `{|n, e⟩, |n+1, g⟩}`.
- Diagonalise one block. Recover the dressed states and the `√(n+1)` scaling of the
  Rabi frequency.
- Vacuum Rabi splitting `= 2g` at resonance is a direct consequence — and is
  asserted in `python/tests/test_hamiltonians.py::test_vacuum_rabi_splitting_equals_2g`.

## 3. Collapse and revival

- A coherent state is a superposition over many `|n⟩`, each rotating at
  `Ω_n = 2g√(n+1)`.
- Dephasing of those incommensurate frequencies → collapse. Their partial
  rephasing → revival.
- Estimate the revival time in terms of `g` and `|α|`. Why is this phenomenon
  *evidence for field quantisation* — what would a classical field predict instead?

## 4. From one cavity to a lattice

- Add photon hopping: `-J Σ⟨ij⟩ (a_i† a_j + h.c.)`.
- Justify the minus sign and the nearest-neighbour restriction.
- Define polaritons and state why the JCH model is a photonic analogue of
  Bose-Hubbard.
- **The phase transition.** `J/g` small → Mott-insulating, photons localised;
  `J/g` large → superfluid, photons delocalised. What is the order parameter?

## 5. Open systems

- Motivate the Lindblad master equation: cavities leak, atoms spontaneously emit.
- State it, and identify the collapse operators used here: `√κ a_i`, `√γ σ₋ⁱ`.
- What assumptions does the Lindblad form require (Born, Markov, secular), and when
  do they break?
- Contrast unitary `sesolve` with dissipative `mesolve` — and note the cost jump
  from evolving a state vector of dimension `D` to a density matrix of `D²`.

## 6. The computational wall

This is the section that motivates the entire project, so make it concrete.

- Hilbert-space dimension: `(n_cutoff · 2)^M` for `M` JCH sites.
- Work the arithmetic for a 3×3 lattice at `n_cutoff = 6`: ≈ 5 × 10⁹ amplitudes.
  Convert to memory at 16 bytes per complex double. State plainly why exact
  diagonalisation is impossible there.
- Explain how truncating the Fock space is itself an approximation, and how you
  chose `n_cutoff` — what convergence check justifies it?

## 7. Mean-field theory

- The Gutzwiller ansatz: `|Ψ⟩ ≈ ⊗ᵢ|ψᵢ⟩`.
- Derive `H_eff = H_local − J(a†Σαⱼ + a Σαⱼ*)` with `αⱼ = ⟨aⱼ⟩`.
- Cost drops from exponential to linear in the number of sites.
- **What is thrown away:** inter-site entanglement, by construction. So where is
  mean field quantitative, and where only qualitative? (Superfluid regime versus
  near the Mott transition.)
- Why must the seed state be coherent rather than Fock? A Fock state has `⟨a⟩ = 0`,
  so neighbours feel nothing and the dynamics never start — this is implemented in
  `julia/src/meanfield.jl` and worth being able to explain on sight.

## References

Pull from your literature review bibliography.
