# Transport-mediated two-qubit gates

Two two-level atoms transported through a single-mode cavity. The atom–cavity
coupling is time-dependent, set by each atom's motion through the Gaussian mode:

$$g_i(t) = g_0 \exp\!\big[-(v_i (t - t_{0i})/w)^2\big]$$

so the **transit velocity** `v` is a physical control knob — it sets the pulse area
of the interaction, and therefore the exchange angle.

Read in order; each notebook builds on the previous one's conclusion.

| # | Notebook | Question | Answer |
|---|---|---|---|
| 01 | `01_two_atoms_entanglement.ipynb` | Does sequential transport entangle the atoms? | Yes — Bell state on the `\|eg⟩` input, optimal delay set by the cavity memory ~1/κ |
| 02 | `02_transport_gate_study.ipynb` | Is that protocol a *gate*? | **No.** Structurally asymmetric; `\|ge⟩` and `\|ee⟩` leak entirely into the cavity |
| 03 | `03_gate_characterisation.ipynb` | What gate is the dispersive protocol? | **iSWAP** at full exchange, **√iSWAP** at half exchange |
| 04 | `04_scaling_law_multicavity.ipynb` | How does the optimal detuning scale with loss? | Δ_opt ∝ g₀√(κ/γ), tested against numerics |

## Key results

**The sequential protocol is not a gate.** A two-qubit gate must act unitarily on
the atomic subspace, which requires the cavity to return to vacuum for *every*
input. It doesn't: atom 1 in the ground state does nothing on its half transit, so
atom 2's full transit dumps a photon into an empty cavity with nobody to collect it.
Unitarity error ‖M†M − I‖ ≈ 1.0 — total failure on two of four inputs. This is a
genuine negative result and worth reporting as one.

**The dispersive protocol is a gate.** Detuning far from resonance (Δ ≫ g₀)
suppresses real photon exchange but leaves an effective atom–atom XY coupling
J ≈ g₀²/Δ mediated by a *virtual* photon. Measured cavity population peaks at
≈ 0.015 against the dispersive prediction (g₀/Δ)² = 0.0156 — confirmation the photon
is genuinely virtual.

| protocol | ‖M†M − I‖ | verdict |
|---|---|---|
| sequential (resonant) | ≈ 1.0 | not a gate |
| dispersive, full exchange | ≈ 0.001 | **iSWAP**, F ≈ 0.998 |
| dispersive, half exchange | ≈ 0.0005 | **√iSWAP** — entangling, universal with single-qubit rotations |

**The loss trade-off.** Larger Δ buys protection from cavity loss (cavity population
∝ (g₀/Δ)²) but costs gate duration (t_gate ∝ Δ), which is paid to atomic decay γ.
The optimum moves to larger detuning as the cavity gets leakier:

| κ | optimal Δ | F at optimum |
|---|---|---|
| 0 | (boundary of sweep) | 0.868 |
| 0.05 | ≈ 2.3 | 0.847 |
| 0.20 | ≈ 2.6 | 0.791 |
| 0.50 | ≈ 3.5 | 0.713 |

**The scaling law.** Holding the exchange angle θ fixed (achieved by recalibrating
transit velocity as v ∝ 1/Δ), the infidelity budget is

$$1 - F \approx A\,\gamma\Delta + B\,\kappa/\Delta \quad\Longrightarrow\quad \Delta_{\rm opt} \propto g_0\sqrt{\kappa/\gamma}$$

predicting exponent p = 0.5. Notebook 04 extracts Δ_opt numerically over ~2 decades
of κ/γ and fits p, with boundary maxima rejected as non-optima.

## Status and honest caveats

- **Metric is average basis-state fidelity**, not full process fidelity. An honest
  proxy, but an upgrade to process/gate fidelity (χ-matrix or average gate fidelity
  over the full Haar measure) is the correct next step.
- **The κ = 0 optimum sits at the sweep boundary** and is therefore not a real
  optimum — notebook 04 explicitly rejects boundary maxima for this reason.
- **The two-cavity extension in notebook 04 is exploratory.** The effective coupling
  falls to ~g₀²J_hop/Δ², much weaker than single-cavity, so transits must be
  correspondingly slower. Not yet a confirmed result.
- **`scripts/` holds work not yet folded into the notebooks.** `scaling_law.jl` is
  the earlier standalone version, superseded by notebook 04 (which is self-contained).
  `followup_cells.jl` contains an extended κ/γ sweep and a restricted-regime refit
  that have been written but not run — these are the pending items for confirming
  the p → 0.5 asymptote.

## Relation to the rest of the repository

This thread uses the same coupled-cavity machinery as the JCH lattice work, applied
to a different question. The JCH side asks how photons behave in a lattice; this side
asks whether transport through a cavity can implement a named quantum gate. The
two-cavity extension in notebook 04 is where they meet.
