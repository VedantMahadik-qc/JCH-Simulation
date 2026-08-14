# =====================================================================
#  FOLLOW-UP CELLS  — paste into scaling_and_multicavity.ipynb
# =====================================================================


# ---------------------------------------------------------------------
# CELL 1 — refit restricted to the asymptotic (dispersive) regime
#   Fit only points where Δ_opt is safely above the crossover, so the
#   dispersive approximation the prediction rests on actually holds.
# ---------------------------------------------------------------------
function fit_exponent(rs, ds; label = "")
    x = log.(rs); y = log.(ds); n = length(x)
    p = (n*sum(x.*y) - sum(x)*sum(y)) / (n*sum(x.^2) - sum(x)^2)
    c = (sum(y) - p*sum(x))/n
    R2 = 1 - sum((y .- (p.*x .+ c)).^2)/sum((y .- sum(y)/n).^2)
    println(rpad(label, 28), "p = ", rpad(round(p, digits=3), 8),
            "R² = ", round(R2, digits=4))
    return p, c
end

println("Fitting over progressively more dispersive subsets:")
for k in 3:length(ratios)
    fit_exponent(ratios[end-k+1:end], Δopts[end-k+1:end];
                 label = "last $k points")
end
println("\n(prediction: p → 0.5 as the fit is restricted to larger κ/γ)")

# cut at a threshold in Δ_opt rather than by point count
thr = 4.0        # only trust points where Δ_opt > 4 g0
sel2 = Δopts .> thr
p2, c2 = fit_exponent(ratios[sel2], Δopts[sel2];
                      label = "Δ_opt > $thr g0")

plot(ratios, Δopts, seriestype=:scatter, m=:circle, ms=6,
     xscale=:log10, yscale=:log10, label="numerics",
     xlabel="κ / γ", ylabel="Δ_opt / g₀",
     title="Optimal detuning: scaling law and its breakdown")
plot!(ratios[sel2], exp(c2).*ratios[sel2].^p2, lw=2,
      label="dispersive fit: p = $(round(p2,digits=3))")
plot!(ratios[sel2], exp(c2).*ratios[sel2].^0.5, ls=:dash, lw=2,
      label="prediction: p = 0.5")
hline!([minimum(Δopts)], ls=:dot, lc=:gray,
       label="saturation ≈ $(round(minimum(Δopts),digits=2)) g₀ (crossover floor)")


# ---------------------------------------------------------------------
# CELL 2 — extend the sweep to larger κ/γ to confirm the asymptote
#   SLOW: gate time grows as Δ, so these are the longest runs.
#   Run overnight or trim the list if short on time.
# ---------------------------------------------------------------------
κ_ext = [3.2, 6.4, 12.8]
println("extending sweep (γ = $γ_fixed):")
for κv in κ_ext
    Δo, Fo, ok = find_Δopt(κv, γ_fixed; Δlo = 5.0, Δhi = 60.0, n = 20)
    if ok
        push!(ratios, κv/γ_fixed); push!(Δopts, Δo)
        push!(Fopts, Fo); push!(interior, true)
    end
    println("  κ = ", rpad(κv,7), "κ/γ = ", rpad(round(κv/γ_fixed,digits=0),9),
            "Δ_opt = ", rpad(round(Δo,digits=3),9), "F = ", rpad(round(Fo,digits=4),9), ok)
end
# then re-run CELL 1 to refit with the new points included


# ---------------------------------------------------------------------
# CELL 3 — two-cavity gate: calibrate the transit velocity
#   The earlier run stopped mid-exchange (0.67/0.33). Scan v to find the
#   FULL exchange (transfer → 1) and the HALF exchange (transfer → 0.5).
#   Effective coupling J_eff ~ g0² J_hop / Δ², so this is much slower
#   than the single-cavity case.
# ---------------------------------------------------------------------
function transfer_two_cav(Δ, Jhop, v; nsteps = 2500)
    _, ρt = run_two_cav(Δ, Jhop; v = v, nsteps = nsteps)
    ρf = ρt[end]
    return real(expect(E2, ρf)), real(maximum(expect(Ntot, ρt)))
end

Δ_2c, Jhop_2c = 8.0, 1.0
vs2 = exp10.(range(log10(0.002), log10(0.05), length = 16))
tr2, nph2 = Float64[], Float64[]
for v in vs2
    t, n = transfer_two_cav(Δ_2c, Jhop_2c, v)
    push!(tr2, t); push!(nph2, n)
end

for (v, t, n) in zip(vs2, tr2, nph2)
    println("v = ", rpad(round(v, digits=4), 9),
            "transfer = ", rpad(round(t, digits=3), 8),
            "peak photons = ", round(n, digits=4))
end

plot(vs2, tr2, m=:circle, lw=2, xscale=:log10, label="transfer to atom 2",
     xlabel="transit velocity v", ylabel="value",
     title="Two-cavity gate: velocity calibration (Δ=$Δ_2c, J_hop=$Jhop_2c)")
plot!(vs2, nph2, m=:square, lw=2, label="peak photon number")
hline!([1.0], ls=:dot,  lc=:black, label="full exchange (iSWAP)")
hline!([0.5], ls=:dash, lc=:black, label="half exchange (√iSWAP)")


# ---------------------------------------------------------------------
# CELL 4 — test the predicted J_eff ∝ g0² J_hop / Δ² scaling
#   At full exchange the required transit time T ∝ 1/J_eff, so the
#   full-exchange velocity should obey  v_full ∝ J_eff ∝ J_hop / Δ².
#   Set v_full below from CELL 3 first.
# ---------------------------------------------------------------------
function v_full_exchange(Δ, Jhop; vlo = 0.001, vhi = 0.08, n = 14)
    vs = exp10.(range(log10(vlo), log10(vhi), length = n))
    ts = [transfer_two_cav(Δ, Jhop, v)[1] for v in vs]
    i  = argmax(ts)                      # velocity giving maximum transfer
    return vs[i], ts[i]
end

println("\nJ_hop dependence at Δ = $Δ_2c:")
Jhops = [0.25, 0.5, 1.0, 2.0, 4.0]
vfs = Float64[]
for J in Jhops
    v, t = v_full_exchange(Δ_2c, J)
    push!(vfs, v)
    println("  J_hop = ", rpad(J,7), "v_full ≈ ", rpad(round(v,digits=4),9),
            "max transfer = ", round(t,digits=3))
end

pJ, cJ = fit_exponent(Jhops, vfs; label = "v_full vs J_hop")
println("(prediction: p = 1, since J_eff ∝ J_hop at fixed Δ)")

plot(Jhops, vfs, m=:circle, lw=2, xscale=:log10, yscale=:log10,
     xlabel="J_hop / g₀", ylabel="full-exchange velocity",
     title="Effective coupling vs inter-cavity hopping", label="numerics")
plot!(Jhops, exp(cJ).*Jhops.^pJ, ls=:dash, lw=2,
      label="fit: p = $(round(pJ,digits=3))")
