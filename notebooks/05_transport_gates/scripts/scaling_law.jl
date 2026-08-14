# =====================================================================
#  Scaling law for the optimal detuning of a transport-mediated
#  dispersive iSWAP gate.
#
#  PREDICTION (derive this analytically in the report):
#    At fixed exchange angle θ:
#      - cavity-loss penalty  = κ θ / Δ        (shape-independent!)
#      - atomic-decay penalty = γ T ∝ γ θ Δ / g0²
#    Total infidelity ≈ A γ Δ + B κ / Δ
#    => Δ_opt ∝ g0 sqrt(κ/γ)      [square-root scaling]
#
#  TEST: extract Δ_opt numerically over ~2 decades of κ/γ and fit the
#  exponent. Expect slope → 0.5 in the dispersive regime, and deviation
#  in the resonant–dispersive crossover (small κ/γ, where Δ_opt ~ g0).
#
#  Run AFTER the cells in transport_gate_study.ipynb that define
#  b_cav, b_at, a, s1m, s2m, exc1, exc2, Hc1, Hc2, gpulse, labels,
#  full_ket, g0, w.
# =====================================================================

using LinearAlgebra, Plots
gr()

# ---------------------------------------------------------------------
# average basis-state fidelity for the full-exchange (iSWAP) gate,
# with the transit velocity recalibrated at each Δ so the exchange
# angle θ is held FIXED as Δ varies (this is what makes the scaling
# argument apply).
# ---------------------------------------------------------------------
function gate_fid(Δ, κv, γv; nsteps = 1500)
    v  = g0^2 * w * sqrt(2/pi) / Δ          # fixed θ  =>  v ∝ 1/Δ
    t0 = 4*(w/v)
    T  = range(0, 8*(w/v), length = nsteps)
    H0 = Δ*(exc1 + exc2)
    Jl = AbstractOperator[]
    κv > 0 && push!(Jl, sqrt(κv)*a)
    γv > 0 && push!(Jl, sqrt(γv)*s1m)
    γv > 0 && push!(Jl, sqrt(γv)*s2m)
    Jd = dagger.(Jl)
    f(t, ρ) = (H0 + gpulse(t, t0, v)*(Hc1 + Hc2), Jl, Jd)
    ideal = ["gg", "eg", "ge", "ee"]        # iSWAP action in populations
    F = 0.0
    for (j, si) in enumerate(labels)
        _, ρt = timeevolution.master_dynamic(T, full_ket(si), f)
        ψid = full_ket(ideal[j])
        F += real(dagger(ψid) * (ρt[end] * ψid))
    end
    return F/4
end

# ---------------------------------------------------------------------
# locate Δ_opt by coarse scan + parabolic refinement about the peak
# ---------------------------------------------------------------------
function find_Δopt(κv, γv; Δlo = 1.5, Δhi = 40.0, n = 22)
    Δs = range(Δlo, Δhi, length = n)
    Fs = [gate_fid(Δ, κv, γv) for Δ in Δs]
    i  = argmax(Fs)
    # reject boundary maxima: the optimum must be interior to be meaningful
    if i == 1 || i == n
        return (Δs[i], Fs[i], false)
    end
    # parabolic interpolation through the three points around the peak
    x1,x2,x3 = Δs[i-1], Δs[i], Δs[i+1]
    y1,y2,y3 = Fs[i-1], Fs[i], Fs[i+1]
    denom = (y1 - 2y2 + y3)
    Δopt  = denom == 0 ? x2 : x2 - 0.5*(x3-x1)*(y3-y1)/(4*denom)
    return (Δopt, y2, true)
end

# ---------------------------------------------------------------------
# MAIN SWEEP: vary κ/γ over ~2 decades
# ---------------------------------------------------------------------
γ_fixed = 0.002                       # small enough that Δ_opt moves out
                                      # into the genuinely dispersive regime
κ_list  = [0.02, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6]

ratios, Δopts, Fopts, interior = Float64[], Float64[], Float64[], Bool[]
println("γ = ", γ_fixed)
println(" κ      κ/γ      Δ_opt     F_opt    interior?")
for κv in κ_list
    Δo, Fo, ok = find_Δopt(κv, γ_fixed)
    push!(ratios, κv/γ_fixed); push!(Δopts, Δo)
    push!(Fopts, Fo); push!(interior, ok)
    println(rpad(κv,7), rpad(round(κv/γ_fixed,digits=1),9),
            rpad(round(Δo,digits=3),10), rpad(round(Fo,digits=4),9), ok)
end

# ---------------------------------------------------------------------
# FIT THE EXPONENT:  Δ_opt ∝ (κ/γ)^p    ->  expect p ≈ 0.5
# (fit only the interior optima; boundary hits are not real optima)
# ---------------------------------------------------------------------
sel = interior
x = log.(ratios[sel]); y = log.(Δopts[sel])
n = length(x)
p = (n*sum(x.*y) - sum(x)*sum(y)) / (n*sum(x.^2) - sum(x)^2)
c = (sum(y) - p*sum(x))/n
ss_res = sum((y .- (p.*x .+ c)).^2)
ss_tot = sum((y .- sum(y)/n).^2)
R2 = 1 - ss_res/ss_tot

println("\n--- FIT ---")
println("fitted exponent p = ", round(p, digits=3), "   (prediction: 0.5)")
println("R² = ", round(R2, digits=4))
println(abs(p-0.5) < 0.08 ? "=> consistent with sqrt(κ/γ) scaling" :
        "=> deviates from sqrt scaling — check whether Δ_opt is still in the crossover")

plot(ratios[sel], Δopts[sel], seriestype=:scatter, m=:circle, ms=6,
     xscale=:log10, yscale=:log10, label="numerics",
     xlabel="κ / γ", ylabel="Δ_opt / g₀",
     title="Scaling of the optimal detuning")
plot!(ratios[sel], exp(c).*ratios[sel].^p, lw=2,
      label="fit: slope = $(round(p,digits=3))")
plot!(ratios[sel], exp(c).*ratios[sel].^0.5, ls=:dash, lw=2,
      label="prediction: slope = 0.5")
