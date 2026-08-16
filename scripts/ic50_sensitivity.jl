# IC50 sensitivity analysis - see RAS/src/simulation/ic50_sensitivity.jl for the pipeline
# (same DoseResponseProblem used by the other dose-curve scripts, differentiated via ForwardDiff).
using RAS, Plots, Printf

const MUTANT = :G12D
const FRACT_MUT = 0.25

result = ic50_sensitivity(MUTANT, FRACT_MUT)

@printf("IC50(%s, fract_mut=%g) = %.4g M\n", MUTANT, FRACT_MUT, result.ic50)

println("\nFD spot-check relative errors (guards against a param landing on a curve kink):")
for (s, e) in result.fd_relerr
    @printf("  %-16s relerr=%.2g%s\n", s, e, e > 1e-2 ? "  <-- check" : "")
end

order = sortperm(abs.(result.elasticities), rev=true)
println("\nElasticities dlog(IC50)/dlog(theta), sorted by magnitude:")
for i in order
    @printf("  %-22s %+.4f\n", result.param_syms[i], result.elasticities[i])
end

sorted_syms  = string.(result.param_syms[order])
sorted_elast = result.elasticities[order]
bar(sorted_syms[end:-1:1], sorted_elast[end:-1:1], orientation=:h,
    xlabel="Elasticity  dlog(IC50)/dlog(theta)", legend=false, size=(800, 900),
    ytickfontsize=6, title="IC50 sensitivity - $MUTANT (fract_mut=$FRACT_MUT)")
vline!([0.0], color=:black, linewidth=1)
