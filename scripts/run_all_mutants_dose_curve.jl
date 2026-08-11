using RAS, Plots, Printf

const FRACT_MUT = 0.25

# Only sweep mutants present in both the base-kinetics and drug-tricomplex sheets -
# build_tricomplex_param_pairs needs a hit in both dicts.
mutant_names = sort(collect(intersect(keys(get_mutant_params()), keys(get_mutant_tri_drug_params()))))

# WT baseline: fract_mut=0 zeroes out the mutant subpopulation's initial condition, so the
# response curve is pure WT regardless of which mutant's kinetics are loaded alongside it.
wt_prob = DoseResponseProblem(first(mutant_names), 0.0)
ic50_wt = drug_dose_response_ic50!(wt_prob)

ic50s = Vector{Float64}(undef, length(mutant_names))
for (i, mutant) in enumerate(mutant_names)
    prob = DoseResponseProblem(mutant, FRACT_MUT)
    ic50s[i] = drug_dose_response_ic50!(prob)
    @printf("%-10s IC50 = %.4g M (%.3gx WT)\n", mutant, ic50s[i], ic50s[i] / ic50_wt)
end

ratios = ic50s ./ ic50_wt
order = sortperm(ratios)
sorted_names  = string.(mutant_names[order])
sorted_ratios = ratios[order]

bar(sorted_names, sorted_ratios, yscale=:log10, xrotation=60, xtickfontsize=6,
    ylabel="IC50 / WT IC50", legend=false,
    title="Mutant Drug IC50s relative to WT (fract_mut=$FRACT_MUT)")
hline!([1.0], label="WT", color=:black, linestyle=:dash, legend=true)
