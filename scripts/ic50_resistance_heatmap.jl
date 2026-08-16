# Resistance-potential heatmap: for every mutant we have kinetics for, and every abundance
# parameter (ABUNDANCE_PARAMS - total RAS/effector pools, GAP/GEF dosage, nucleotide pool,
# mutant allele fraction), the best achievable IC50 fold-change within a fold-change box,
# one parameter at a time. See RAS/src/simulation/ic50_resistance.jl.
using RAS, Printf

const FRACT_MUT = 0.25
const FOLD = 10.0

mutant_names = sort(collect(intersect(keys(get_mutant_params()), keys(get_mutant_tri_drug_params()))))
println("Mutants: ", mutant_names)

result = ic50_resistance_scan_matrix(mutant_names, FRACT_MUT, ABUNDANCE_PARAMS; fold=FOLD)

println("\nIC50 fold-change matrix (rows=mutant, cols=parameter):")
print(rpad("", 12))
for s in result.param_syms
    print(rpad(string(s), 11))
end
println()
for (i, m) in enumerate(result.mutant_names)
    print(rpad(string(m), 12))
    for j in eachindex(result.param_syms)
        @printf("%-11.4f", result.fold_changes[i, j])
    end
    println()
end
