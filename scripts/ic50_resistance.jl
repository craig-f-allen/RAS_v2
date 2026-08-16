# Screens abundance-type parameters (total RAS/effector pools, GAP/GEF dosage, nucleotide
# pool, mutant allele fraction - not protein-specific kinetics, not drug binding) for potential
# resistance mechanisms: which ones, if shifted within a plausible fold-change range, raise
# IC50 the most - one at a time, and jointly. See RAS/src/simulation/ic50_resistance.jl.
using RAS, Printf

const MUTANT = :G12D
const FRACT_MUT = 0.25
const FOLD = 10.0

println("=== One parameter at a time ===")
scan = ic50_resistance_scan(MUTANT, FRACT_MUT, ABUNDANCE_PARAMS; fold=FOLD)
@printf("Baseline IC50 = %.4g M\n\n", scan.ic50_baseline)
for r in scan.rows
    @printf("  %-10s %.4g -> %.4g (%.2fx)   IC50 %.3fx baseline\n",
            r.param, r.baseline, r.best_value, r.best_fold, r.ic50_fold_change)
end

println("\n=== Joint optimization (all params together) ===")
opt = ic50_resistance_optimize(MUTANT, FRACT_MUT, ABUNDANCE_PARAMS; fold=FOLD)
@printf("converged=%s\n", opt.converged)
@printf("IC50 %.4g M -> %.4g M  (%.3fx baseline)\n\n", opt.ic50_baseline, opt.ic50_opt, opt.ic50_fold_change)
for (s, x0, xopt, fc) in zip(opt.param_syms, opt.theta0, opt.theta_opt, opt.fold_change)
    @printf("  %-10s %.4g -> %.4g  (%.2fx)\n", s, x0, xopt, fc)
end
