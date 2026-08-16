# Screens for potential drug resistance mechanisms: which parameters, if shifted within a
# plausible fold-change range, push IC50 up the most. Two disjoint scopes are provided:
#   - MUTANT_KINETIC_PARAMS: the mutant RAS's own rate/Michaelis constants - what a resistance
#     mutation in the mutant allele itself could alter.
#   - ABUNDANCE_PARAMS: total protein/nucleotide pool sizes and the mutant allele fraction -
#     what a resistance mechanism acting through expression/copy-number changes (RAS or
#     effector amplification, GAP/GEF dosage, clonal expansion of the mutant allele) could
#     alter. Neither touches drug-binding chemistry (K_D_1, WT/Mut_K_D_2, CYPA_0, Drug_0).
#
# Both functions below search the same box (each parameter constrained to [theta0/fold,
# theta0*fold] - searched in log-space so the bound is symmetric in fold-change, with fract_mut
# additionally clamped to stay inside (0,1) since it's a fraction, not a free-ranging rate or
# concentration) and reuse the same DoseResponseProblem/drug_dose_response_ic50! pipeline as
# everything else in the package:
#   - ic50_resistance_scan:     one parameter at a time - "which single mechanism"
#   - ic50_resistance_optimize: all of them jointly, via ForwardDiff-gradient-driven Fminbox -
#                                "the best combination", which may compound beyond any one hit
const MUTANT_KINETIC_PARAMS = [:Mut_kint, :Mut_kdGDP, :Mut_kdGTP, :Mut_kaGDP, :Mut_kcat,
                                :Mut_KM, :Mut_kGDP, :Mut_KMGDP, :Mut_KMGTP, :Mut_kGTP,
                                :Mut_KD, :Mut_kaEff]
const ABUNDANCE_PARAMS = [:TotalRAS, :TotalEff, :fract_mut, :GAP, :GEF, :GTP, :GDP]

function _param_log_bounds(sym::Symbol, theta0::Real, fold::Real)
    if sym === :fract_mut
        return log(max(theta0 / fold, 0.01)), log(min(theta0 * fold, 0.99))
    else
        return log(theta0 / fold), log(theta0 * fold)
    end
end

function ic50_resistance_scan(mutant::Symbol, fract_mut::Real,
        param_syms::AbstractVector{Symbol}=MUTANT_KINETIC_PARAMS; fold::Real=10.0, kwargs...)
    base_dict = Dict(build_tricomplex_param_pairs(mutant, fract_mut; kwargs...))
    ic50_0 = drug_dose_response_ic50!(DoseResponseProblem(mutant, fract_mut; kwargs...))

    rows = map(param_syms) do s
        theta0 = base_dict[s]
        lo, hi = _param_log_bounds(s, theta0, fold)
        ic50_at(logtheta) = drug_dose_response_ic50!(
            DoseResponseProblem(mutant, fract_mut; overrides=Dict(s => exp(logtheta)), kwargs...))
        res = Optim.optimize(lt -> -ic50_at(lt), lo, hi)  # 1D bounded, no gradient needed (Brent)

        best_logtheta = Optim.minimizer(res)
        best_ic50 = -Optim.minimum(res)
        (; param=s, baseline=theta0, best_value=exp(best_logtheta),
           best_fold=exp(best_logtheta) / theta0, ic50=best_ic50, ic50_fold_change=best_ic50 / ic50_0)
    end
    return (; ic50_baseline=ic50_0, rows=sort(collect(rows); by=r -> -r.ic50_fold_change))
end

function ic50_resistance_optimize(mutant::Symbol, fract_mut::Real,
        param_syms::AbstractVector{Symbol}=MUTANT_KINETIC_PARAMS; fold::Real=10.0,
        search_n_doses::Int=300, refine_n_doses::Int=1000,
        optim_options=Optim.Options(outer_iterations=6, iterations=25), kwargs...)
    base_dict = Dict(build_tricomplex_param_pairs(mutant, fract_mut; kwargs...))
    theta0 = Float64[base_dict[s] for s in param_syms]
    logtheta0 = log.(theta0)
    bounds = _param_log_bounds.(param_syms, theta0, fold)
    lo = first.(bounds)
    hi = last.(bounds)

    # Search at a coarser dose resolution (gradient evals are the expensive part - each one is
    # a full dose sweep in ForwardDiff.Dual arithmetic), then refine the final answer at full
    # resolution once the optimizer has converged.
    neg_ic50(logtheta; n_doses) = -drug_dose_response_ic50!(
        DoseResponseProblem(mutant, fract_mut; overrides=Dict(param_syms .=> exp.(logtheta)), n_doses, kwargs...))

    res = Optim.optimize(lt -> neg_ic50(lt; n_doses=search_n_doses), lo, hi, logtheta0,
                          Optim.Fminbox(Optim.LBFGS()), optim_options; autodiff=:forward)

    logtheta_opt = Optim.minimizer(res)
    ic50_opt = -neg_ic50(logtheta_opt; n_doses=refine_n_doses)
    ic50_0 = drug_dose_response_ic50!(DoseResponseProblem(mutant, fract_mut; n_doses=refine_n_doses, kwargs...))

    return (; param_syms, theta0, theta_opt=exp.(logtheta_opt),
              fold_change=exp.(logtheta_opt .- logtheta0),
              ic50_baseline=ic50_0, ic50_opt, ic50_fold_change=ic50_opt / ic50_0,
              converged=Optim.converged(res))
end

# Runs ic50_resistance_scan across multiple mutants and assembles the one-at-a-time results
# into a mutants x parameters matrix of best achievable IC50 fold-change - the shape a
# mutant-vs-parameter heatmap needs.
function ic50_resistance_scan_matrix(mutant_names::AbstractVector{Symbol}, fract_mut::Real,
        param_syms::AbstractVector{Symbol}=MUTANT_KINETIC_PARAMS; fold::Real=10.0, kwargs...)
    fold_changes = Matrix{Float64}(undef, length(mutant_names), length(param_syms))
    for (i, mutant) in enumerate(mutant_names)
        scan = ic50_resistance_scan(mutant, fract_mut, param_syms; fold, kwargs...)
        by_param = Dict(r.param => r.ic50_fold_change for r in scan.rows)
        fold_changes[i, :] .= [by_param[s] for s in param_syms]
    end
    return (; mutant_names=collect(mutant_names), param_syms=collect(param_syms), fold_changes)
end
