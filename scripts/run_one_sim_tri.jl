using RAS, DifferentialEquations, SymbolicIndexingInterface
using Plots, Printf, SteadyStateDiffEq, ModelingToolkit, BenchmarkTools
using NonlinearSolve # For direct Newton solves

# -------------------------------------------------------------
# 1. SETUP PHASE (Script/Global Scope)
# -------------------------------------------------------------
const fract_mut = 0.25
const mut = :G12V

mutants_base = get_mutant_params("/home/craig/Documents/Projects/RAS_v2/RAS/data/kinetic_parameter_multipliers.xlsx")
mutants_tri_drug = get_mutant_tri_drug_params("/home/craig/Documents/Projects/RAS_v2/RAS/data/kinetic_parameter_multipliers.xlsx")

param_dict_base = build_params_for_base_ras(
    WT, mutants_base[mut], 
    6e-11, 2e-10, 18e-6, 180e-6, 4e-7, 4e-7, 
    fract_mut
)

param_dict_tri_drug = build_params_for_ras_tricomplex(
    WT, mutants_tri_drug[:WT], mutants_base[mut], mutants_tri_drug[mut],
    6e-11, 2e-10, 18e-6, 180e-6, 4e-7, 4e-7, 
    fract_mut,
    4e-7, 1e-6, 
)

# Merge into a concrete vector of pairs
param_pairs = Pair[k => v for (k, v) in merge(param_dict_base, param_dict_tri_drug)]

# Build MTK System
@mtkbuild sys = RAS_Tricomplex()

# Create a baseline problem ONCE (this allocates during symbolic compilation)
base_prob = SteadyStateProblem(sys, param_pairs)

function solve_steady_state(prob)
    return solve(prob, DynamicSS(Rosenbrock23()); abstol=1e-12, reltol=1e-10)
end
