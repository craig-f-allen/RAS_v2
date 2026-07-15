using RAS
using DifferentialEquations
using SymbolicIndexingInterface
using Plots
using Printf
# using SteadyStateDiffEq

# Get mutant parameters
mutants = get_mutant_params("/home/craig/Documents/Projects/RAS_v2/RAS/data/kinetic_parameter_multipliers.xlsx")

# Build WT + Mutant + GEF + GAP + GDP + GTP + TotalRAS + TotalEff + fract_mut param Dict.
param_dict = build_params_for_mtk(WT, mutants[:G12C], 6e-11, 2e-10, 18e-6, 180e-6, 4e-7, 4e-7, 0.0)

# Build core RAS mtk model into ODE system
@mtkbuild sys = RAS_Base()

# Map MTK model symbols to param dict keys for efficient lookup.
param_symbols = [getproperty(sys, k) for k in keys(param_dict)]

# Create initial problem. Will take most compilation speed and ready prob_template for remake()
initial_pairs = [k => v for (k, v) in param_dict]
prob_template = SteadyStateProblem(sys, initial_pairs)

# Precompile setter targeting param_symbols order
setter! = setp(sys, param_symbols)

# Storage for solutions
sols = Dict()

# Iteration over mutant parameters, remake makes it safe and efficient.
for (mutant_name, mutant_params) in mutants
    # Rebuild full parameter dictionary for this specific mutant
    p_dict = build_params_for_mtk(WT, mutants[:G12C], 6e-11, 2e-10, 18e-6, 180e-6, 4e-7, 4e-7, 0.0)
    
    # 1. Create shallow copy
    prob = remake(prob_template)
    
    # 2. Extract values strictly ordered by param_symbols keys
    ordered_values = [p_dict[Symbol(p)] for p in param_symbols]
    
    # 3. Mutate parameters in-place safely
    setter!(prob, ordered_values)
    
    # 4. Solve using SteadyState solver
    sols[mutant_name] = solve(prob, DynamicSS(Rosenbrock23()))
end