using RAS
using DifferentialEquations
using SymbolicIndexingInterface
using Plots
using Printf
# using SteadyStateDiffEq

# Set fract_mut
fract_mut = 0.5

# Get mutant parameters
mutants = get_mutant_params("/home/craig/Documents/Projects/RAS_v2/RAS/data/kinetic_parameter_multipliers.xlsx")

# Build WT + Mutant + GEF + GAP + GDP + GTP + TotalRAS + TotalEff + fract_mut param Dict.
param_dict = build_params_for_mtk(WT, mutants[:G13D], 6e-11, 2e-10, 18e-6, 180e-6, 4e-7, 4e-7, fract_mut)

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
    # 1. Rebuild parameter dictionary dynamically for the CURRENT mutant
    p_dict = build_params_for_mtk(WT, mutant_params, 6e-11, 2e-10, 18e-6, 180e-6, 4e-7, 4e-7, fract_mut)
    
    # 2. Extract ordered values using matching keys from keys(param_dict)
    ordered_values = [p_dict[k] for k in keys(param_dict)]
    
    # 3. Create a clean problem instance via remake
    prob = remake(prob_template)
    
    # 4. Safely mutate parameters in-place
    setter!(prob, ordered_values)
    
    # 5. Solve steady state
    sols[mutant_name] = solve(prob, DynamicSS(Rosenbrock23()); abstol=1e-12, reltol=1e-10)
end