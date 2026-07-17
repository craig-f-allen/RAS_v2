using RAS
using DifferentialEquations
using SymbolicIndexingInterface
using Plots
using Printf
using SteadyStateDiffEq
using ModelingToolkit

# Set fract_mut
fract_mut = 0.25

# Get mutant parameters
mutants = get_mutant_params("/home/craig/Documents/Projects/RAS_v2/RAS/data/kinetic_parameter_multipliers.xlsx")

# Build parameter dictionary
param_dict = build_params_for_mtk(
    WT, mutants[:G13D], 
    6e-11, 2e-10, 18e-6, 180e-6, 4e-7, 4e-7, 
    fract_mut
)

# Build core RAS ModelingToolkit model into an ODE system
@mtkbuild sys = RAS_Base()

# Pass param_dict directly to the ODEProblem definition
# Note: Ensure tspan starts at 0.0 unless your system specifically requires t=1.0 as initial time
# tspan = (0.0, 10000.0)
# prob = ODEProblem(sys, param_dict, tspan)

initial_pairs = [k => v for (k, v) in param_dict]
prob = SteadyStateProblem(sys, initial_pairs)

# Solve the ODE system
sol = solve(prob, DynamicSS(Rosenbrock23()); abstol=1-12, reltol=1e-10)
#sol = solve(prob, Rosenbrock23())