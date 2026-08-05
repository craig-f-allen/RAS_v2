using RAS
using DifferentialEquations
using SymbolicIndexingInterface
using Plots
using Printf
using SteadyStateDiffEq
using ModelingToolkit
using BenchmarkTools
using SciMLStructures

# Set fract_mut
const fract_mut = 0.25
const mut = :G12V

# Get mutant parameters
mutants_base = get_mutant_params("/home/craig/Documents/Projects/RAS_v2/RAS/data/kinetic_parameter_multipliers.xlsx")
mutants_tri_drug = get_mutant_tri_drug_params("/home/craig/Documents/Projects/RAS_v2/RAS/data/kinetic_parameter_multipliers.xlsx")

# Build parameter dictionary
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

param_dict = merge(param_dict_base, param_dict_tri_drug)

# 1. Build MTK System
@mtkbuild sys = RAS_Tricomplex()

# 2. Let MTK compile the system ONCE into a standard ODEProblem
tspan = (0.0, 100.0)
base_ode_prob = ODEProblem(sys, param_dict, tspan)

# 3. Extract the raw numerical function and state vector
raw_f  = base_ode_prob.f.f
u0_raw = copy(base_ode_prob.u0)

# Extract raw Float64 vector from MTKParameters wrapper
p_vec, _ = SciMLStructures.canonicalize(SciMLStructures.Tunable(), base_ode_prob.p)
p_raw = copy(p_vec) # Guaranteed Vector{Float64}

# 4. Get integer array offset for Drug_0
p_idx = variable_index(base_ode_prob, sys.Drug_0)

# 5. Create raw SteadyStateProblem
ss_prob = SteadyStateProblem(raw_f, u0_raw, p_raw)
drug_doses = exp10.(range(-4.0, 2.0, length=100))

# 6. Pre-allocate parameter buffers (Now plain Vector{Float64})
p_buffers = [copy(p_raw) for _ in 1:length(drug_doses)]
for i in 1:length(drug_doses)
    p_buffers[i][p_idx] = drug_doses[i]
end

results = Vector{Float64}(undef, length(drug_doses))
alg = DynamicSS(Rosenbrock23(autodiff=false))

# Warmup run to compile JIT off the clock
sol_warm = solve(remake(ss_prob, p=p_buffers[1]), alg; abstol=1e-12, reltol=1e-10)

# 7. Fast Timed Execution
println("Solving problems...")
@time Threads.@threads for i in 1:length(drug_doses)
    sol = solve(remake(ss_prob, p=p_buffers[i]), alg; abstol=1e-12, reltol=1e-10)
    results[i] = sol[sys.RAS_GTP_Total, end]
end

println("Plotting...")
plot(drug_doses, results; xscale=:log10)