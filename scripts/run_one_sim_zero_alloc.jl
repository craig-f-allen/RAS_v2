using RAS, DifferentialEquations, SymbolicIndexingInterface
using Plots, Printf, SteadyStateDiffEq, ModelingToolkit, BenchmarkTools
using NonlinearSolve # For direct Newton solves

# -------------------------------------------------------------
# 1. SETUP PHASE (Script/Global Scope)
# -------------------------------------------------------------
const fract_mut = 0.25
const mut = :G12V

mutants_base = get_mutant_params("C:\\Users\\cfa13\\Projects\\RAS_v2\\RAS\\data\\kinetic_parameter_multipliers.xlsx")
mutants_tri_drug = get_mutant_tri_drug_params("C:\\Users\\cfa13\\Projects\\RAS_v2\\RAS\\data\\kinetic_parameter_multipliers.xlsx")

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

# Setup ODEProblem.
const SS_ABSTOL = 1e-8              # Tolerance for steady state solver. Const for zero alloc.
tspan = (0.0, 1e6)                  # Needed for ODE solver.
odeprob = ODEProblem(sys, param_pairs, tspan)

# Setup integrator. Instead fo solving the ODEProblem we initialize an integrator and use it in place to avoid allocations. DynamicSS cannot do this w/ zero alloc.
integrator = init(odeprob, Rodas5P(); save_everystep=false, save_start=false, dense = false, abstol = 1e-12, reltol = 1e-10)

# Preallocate buffers in memory for saving du. This prevents integrator from allocating memory during the solve.
const n_states = length(integrator.u)  # Number of states in the system. Const for zero alloc.
const du_buf = zeros(n_states)

# Preallocate the drug index.
const drug_idx = variable_index(sys, :Drug)  # Index of the drug variable in the state vector. Const for zero alloc.

# One alllocation steady state solver. This is a custom implementation of DynamicSS that uses the preallocated buffers to avoid allocations during the solve.
function run_ss!(integrator, drug_dose, du_buf)
    #integrator.u[drug_idx] = drug_dose # Set the drug dose in the integrator state vector.
    reinit!(integrator, integrator.u; reinit_dae = false, reinit_cache = false, reset_dt = false, reinit_callbacks = false, reinit_retcode = true) # Update integrator u0.

    # Run the integrator until steady state is reached. This is a custom implementation of DynamicSS that uses the preallocated buffers to avoid allocations during the solve.
    for _ in 1:10_000
        step!(integrator) # take integrator step
        integrator.f(du_buf, integrator.u, integrator.p, integrator.t) # calculate du and store in preallocated buffer
        if maximum(abs, du_buf) < SS_ABSTOL # if du is below tolerance, or approxumately zero, stop.
            break
        end
    end
    return nothing
end


