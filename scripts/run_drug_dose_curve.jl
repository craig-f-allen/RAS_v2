using RAS, DifferentialEquations, SymbolicIndexingInterface
using Plots, Printf, SteadyStateDiffEq, ModelingToolkit, BenchmarkTools
using NonlinearSolve # For direct Newton solves

# -------------------------------------------------------------
# 1. SETUP PHASE (Script/Global Scope)
# -------------------------------------------------------------
const fract_mut = 0.25
const mut = :G12V
const N_DRUG_DOSES = 10000

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
const SS_ABSTOL = 1e-9                  # Tolerance for steady state solver. Const for zero alloc.
tspan = (0.0, 1e6)                      # Needed for ODE solver.
const n_threads = Threads.nthreads()    # Number of threads available for parallel execution. Const for zero alloc.

# Setup vector of prealloc integrators. Using let block to show odeprob is locally scoped, to ensure all integrators are independent to avoid heap alloc and GC.
integrators = [
    let odeprob = ODEProblem(sys, param_pairs, tspan)
        integrator = init(odeprob, Rosenbrock23(); save_everystep=false, save_start=false, dense = false, abstol = 1e-12, reltol = 1e-10)
    end
    for _ in 1:n_threads
]

# Preallocate buffers in memory for saving du. This prevents integrators from allocating memory during the solve.
const n_states = length(integrator.u)  # Number of states in the system. Const for zero alloc.
du_bufs = [zeros(n_states) for _ in 1:n_threads]
u_resets = [copy(integrators[1].u) for _ in 1:n_threads] # Preallocate u resets for each thread to avoid allocations during reinit.

# # Preallocate the drug index.
# const drug_idx = variable_index(sys, :Drug)  # Index of the drug variable in the state vector. Const for zero alloc.

# One alllocation steady state solver. This is a custom implementation of DynamicSS that uses the preallocated buffers to avoid allocations during the solve. Can get to 0 alloc but increases time 1000x
function run_ss!(integrator, drug_dose, du_buf, u_reset)
    u_reset[drug_idx] = drug_dose
    reinit!(integrator, u_reset; reinit_dae = false, reinit_cache = true, reset_dt = true, reinit_callbacks = false, reinit_retcode = true) # Update integrator u0.

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

# Preallocate results vector.
results = [zeros(n_states) for _ in 1:N_DRUG_DOSES] 
drug_doses = exp10.(range(-2, stop=2, length=N_DRUG_DOSES)) # Log spaced drug doses from 10^-2 to 10^2

function ensemble_run_ss!(results, drug_doses, integrators, du_bufs, u_resets)
    n = length(drug_doses)
    nt = length(integrators)
    chunks = Iterators.partition(1:n, cld(n, nt)) # Partition the drug doses into chunks for each thread.
    tasks = map(enumerate(chunks)) do (tid, chunk)
        Threads.@spawn begin
            for i in chunk
                run_ss!(integrators[tid], drug_doses[i], du_bufs[tid], u_resets[tid])
                results[i] .= integrators[tid].u # Save the steady state results for this drug dose.
            end
        end
    end
    foreach(wait, tasks)
end

println("Running ensemble steady state solver...")
ensemble_run_ss!(results, drug_doses, integrators, du_bufs, u_resets) # Run the steady state solver in parallel for all drug doses.

println("Done. Plotting results...")
WT_RAS_GTP_Eff_idx = variable_index(sys, :WT_RAS_GTP_Eff) # Index of the RAS_GTP_Eff_Total variable in the state vector. Const for zero alloc.
Mut_RAS_GTP_Eff_idx = variable_index(sys, :Mut_RAS_GTP_Eff) # Index of the RAS_GTP_Eff_Total variable in the state vector. Const for zero alloc.
ys = [result[WT_RAS_GTP_Eff_idx]+result[Mut_RAS_GTP_Eff_idx] for result in results] # Extract the RAS_GTP_Eff_Total values from the results.