using RAS, DifferentialEquations, SymbolicIndexingInterface
using Plots, Printf, SteadyStateDiffEq, ModelingToolkit, BenchmarkTools
using NonlinearSolve # For direct Newton solves

# -------------------------------------------------------------
# 1. SETUP PHASE (Script/Global Scope)
# -------------------------------------------------------------
const fract_mut = 0.25
const mut = :G12V
const N_DRUG_DOSES = 1000

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
    1e2, 1e-6,
)

# Merge into a concrete vector of pairs
param_pairs = Pair[k => v for (k, v) in merge(param_dict_base, param_dict_tri_drug)]

# Build MTK System
@mtkbuild sys = RAS_Tricomplex()

# Setup ODEProblem.
const SS_ABSTOL = 1e-16                 # Absolute tolerance for steady state solver. Const for zero alloc.
const SS_RELTOL = 1e-10                 # Relative tolerance for steady state solver. Const for zero alloc.
tspan = (0.0, 1e6)                      # Needed for ODE solver.
const n_threads = Threads.nthreads()    # Number of threads available for parallel execution. Const for zero alloc.

# Setup vector of prealloc integrators. Using let block to show odeprob is locally scoped, to ensure all integrators are independent to avoid heap alloc and GC.
integrators = [
    let odeprob = ODEProblem(sys, param_pairs, tspan)
        integrator = init(odeprob, Rodas5P(); save_everystep=false, save_start=false, dense = false, abstol = 1e-8, reltol = 1e-6)
    end
    for _ in 1:n_threads
]

# Preallocate buffers in memory for saving du. This prevents integrators from allocating memory during the solve.
const n_states = length(integrators[1].u)   # Number of states in the system. Const for zero alloc.
const drug_idx = variable_index(sys, :Drug) # Index of the drug variable in the state vector. Const for zero alloc.
du_bufs = [zeros(n_states) for _ in 1:n_threads]
u_resets = [copy(integrators[1].u) for _ in 1:n_threads] # Preallocate u resets for each thread to avoid allocations during reinit.

# Preallocate results vector.
results = [zeros(n_states) for _ in 1:N_DRUG_DOSES]
drug_doses = exp10.(range(-12, stop=4, length=N_DRUG_DOSES)) # Log spaced drug doses from 10^-12 to 10^4

function ensemble_run_ss!(results, drug_doses, integrators, du_bufs, u_resets, drug_idx)
    n = length(drug_doses)
    nt = length(integrators)
    chunks = Iterators.partition(1:n, cld(n, nt)) # Partition the drug doses into chunks for each thread.
    tasks = map(enumerate(chunks)) do (tid, chunk)
        Threads.@spawn begin
            for i in chunk
                u_resets[tid][drug_idx] = drug_doses[i] # Set this dose's drug concentration before reinit.
                run_ss!(integrators[tid], du_bufs[tid], u_resets[tid], SS_ABSTOL, SS_RELTOL)
                results[i] .= integrators[tid].u # Save the steady state results for this drug dose.
            end
        end
    end
    foreach(wait, tasks)
end

println("Running ensemble steady state solver...")
ensemble_run_ss!(results, drug_doses, integrators, du_bufs, u_resets, drug_idx) # Run the steady state solver in parallel for all drug doses.

println("Done. Plotting results...")
WT_RAS_GTP_Eff_idx = variable_index(sys, :WT_RAS_GTP_Eff) # Index of the RAS_GTP_Eff_Total variable in the state vector. Const for zero alloc.
Mut_RAS_GTP_Eff_idx = variable_index(sys, :Mut_RAS_GTP_Eff) # Index of the RAS_GTP_Eff_Total variable in the state vector. Const for zero alloc.
ys = [result[WT_RAS_GTP_Eff_idx]+result[Mut_RAS_GTP_Eff_idx] for result in results] # Extract the RAS_GTP_Eff_Total values from the results.
plot(drug_doses, ys, xscale=:log10, xlabel="Drug Dose (M)", ylabel="RAS_GTP_Eff_Total", title="Drug Dose Response Curve", legend=false)


