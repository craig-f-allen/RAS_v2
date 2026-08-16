"""
RAS-GTP/effector ODE model for wild-type and mutant RAS, including a tricomplex drug-binding
model used for dose-response curves and IC50s. See `DoseResponseProblem`/`SingleSimProblem` for
running the model, and `ic50_sensitivity`/`ic50_resistance_scan` for downstream analysis.
"""
module RAS

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using SciCompDSL
using DifferentialEquations: reinit!, step!, ODEProblem, init, Rodas5P, Rosenbrock23
using SymbolicIndexingInterface: variable_index
using ForwardDiff
using Optim

include("model/kinetic_params.jl")
include("model/ras.jl")
include("model/tricomplex.jl")
include("model/build.jl")
include("simulation/run_ss.jl")
include("simulation/drug_dose_response.jl")
include("simulation/dose_response_problem.jl")
include("simulation/single_sim_problem.jl")
include("simulation/ic50_sensitivity.jl")
include("simulation/ic50_resistance.jl")

export get_mutant_params
export get_mutant_tri_drug_params

export WT
export RAS_Base
export RAS_Tricomplex
export build_params_for_base_ras
export build_params_for_ras_tricomplex
export run_ss!
export ensemble_run_ss!
export find_ic50
export smooth3!
export drug_dose_response_ic50!
export build_tricomplex_system
export build_tricomplex_param_pairs
export DoseResponseProblem
export SingleSimProblem
export ic50_sensitivity
export MUTANT_KINETIC_PARAMS
export ABUNDANCE_PARAMS
export ic50_resistance_scan
export ic50_resistance_optimize
export ic50_resistance_scan_matrix

end # module RAS
