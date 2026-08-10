# Setup/construction code. Allocation here is fine and expected - this runs once per problem,
# not once per dose. Function-local scope keeps it type-stable without needing `const` anywhere
# (unlike script/global scope, where non-const bindings are type-unstable and allocate on every
# read - that's the whole reason the old scripts marked everything `const`).

# RAS_Tricomplex() takes no arguments, so @mtkbuild produces the same structure every call -
# cache it once instead of re-running structural_simplify per mutant/problem.
const _TRICOMPLEX_SYS_CACHE = Ref{Any}(nothing)

function build_tricomplex_system(; force_rebuild::Bool=false)
    if force_rebuild || _TRICOMPLEX_SYS_CACHE[] === nothing
        @mtkbuild sys = RAS_Tricomplex()
        _TRICOMPLEX_SYS_CACHE[] = sys
    end
    return _TRICOMPLEX_SYS_CACHE[]
end

# Bundles everything ensemble_run_ss!/drug_dose_response_ic50! need: one integrator/buffer set
# per thread, the dose sweep, and scratch space for the response curve. Parametric on {S,Integ,
# RespIdx} so every field keeps its concrete type for a given instance - that's what keeps field
# access (p.ys, p.integrators, ...) type-stable, the same reason SciML's own ODEProblem/
# ODEIntegrator types are parametric. Downstream calls then allocate only what run_ss!'s step!
# calls do (~1 small allocation per step - see run_ss.jl) - nothing extra from this layer.
struct DoseResponseProblem{S, Integ, RespIdx}
    sys::S
    integrators::Vector{Integ}
    du_bufs::Vector{Vector{Float64}}
    u_resets::Vector{Vector{Float64}}
    results::Vector{Vector{Float64}}
    doses::Vector{Float64}
    ys::Vector{Float64}
    drug_idx::Int
    response_idxs::RespIdx
    ss_abstol::Float64
    ss_reltol::Float64
end

function DoseResponseProblem(mutant::Symbol, fract_mut::Real;
        sys = build_tricomplex_system(),
        GAP::Float64=6e-11, GEF::Float64=2e-10, GDP::Float64=18e-6, GTP::Float64=180e-6,
        TotalRAS::Float64=4e-7, TotalEff::Float64=4e-7, Drug0::Float64=1e2, CYPA::Float64=1e-6,
        dose_range::Tuple{<:Real,<:Real}=(1e-12, 1e4), n_doses::Int=1000,
        response_vars = (:WT_RAS_GTP_Eff, :Mut_RAS_GTP_Eff),
        alg = Rodas5P(), ode_abstol::Float64=1e-8, ode_reltol::Float64=1e-6,
        ss_abstol::Float64=1e-16, ss_reltol::Float64=1e-10,
        tspan::Tuple{Float64,Float64}=(0.0, 1e6), n_threads::Int=Threads.nthreads(),
        xlsx_path::String=DEFAULT_KINETIC_PARAMS_PATH)

    param_pairs = build_tricomplex_param_pairs(mutant, fract_mut;
        GAP, GEF, GDP, GTP, TotalRAS, TotalEff, Drug0, CYPA, xlsx_path)

    # One independent ODEProblem/integrator per thread - let-bound so each is locally scoped,
    # keeping threads from sharing mutable heap state.
    integrators = [
        let odeprob = ODEProblem(sys, param_pairs, tspan)
            init(odeprob, alg; save_everystep=false, save_start=false, dense=false,
                 abstol=ode_abstol, reltol=ode_reltol)
        end
        for _ in 1:n_threads
    ]

    n_states = length(integrators[1].u)
    du_bufs  = [zeros(n_states) for _ in 1:n_threads]
    u_resets = [copy(integrators[1].u) for _ in 1:n_threads]

    lo, hi  = dose_range
    doses   = exp10.(range(log10(Float64(lo)), stop=log10(Float64(hi)), length=n_doses))
    results = [zeros(n_states) for _ in 1:n_doses]
    ys      = zeros(n_doses)

    drug_idx = variable_index(sys, :Drug)
    rv = response_vars isa Symbol ? (response_vars,) : response_vars
    response_idxs = length(rv) == 1 ? variable_index(sys, rv[1]) :
                                       Tuple(variable_index(sys, v) for v in rv)

    return DoseResponseProblem(sys, integrators, du_bufs, u_resets, results, doses, ys,
                                drug_idx, response_idxs, ss_abstol, ss_reltol)
end

# Forwards to the existing drug_dose_response_ic50! method - concrete struct field types mean
# this unpacking is a type-stable field load, not a dynamic dispatch, so the only allocation is
# the ~1 per-step! cost inherited from run_ss! (see run_ss.jl).
drug_dose_response_ic50!(p::DoseResponseProblem) =
    drug_dose_response_ic50!(p.ys, p.doses, p.results, p.response_idxs,
                              p.integrators, p.du_bufs, p.u_resets, p.drug_idx,
                              p.ss_abstol, p.ss_reltol)
