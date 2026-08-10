# Single-integrator counterpart to DoseResponseProblem - same setup-once/mutate-in-place shape,
# just one integrator instead of a per-thread vector, for the "run one dose, get steady state" case.
struct SingleSimProblem{S, Integ}
    sys::S
    integrator::Integ
    du_buf::Vector{Float64}
    u_reset::Vector{Float64}
    drug_idx::Int
    ss_abstol::Float64
    ss_reltol::Float64
end

function SingleSimProblem(mutant::Symbol, fract_mut::Real;
        sys = build_tricomplex_system(),
        GAP::Float64=6e-11, GEF::Float64=2e-10, GDP::Float64=18e-6, GTP::Float64=180e-6,
        TotalRAS::Float64=4e-7, TotalEff::Float64=4e-7, Drug0::Float64=1e2, CYPA::Float64=1e-6,
        alg = Rosenbrock23(), ode_abstol::Float64=1e-12, ode_reltol::Float64=1e-10,
        ss_abstol::Float64=1e-16, ss_reltol::Float64=1e-10,
        tspan::Tuple{Float64,Float64}=(0.0, 1e6),
        xlsx_path::String=DEFAULT_KINETIC_PARAMS_PATH)

    param_pairs = build_tricomplex_param_pairs(mutant, fract_mut;
        GAP, GEF, GDP, GTP, TotalRAS, TotalEff, Drug0, CYPA, xlsx_path)

    odeprob = ODEProblem(sys, param_pairs, tspan)
    integrator = init(odeprob, alg; save_everystep=false, save_start=false, dense=false,
                       abstol=ode_abstol, reltol=ode_reltol)
    n_states = length(integrator.u)
    du_buf  = zeros(n_states)
    u_reset = copy(integrator.u)
    drug_idx = variable_index(sys, :Drug)

    return SingleSimProblem(sys, integrator, du_buf, u_reset, drug_idx, ss_abstol, ss_reltol)
end

# Extends the already-exported run_ss! via dispatch, same pattern as drug_dose_response_ic50!(::DoseResponseProblem).
function run_ss!(p::SingleSimProblem, drug_dose::Real)
    p.u_reset[p.drug_idx] = drug_dose
    run_ss!(p.integrator, p.du_buf, p.u_reset, p.ss_abstol, p.ss_reltol)
    return p.integrator.u
end
