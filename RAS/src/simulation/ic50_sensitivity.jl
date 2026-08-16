# Sensitivity of the IC50 output to named model parameters, via forward-mode AD. Reuses the
# exact same DoseResponseProblem/drug_dose_response_ic50! pipeline used for production dose
# curves - T is inferred from whatever ForwardDiff seeds through `overrides`, so there is no
# separate/duplicated "AD version" of the IC50 pipeline to keep in sync.
#
# `param_syms` defaults to every parameter build_tricomplex_param_pairs exposes. Pass a subset
# to only differentiate w.r.t. those (cheaper - fewer ForwardDiff chunks).
function ic50_sensitivity(mutant::Symbol, fract_mut::Real,
        param_syms::Union{Nothing,AbstractVector{Symbol}}=nothing;
        fd_check::Int=5, kwargs...)
    base_pairs = build_tricomplex_param_pairs(mutant, fract_mut; kwargs...)
    syms = param_syms === nothing ? first.(base_pairs) : collect(param_syms)
    base_dict = Dict(base_pairs)
    theta0 = Float64[base_dict[s] for s in syms]

    f(theta) = drug_dose_response_ic50!(
        DoseResponseProblem(mutant, fract_mut; overrides=Dict(syms .=> theta), kwargs...))

    ic50_0 = f(theta0)
    grad = ForwardDiff.gradient(f, theta0)
    elasticities = grad .* theta0 ./ ic50_0

    # find_ic50's linear interpolation between dose-grid points has a kink at the half-max
    # crossing - spot-check the largest-magnitude elasticities against central finite
    # differences as a cheap guard against a parameter landing exactly on one.
    fd_relerr = Dict{Symbol,Float64}()
    if fd_check > 0
        order = sortperm(abs.(elasticities), rev=true)
        for i in order[1:min(fd_check, end)]
            h = theta0[i] * 1e-4
            thetap = copy(theta0); thetap[i] += h
            thetam = copy(theta0); thetam[i] -= h
            fd = (f(thetap) - f(thetam)) / (2h)
            fd_relerr[syms[i]] = abs(grad[i] - fd) / max(abs(fd), eps())
        end
    end

    return (; ic50=ic50_0, param_syms=syms, theta=theta0, grad, elasticities, fd_relerr)
end
