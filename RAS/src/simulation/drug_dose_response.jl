# Sweeps drug dose across the preallocated integrators/threads (same pattern as run_ss!) and stores
# the steady-state vector for each dose in `results`. Threaded; this function's own logic allocates
# nothing beyond what run_ss! does (1 small allocation per integrator step! call - see run_ss.jl),
# plus the usual Task overhead from Threads.@spawn.
function ensemble_run_ss!(results, doses, integrators, du_bufs, u_resets, drug_idx, ss_abstol, ss_reltol)
    n = length(doses)
    nt = length(integrators)
    chunks = Iterators.partition(1:n, cld(n, nt)) # Partition the drug doses into chunks for each thread.
    tasks = map(enumerate(chunks)) do (tid, chunk)
        Threads.@spawn begin
            for i in chunk
                u_resets[tid][drug_idx] = doses[i] # Set this dose's drug concentration before reinit.
                run_ss!(integrators[tid], du_bufs[tid], u_resets[tid], ss_abstol, ss_reltol)
                results[i] .= integrators[tid].u # Save the steady state result for this drug dose.
            end
        end
    end
    foreach(wait, tasks)
    return results
end

# In-place 3-point moving average, single pass, O(1) extra memory (no second buffer needed - the trick
# is carrying the pre-overwrite value of ys[i-1] forward in `prev`). Steady-state solver tolerance leaves
# a little non-monotonic jitter between adjacent doses; this damps it out before IC50 interpolation so a
# spurious local wiggle doesn't get mistaken for the half-max crossing. Endpoints are left untouched since
# find_ic50 anchors y_top/y_bot on them and they already sit on the curve's flat plateaus.
function smooth3!(ys::AbstractVector{<:Real})
    n = length(ys)
    n < 3 && return ys
    prev = ys[1]
    @inbounds for i in 2:n-1
        cur = ys[i]
        ys[i] = (prev + cur + ys[i+1]) / 3
        prev = cur
    end
    return ys
end

# IC50 finder for dose-response curves. Assumes the response is monotonic between its first and last
# points (true of a smoothed, well-formed sigmoidal dose-response curve) - lets us do a single linear
# scan for the half-max crossing instead of a real root-find. Doses don't need to be ascending - only
# consecutive entries need to be adjacent doses, which holds whichever direction the sweep was built in.
# Interpolates in log10(dose) space since doses are log-spaced and the curve is ~sigmoidal there, so a
# straight-line interpolation between the two bracketing points is a good, cheap approximation.

# Core routine: single pass, no allocations, no bounds checks in the hot loop.
function find_ic50(doses::AbstractVector{<:Real}, ys::AbstractVector{<:Real})
    n = length(doses)
    y_top = ys[1]
    y_bot = ys[n]
    half = (y_top + y_bot) / 2
    decreasing = y_bot < y_top # inhibitors: response falls as dose rises; handle the opposite too.

    @inbounds for i in 1:n-1
        y1 = ys[i]
        y2 = ys[i+1]
        crossed = decreasing ? (y1 >= half >= y2) : (y1 <= half <= y2)
        if crossed
            y2 == y1 && return doses[i] # flat segment straddling half; avoid 0/0.
            frac = (half - y1) / (y2 - y1)
            log_d1 = log10(doses[i])
            log_d2 = log10(doses[i+1])
            return exp10(log_d1 + frac * (log_d2 - log_d1))
        end
    end
    return doses[n] * NaN # curve never crosses half-max (not monotonic / bad data).
end

@inline response_at(result, idx::Int) = @inbounds result[idx]
@inline function response_at(result, idxs::NTuple{N,Int}) where {N}
    # Explicit loop instead of sum(f, idxs) - the closure form boxes `result` (abstract
    # element type) and allocates; this doesn't.
    s = zero(eltype(result))
    @inbounds for i in idxs
        s += result[i]
    end
    return s
end

# Convenience overload for the `results` vector-of-state-vectors produced by ensemble_run_ss! -
# reads the response (single state index, or a tuple of indices to sum, e.g. WT+Mut RAS-GTP-Eff)
# straight out of each state vector instead of materializing a separate ys array first.
function find_ic50(doses::AbstractVector{<:Real}, results::AbstractVector{<:AbstractVector{<:Real}}, idxs)
    n = length(doses)
    y_top = response_at(results[1], idxs)
    y_bot = response_at(results[n], idxs)
    half = (y_top + y_bot) / 2
    decreasing = y_bot < y_top

    @inbounds for i in 1:n-1
        y1 = response_at(results[i], idxs)
        y2 = response_at(results[i+1], idxs)
        crossed = decreasing ? (y1 >= half >= y2) : (y1 <= half <= y2)
        if crossed
            y2 == y1 && return doses[i]
            frac = (half - y1) / (y2 - y1)
            log_d1 = log10(doses[i])
            log_d2 = log10(doses[i+1])
            return exp10(log_d1 + frac * (log_d2 - log_d1))
        end
    end
    return doses[n] * NaN
end

# End-to-end pipeline: sweep the doses (creating the curve via ensemble_run_ss!), extract the response
# (single state index, or a tuple of indices to sum) from each resulting state vector into the caller's
# `ys` buffer, smooth it, then interpolate the IC50. `ys` is reused as scratch space throughout, so
# nothing here allocates beyond ensemble_run_ss!'s own per-step! cost (see run_ss.jl) - the
# extraction/smoothing/interpolation stages themselves are genuinely zero allocation.
function drug_dose_response_ic50!(ys::AbstractVector{<:Real}, doses::AbstractVector{<:Real},
                                   results::AbstractVector{<:AbstractVector{<:Real}}, idxs,
                                   integrators, du_bufs, u_resets, drug_idx, ss_abstol, ss_reltol)
    ensemble_run_ss!(results, doses, integrators, du_bufs, u_resets, drug_idx, ss_abstol, ss_reltol)
    @inbounds for i in eachindex(ys)
        ys[i] = response_at(results[i], idxs)
    end
    smooth3!(ys)
    return find_ic50(doses, ys)
end
