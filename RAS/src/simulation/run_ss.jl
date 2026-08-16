# One allocation steady state solver. This is a custom implementation of DynamicSS that uses the preallocated buffers to avoid allocations during the solve.
# Convergence uses a combined absolute+relative tolerance per state (like a solver's own error control), because
# state magnitudes here span ~15 orders of magnitude (1e-15 to 1e2 M) - a pure absolute tolerance is either far too
# loose for the large-magnitude states (declaring "steady state" mid-transient) or far too tight for the small ones.
function run_ss!(integrator, du_buf, u0, ss_abstol, ss_reltol)
    reinit!(integrator, u0; reinit_dae = false, reinit_cache = true, reset_dt = true, reinit_callbacks = false, reinit_retcode = true) # Update integrator u0.

    # Run the integrator until steady state is reached.
    for _ in 1:10_000
        step!(integrator) # take integrator step
        integrator.f(du_buf, integrator.u, integrator.p, integrator.t) # calculate du and store in preallocated buffer
        converged = true
        @inbounds for i in eachindex(du_buf)
            if abs(du_buf[i]) >= ss_abstol + ss_reltol * abs(integrator.u[i])
                converged = false
                break
            end
        end
        if converged
            break
        end
    end
    return nothing
end