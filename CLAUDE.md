# RAS — Design Philosophy

This is a performance-sensitive scientific computing package (Julia/ModelingToolkit ODE model of
RAS-GTP signaling, plus a drug dose-response/IC50 pipeline). A full dose-response curve solves in
~20ms, threaded, with near-zero heap allocation on the hot path. That performance is a
first-class design constraint, not an afterthought. Read this before touching
`RAS/src/simulation/` or `RAS/src/model/`.

## Core principles (why the code looks the way it does)

1. **Zero/near-zero allocation on the hot path.**
   `run_ss!`, `ensemble_run_ss!`, and `drug_dose_response_ic50!` run per dose, per thread,
   potentially thousands of times per call. Everything they touch is preallocated up front
   (integrators, `du_bufs`, `u_resets`, `results`) and mutated in place. New code on this path
   must not allocate — no new `Dict`s, no comprehensions, no closures that box variables, no
   string formatting, inside these functions or anything they call per-iteration.

2. **Allocation is fine — expected — in setup code.**
   Building a `DoseResponseProblem`/`SingleSimProblem` (param dicts, integrators, xlsx parsing)
   runs once per *problem*, not once per *dose*. Don't over-optimize this path at the cost of
   clarity; the allocation budget that matters is the per-step/per-dose one.

3. **Type stability over generality.**
   `DoseResponseProblem{T,S,Integ,RespIdx}` etc. are parametric so field access stays a
   concrete-type load, not a dynamic dispatch. Module-level caches (`_TRICOMPLEX_SYS_CACHE`) are
   `Ref`s behind `const` bindings — a non-`const` global in Julia is type-unstable and allocates
   on every read. Preserve this; don't introduce non-const globals or loosely-typed containers on
   a path that gets called repeatedly.

4. **One pipeline, not a production path plus a separate AD path.**
   Kinetic/abundance parameters are typed `Real` (not `Float64`) through the constructor chain
   specifically so a caller can seed a `ForwardDiff.Dual` through `overrides` and get a
   differentiable result out of the *exact same* `DoseResponseProblem`/`drug_dose_response_ic50!`
   code used for production dose curves (see `ic50_sensitivity.jl`, `RAS/src/model/build.jl`).
   Don't narrow these back to `Float64` for convenience, and don't add a parallel "AD version" of
   a function — extend the generic one.

5. **Custom steady-state solver, not `DynamicSS`.**
   `run_ss!` is a hand-rolled steady-state loop, written specifically because `DynamicSS` didn't
   hit zero allocations at the tolerances needed (see `RAS/docs/notes.md`). It checks combined
   abstol+reltol convergence per state because state magnitudes span ~15 orders of magnitude.
   Don't swap this for a solver-provided steady-state method without discussing the allocation
   and tolerance tradeoff first.

6. **Threading is manual and explicit.**
   One integrator + buffer set per thread, doses partitioned with `Iterators.partition`,
   `Threads.@spawn` + `wait`. This is deliberate — it keeps allocation and scheduling visible and
   controllable. Don't swap in a higher-level parallel-map abstraction without checking it
   preserves the allocation profile.

7. **No machine-specific paths.**
   Data file paths are derived via `pkgdir(@__MODULE__)` (`DEFAULT_KINETIC_PARAMS_PATH`), never
   hardcoded absolute paths.

8. **Simple, direct code — no speculative abstraction.**
   The codebase favors explicit, slightly repetitive code (e.g. the near-duplicate `RAS_Base`/
   `RAS_Tricomplex` models) over premature generalization. Don't introduce abstraction layers,
   config systems, or "flexibility" that current use cases don't need.

## When to speak up

If a change I'm asking for would work against any of the above — e.g. it adds allocation to
`run_ss!`/`ensemble_run_ss!`/`drug_dose_response_ic50!` or anything in their call chain, breaks
type stability, forks the AD-compatible pipeline into a separate path, replaces the custom
steady-state solver or manual threading, or adds abstraction beyond what's asked — **say so in
chat before making the change**, explain the tradeoff, and let me decide. Don't silently apply a
"cleaner" pattern that costs performance or reintroduces something this codebase deliberately
moved away from.
