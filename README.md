# RAS

An ODE model of RAS-GTP/effector signaling, for wild-type RAS and RAS mutants, including a
tricomplex drug-binding model (drug + CYPA + RAS-GTP) used to compute dose-response curves and
IC50s. Also includes IC50 sensitivity analysis (via automatic differentiation) and resistance
screening (which parameter shifts raise IC50 the most).

A single drug dose-response curve solves in ~20ms, roughly 50x faster than an equivalent SciPy
implementation, by reusing preallocated integrators across the dose sweep instead of building a
new ODE problem per dose.

This is a [Julia](https://julialang.org/) package. If you've never used Julia before, see
[New to Julia?](#new-to-julia) below before continuing.

## Requirements

- Julia (developed and tested with 1.12; see [downloads](https://julialang.org/downloads/))
- Git, to clone this repository

## First-time setup

1. Clone this repository and open a terminal in it.

2. Start Julia with the `RAS/` folder activated as the project environment, and threading
   enabled (the dose-sweep code uses multiple threads):

   ```bash
   julia -t auto --project=RAS
   ```

   On macOS/Linux you can also run `./start_julia.sh`, which runs the same command.

3. Install the package's dependencies. In the Julia REPL, press `]` to enter [package mode](https://docs.julialang.org/en/v1/stdlib/Pkg/#Getting-Started)
   (the prompt changes to `(RAS) pkg>`), then run:

   ```
   instantiate
   ```

   This installs the exact package versions recorded in `RAS/Manifest.toml`. It only needs to be
   run once (or again after `RAS/Project.toml`/`RAS/Manifest.toml` change). Press Backspace to
   return to the normal `julia>` prompt.

4. Load the package:

   ```julia
   using Revise
   using RAS
   ```

   [`Revise`](https://timholy.github.io/Revise.jl/stable/) watches the source files and reloads
   your edits into the running session automatically, so you don't need to restart the REPL while
   developing. Always `using Revise` before `using RAS` if you plan to edit the code.

Repeat steps 2 and 4 (`julia -t auto --project=RAS`, then `using Revise; using RAS`) every time
you start a new session. Step 3 (`instantiate`) is only needed once per machine.

## Running the example scripts

`scripts/` has runnable examples of the package's API. With the package loaded as above, run one
from the REPL with `include`:

```julia
include("scripts/run_one_sim.jl")
```

- `run_one_sim.jl` — run a single mutant to steady state at one drug dose.
- `run_drug_dose_curve.jl` — full dose-response curve and IC50 for one mutant; plots the curve.
- `run_all_mutants_dose_curve.jl` — IC50 for every mutant in the dataset, relative to WT; bar plot.
- `ic50_sensitivity.jl` — how sensitive a mutant's IC50 is to each model parameter.
- `ic50_resistance.jl` — screens abundance parameters (RAS/effector levels, GAP/GEF dosage, ...)
  for the biggest potential IC50 shift, one at a time and jointly.
- `ic50_resistance_heatmap.jl` — the same resistance screen across every mutant, as a heatmap.

## Basic usage

```julia
using RAS

# Run one mutant to steady state at a fixed drug dose.
prob = SingleSimProblem(:G12V, 0.25)   # mutant, fraction of RAS pool that is mutant
u_ss = run_ss!(prob, 1e5)              # steady state at drug dose = 1e5

# Full dose-response curve and IC50.
prob = DoseResponseProblem(:G12V, 0.25)
ic50 = drug_dose_response_ic50!(prob)  # prob.doses / prob.ys hold the curve afterward
```

## Repository layout

```
RAS/
  src/
    model/        # ModelingToolkit ODE system definitions and kinetic parameter loading
    simulation/    # steady-state solving, dose-response sweeps, IC50 finding/sensitivity/resistance
  data/            # kinetic parameters (xlsx)
scripts/           # runnable examples of the package API (see above)
```

## New to Julia?

If you're coming from another language, the official
[Getting Started](https://docs.julialang.org/en/v1/manual/getting-started/) guide covers running
the REPL and basic syntax, and the [Pkg manual](https://docs.julialang.org/en/v1/stdlib/Pkg/)
covers environments and package installation (the `]`-prompt used in step 3 above). A `--project=X`
flag tells Julia to use `X`'s `Project.toml`/`Manifest.toml` as the active environment instead of
your global one, which is why all commands above are run with `--project=RAS`.
