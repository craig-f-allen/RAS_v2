Before running code, start Julia REPL with the following commands:

```bash
julia -t auto --project=RAS
```

Then in the Julia REPL:

```julia
using Revise
using RAS
```

This will properly load packages and ensure all threads will be used in CPU.