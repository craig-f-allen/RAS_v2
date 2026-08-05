When using this the first time, make sure to install RAS package properly. It is highly reccomended to understand package activation and instantiation and the REPL before running this code. Code is optimized for speed of simulation, Julia has first compilation time that can be quite long. This is a tradeoff, so for larger simulations we beat Python or other more dynamic language speeds. A single drug dose curve costs 20ms compared to 1 second with SciPY on a 10 core, 20 threaded AMD Ryzen AI 9 card. 

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