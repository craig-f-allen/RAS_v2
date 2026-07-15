module RAS

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using SciCompDSL

include("model/kinetic_params.jl")
include("model/ras.jl")
include("model/build.jl")

export get_mutant_params
export WT
export RAS_Base
export build_params_for_mtk

end # module RAS
