module RAS

include("model/kinetic_params.jl")
include("model/ras.jl")
include("model/build.jl")

export get_mutant_params
export WT
export RAS_Base
export build_params_for_mtk

end # module RAS
