module RAS

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using SciCompDSL

include("model/kinetic_params.jl")
include("model/ras.jl")
include("model/tricomplex.jl")
include("model/build.jl")

export get_mutant_params
export get_mutant_tri_drug_params

export WT
export RAS_Base
export RAS_Tricomplex
export build_params_for_base_ras
export build_params_for_ras_tricomplex

end # module RAS
