include("kinetic_params.jl")
include("ras.jl")


function build_params_for_mtk(WT::RAS.WTKineticParams, Mut::RAS.MutantKineticParams, GAP::Float64, GEF::Float64, GDP::Float64, GTP::Float64, TotalRAS::Float64, TotalEff::Float64, fract_mut::Float64)::Dict{Symbol, Float64}

    ps = Dict{Symbol, Float64}()

    ps[:TotalRAS] = TotalRAS
    ps[:TotalEff] = TotalEff
    ps[:fract_mut] = fract_mut

    ps[:GAP] = GAP
    ps[:GEF] = GEF
    ps[:GDP] = GDP
    ps[:GTP] = GTP
    
    # WT params
    ps[:WT_kint]   = WT.kint
    ps[:WT_kdGDP]  = WT.kdGDP
    ps[:WT_kdGTP]  = WT.kdGTP
    ps[:WT_kaGDP]  = WT.kaGDP
    ps[:WT_kaGTP]  = WT.kaGTP
    ps[:WT_kcat]   = WT.kcat
    ps[:WT_KM]     = WT.KM
    ps[:WT_kGDP]   = WT.kGDP
    ps[:WT_KMGDP]  = WT.KMGDP
    ps[:WT_KMGTP]  = WT.KMGTP
    # WT_kGTP is virtual
    ps[:WT_KD]     = WT.KD
    ps[:WT_kaEff]  = WT.kaEff
    # WT_kdEff is virtual

    # Mut params
    ps[:Mut_kint]  = Mut.kint
    ps[:Mut_kdGDP] = Mut.kdGDP
    ps[:Mut_kdGTP] = Mut.kdGTP
    ps[:Mut_kaGDP] = Mut.kaGDP
    # kaGTP is virtual
    ps[:Mut_kcat]  = Mut.kcat
    ps[:Mut_KM]    = Mut.KM
    ps[:Mut_kGDP]  = Mut.kGDP
    ps[:Mut_KMGDP] = Mut.KMGDP
    ps[:Mut_KMGTP] = Mut.KMGTP
    ps[:Mut_kGTP]  = Mut.kGTP     
    ps[:Mut_KD]    = Mut.KD
    ps[:Mut_kaEff] = Mut.kaEff
    # kdEff is virtual

    return ps
end

#TODO: get ode simulated and compared nicely with Eds OG results.