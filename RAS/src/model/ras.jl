#TODO: MTK model, use algebraic equations for dynamic params so that you can use strucural anaysis, AD for optimizaiton or Turing, etc. parameter loading was still necessarily - ensured that is done properly. May need to consider hhow to do this properly when doing mutant mutant vs WT mutant vs WT WT

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@mtkmodel RAS begin
    @parameters begin

        # WT params
        WT_kint::Float64
        WT_kdGDP::Float64
        WT_kdGTP::Float64
        WT_kaGDP::Float64
        WT_kaGTP::Float64   
        WT_kcat::Float64
        WT_KM::Float64
        WT_kGDP::Float64
        WT_KMGDP::Float64
        WT_KMGTP::Float64
        # WT_kGTP is virtual
        WT_KD::Float64
        WT_kaEff::Float64
        # WT_kdEff is virtual

        # Mutant params
        kint::Float64
        kdGDP::Float64
        kdGTP::Float64
        kaGDP::Float64
        # kaGTP is virtual
        kcat::Float64
        KM::Float64
        kGDP::Float64
        KMGDP::Float64
        KMGTP::Float64
        kGTP::Float64     
        KD::Float64
        kaEff::Float64
        # kdEff is virtual

    end
    @variables begin
        
        
    end
    @equations begin
        x ~ k * D(x)
        D(y) ~ -k
        [D(z[i]) ~ 1]...
    end
end