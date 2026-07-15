function build_params_for_mtk(WT::RAS.WTKineticParams, Mut::RAS.MutantKineticParams, GAP::Float64, GEF::Float64, GDP::Float64, GTP::Float64, TotalRAS::Float64, TotalEff::Float64, fract_mut::Float64)

    return Dict{Symbol, Float64}(
        :TotalRAS  => TotalRAS,
        :TotalEff  => TotalEff,
        :fract_mut => fract_mut,
        :GAP       => GAP,
        :GEF       => GEF,
        :GDP       => GDP,
        :GTP       => GTP,
        
        # WT params
        :WT_kint   => WT.kint,
        :WT_kdGDP  => WT.kdGDP,
        :WT_kdGTP  => WT.kdGTP,
        :WT_kaGDP  => WT.kaGDP,
        :WT_kaGTP  => WT.kaGTP,
        :WT_kcat   => WT.kcat,
        :WT_KM     => WT.KM,
        :WT_kGDP   => WT.kGDP,
        :WT_KMGDP  => WT.KMGDP,
        :WT_KMGTP  => WT.KMGTP,
        :WT_KD     => WT.KD,
        :WT_kaEff  => WT.kaEff,

        # Mut params
        :Mut_kint  => Mut.kint,
        :Mut_kdGDP => Mut.kdGDP,
        :Mut_kdGTP => Mut.kdGTP,
        :Mut_kaGDP => Mut.kaGDP,
        :Mut_kcat  => Mut.kcat,
        :Mut_KM    => Mut.KM,
        :Mut_kGDP  => Mut.kGDP,
        :Mut_KMGDP => Mut.KMGDP,
        :Mut_KMGTP => Mut.KMGTP,
        :Mut_kGTP  => Mut.kGTP,
        :Mut_KD    => Mut.KD,
        :Mut_kaEff => Mut.kaEff
    )
end

#TODO: get ode simulated and compared nicely with Eds OG results.