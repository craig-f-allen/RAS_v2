function build_params_for_base_ras(WT::RAS.WTKineticParams, Mut::RAS.MutantKineticParams, GAP::Float64, GEF::Float64, GDP::Float64, GTP::Float64, TotalRAS::Float64, TotalEff::Float64, fract_mut::Float64)

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

function build_params_for_ras_tricomplex(WT::RAS.WTKineticParams, WT_Drug::TriDrugKineticParams, Mut::RAS.MutantKineticParams, Mut_Drug::TriDrugKineticParams, GAP::Float64, GEF::Float64, GDP::Float64, GTP::Float64, TotalRAS::Float64, TotalEff::Float64, fract_mut::Float64, Drug::Float64, CYPA::Float64)

    return Dict{Symbol, Float64}(

        # Abudences (ICs)
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
        :Mut_kaEff => Mut.kaEff,

        # Drug params
        :Drug_0 => Drug,
        :CYPA_0 => CYPA,
        :K_D_1 => 55.3e-9, # M, from PMID: 38593348
        :WT_K_D_2   => WT_Drug.K_D_2,
        :Mut_K_D_2  => Mut_Drug.K_D_2

    )
end

# Loads mutant kinetics for `mutant` from `xlsx_path`, builds and merges the base + tricomplex
# param dicts, and returns them as a concrete Pair vector ready for ODEProblem(sys, param_pairs, tspan).
# Shared by DoseResponseProblem/SingleSimProblem so both constructors don't duplicate this.
function build_tricomplex_param_pairs(mutant::Symbol, fract_mut::Real;
        GAP::Float64=6e-11, GEF::Float64=2e-10, GDP::Float64=18e-6, GTP::Float64=180e-6,
        TotalRAS::Float64=4e-7, TotalEff::Float64=4e-7, Drug0::Float64=1e2, CYPA::Float64=1e-6,
        xlsx_path::String=DEFAULT_KINETIC_PARAMS_PATH)
    fm = Float64(fract_mut)
    mutants_base     = get_mutant_params(xlsx_path)
    mutants_tri_drug = get_mutant_tri_drug_params(xlsx_path)

    pd_base = build_params_for_base_ras(WT, mutants_base[mutant], GAP, GEF, GDP, GTP, TotalRAS, TotalEff, fm)
    pd_tri  = build_params_for_ras_tricomplex(WT, mutants_tri_drug[:WT], mutants_base[mutant], mutants_tri_drug[mutant],
                                               GAP, GEF, GDP, GTP, TotalRAS, TotalEff, fm, Drug0, CYPA)

    return Pair[k => v for (k, v) in merge(pd_base, pd_tri)]
end