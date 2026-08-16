function build_params_for_base_ras(WT::RAS.WTKineticParams, Mut::RAS.MutantKineticParams, GAP::Real, GEF::Real, GDP::Real, GTP::Real, TotalRAS::Real, TotalEff::Real, fract_mut::Real)

    # T covers whatever real type the caller passed in for the varying args (Float64 for
    # a normal solve, a ForwardDiff.Dual when this is called from inside a gradient) - the
    # WT/Mut struct fields are always Float64 and convert into T for free.
    T = promote_type(Float64, typeof(GAP), typeof(GEF), typeof(GDP), typeof(GTP), typeof(TotalRAS), typeof(TotalEff), typeof(fract_mut))
    return Dict{Symbol, T}(
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

function build_params_for_ras_tricomplex(WT::RAS.WTKineticParams, WT_Drug::TriDrugKineticParams, Mut::RAS.MutantKineticParams, Mut_Drug::TriDrugKineticParams, GAP::Real, GEF::Real, GDP::Real, GTP::Real, TotalRAS::Real, TotalEff::Real, fract_mut::Real, Drug::Real, CYPA::Real)

    T = promote_type(Float64, typeof(GAP), typeof(GEF), typeof(GDP), typeof(GTP), typeof(TotalRAS), typeof(TotalEff), typeof(fract_mut), typeof(Drug), typeof(CYPA))
    return Dict{Symbol, T}(

        # Abundances (ICs)
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
#
# GAP/GEF/.../CYPA are `Real` (not `Float64`) so a caller can hand in a ForwardDiff.Dual for
# one of them and get a differentiable parameter set back - that's what makes DoseResponseProblem
# AD-ready with no separate code path. `overrides` covers every *other* named parameter (kinetic
# rate constants, K_D's, etc.) the same way, without needing a dedicated kwarg per parameter.
function build_tricomplex_param_pairs(mutant::Symbol, fract_mut::Real;
        GAP::Real=6e-11, GEF::Real=2e-10, GDP::Real=18e-6, GTP::Real=180e-6,
        TotalRAS::Real=4e-7, TotalEff::Real=4e-7, Drug0::Real=1e2, CYPA::Real=1e-6,
        overrides::AbstractDict{Symbol}=Dict{Symbol,Float64}(),
        xlsx_path::String=DEFAULT_KINETIC_PARAMS_PATH)
    mutants_base     = get_mutant_params(xlsx_path)
    mutants_tri_drug = get_mutant_tri_drug_params(xlsx_path)

    pd_base = build_params_for_base_ras(WT, mutants_base[mutant], GAP, GEF, GDP, GTP, TotalRAS, TotalEff, fract_mut)
    pd_tri  = build_params_for_ras_tricomplex(WT, mutants_tri_drug[:WT], mutants_base[mutant], mutants_tri_drug[mutant],
                                               GAP, GEF, GDP, GTP, TotalRAS, TotalEff, fract_mut, Drug0, CYPA)

    return Pair[k => get(overrides, k, v) for (k, v) in merge(pd_base, pd_tri)]
end