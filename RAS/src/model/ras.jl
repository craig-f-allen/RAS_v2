using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using SciCompDSL

@mtkmodel RAS_Base begin
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

        # Mut params
        Mut_kint::Float64
        Mut_kdGDP::Float64
        Mut_kdGTP::Float64
        Mut_kaGDP::Float64
        # kaGTP is virtual
        Mut_kcat::Float64
        Mut_KM::Float64
        Mut_kGDP::Float64
        Mut_KMGDP::Float64
        Mut_KMGTP::Float64
        Mut_kGTP::Float64     
        Mut_KD::Float64
        Mut_kaEff::Float64
        # kdEff is virtual

        # Static concentrations.
        GEF::Float64
        GAP::Float64
        GTP::Float64
        GDP::Float64

        # Protein totals
        TotalRAS::Float64
        TotalEff::Float64
        fract_mut::Float64

    end
    @variables begin

        # dyanmically generated parameters.
        WT_kGTP(t)
        WT_kdEff(t)
        
        Mut_kaGTP(t)
        Mut_kdEff(t)

        # state variables
        WT_RAS_GDP(t) => (1.0-fract_mut)*TotalRAS   # WT RAS bound to GDP
        WT_RAS_GTP(t) => 0                          #        WT RAS bounds to GTP
        WT_RAS_0(t) => 0                            # WT RAS unbound
        
        Mut_RAS_GDP(t) => fract_mut*TotalRAS        # Mut RAS bound to GDP
        Mut_RAS_GTP(t) => 0                         # Mut RAS bound to GTP
        Mut_RAS_0(t) => 0                           # Mut RAS unbound

        Eff(t) => TotalEff                          # Effector
        WT_RAS_GTP_Eff(t) => 0                      # WT RAS bound to GTP and Eff
        Mut_RAS_GTP_Eff(t) => 0                     # Mut RAS bound to GTP and Eff

        # observable variables
        RAS_GTP_Eff_Total(t)
        RAS_GTP_Total(t)

        # rate expressions
        R1(t)
        R2(t)
        R3(t)
        R4(t)
        R5(t)
        R6(t)
        R7(t)
        R8(t)
        R9(t)
        R10(t)
        R11(t)
        R12(t)
        R13(t)
        R14(t)
        
    end
    @equations begin

        # dynamically generate WT_kT and WT_kdEff
        WT_kGTP ~ WT_kGDP * WT_KMGTP * ((WT_kaGDP * WT_kdGTP) / (WT_kdGDP * WT_kaGTP)) / WT_KMGDP
        WT_kdEff ~ WT_KD * WT_kaEff

        # dynamically generate Mut_kaGTP Mut_kdEff
        Mut_kaGTP ~ Mut_kGDP * Mut_KMGTP * ((Mut_kaGDP * Mut_kdGTP) / (Mut_kdGDP * Mut_kGTP)) / Mut_KMGDP
        Mut_kdEff ~ Mut_KD * Mut_kaEff

        # rate expressions
        R1 ~ ((WT_kGDP * GEF * WT_RAS_GDP / WT_KMGDP) - (WT_kGTP * GEF * (GDP / GTP) * WT_RAS_GTP / WT_KMGTP)) / (1 + WT_RAS_GDP / WT_KMGDP + WT_RAS_GTP / WT_KMGTP + Mut_RAS_GDP / Mut_KMGDP + Mut_RAS_GTP / Mut_KMGTP)
        R2 ~ (WT_kcat * GAP * WT_RAS_GTP) / (WT_KM * (1 + Mut_RAS_GTP / Mut_KM) + WT_RAS_GTP)
        R3 ~ WT_kint * WT_RAS_GTP
        R4 ~ WT_kdGDP * WT_RAS_GDP - (WT_kaGDP * GDP) * WT_RAS_0
        R5 ~ WT_kdGTP * WT_RAS_GTP - (WT_kaGTP * GTP) * WT_RAS_0
        R6 ~ WT_kaEff * WT_RAS_GTP * Eff - WT_kdEff * WT_RAS_GTP_Eff
        R7 ~ WT_kint * WT_RAS_GTP_Eff
        R8 ~ ((Mut_kGDP * GEF * Mut_RAS_GDP / Mut_KMGDP) - (Mut_kGTP * GEF * (GDP / GTP) * Mut_RAS_GTP / Mut_KMGTP)) / (1 + WT_RAS_GDP / WT_KMGDP + WT_RAS_GTP / WT_KMGTP + Mut_RAS_GDP / Mut_KMGDP + Mut_RAS_GTP / Mut_KMGTP)
        R9 ~ (Mut_kcat * GAP * Mut_RAS_GTP) / (Mut_KM * (1 + WT_RAS_GTP / WT_KM) + Mut_RAS_GTP)
        R10 ~ Mut_kint * Mut_RAS_GTP
        R11 ~ Mut_kdGDP * Mut_RAS_GDP - (Mut_kaGDP * GDP) * Mut_RAS_0
        R12 ~ Mut_kdGTP * Mut_RAS_GTP - (Mut_kaGTP * GTP) * Mut_RAS_0
        R13 ~ Mut_kaEff * Mut_RAS_GTP * Eff - Mut_kdEff * Mut_RAS_GTP_Eff
        R14 ~ Mut_kint * Mut_RAS_GTP_Eff

        # rate equations
        D(WT_RAS_GDP) ~ -R1+R2+R3-R4+R7
        D(WT_RAS_GTP) ~ R1-R2-R3-R5-R6
        D(WT_RAS_0) ~ R4+R5
        D(Eff) ~ -R6+R7-R13+R14
        D(WT_RAS_GTP_Eff) ~ R6-R7
        D(Mut_RAS_GDP) ~ -R8+R9+R10-R11+R14
        D(Mut_RAS_GTP) ~ R8-R9-R10-R12-R13
        D(Mut_RAS_0) ~ R11+R12
        D(Mut_RAS_GTP_Eff) ~ R13-R14

        # observable equations
        RAS_GTP_Eff_Total ~ WT_RAS_GTP_Eff + Mut_RAS_GTP_Eff
        RAS_GTP_Total ~ WT_RAS_GTP_Eff + Mut_RAS_GTP_Eff + WT_RAS_GTP + Mut_RAS_GTP
    end
end

