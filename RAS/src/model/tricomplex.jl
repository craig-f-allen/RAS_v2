@mtkmodel RAS_Tricomplex begin
    @parameters begin

        # WT params
        WT_kint
        WT_kdGDP
        WT_kdGTP
        WT_kaGDP
        WT_kaGTP   
        WT_kcat
        WT_KM
        WT_kGDP
        WT_KMGDP
        WT_KMGTP
        # WT_kGTP is virtual
        WT_KD
        WT_kaEff
        # WT_kdEff is virtual

        # Mut params
        Mut_kint
        Mut_kdGDP
        Mut_kdGTP
        Mut_kaGDP
        # kaGTP is virtual
        Mut_kcat
        Mut_KM
        Mut_kGDP
        Mut_KMGDP
        Mut_KMGTP
        Mut_kGTP     
        Mut_KD
        Mut_kaEff
        # kdEff is virtual

        # Static concentrations.
        GEF
        GAP
        GTP
        GDP

        # Protein totals
        TotalRAS
        TotalEff
        fract_mut

        # Drugging parameters
        Drug_0
        CYPA_0
        k_on_1 = 1e7
        K_D_1 = 195e-9
        WT_K_D_2
        Mut_K_D_2
        WT_k_on_tricomplex = 1e7
        Mut_k_on_tricomplex = 1e7
        

    end
    @variables begin

        # dynamically generated parameters.
        WT_kGTP(t)
        WT_kdEff(t)
        
        Mut_kaGTP(t)
        Mut_kdEff(t)

        k_off_1(t)
        WT_k_off_tricomplex(t)
        Mut_k_off_tricomplex(t)

        # state variables
        WT_RAS_GDP(t) = (1.0-fract_mut)*TotalRAS   # WT RAS bound to GDP
        WT_RAS_GTP(t) = 0                          # WT RAS bounds to GTP
        WT_RAS_0(t) =  0                            # WT RAS unbound
        
        Mut_RAS_GDP(t) = fract_mut*TotalRAS        # Mut RAS bound to GDP
        Mut_RAS_GTP(t) = 0                         # Mut RAS bound to GTP
        Mut_RAS_0(t) = 0                           # Mut RAS unbound

        Eff(t) = TotalEff                          # Effector
        WT_RAS_GTP_Eff(t) = 0                      # WT RAS bound to GTP and Eff
        Mut_RAS_GTP_Eff(t) = 0                     # Mut RAS bound to GTP and Eff

        Drug(t) = Drug_0
        CYPA(t) = CYPA_0
        Bicomplex(t) = 0
        WT_Tricomplex(t) = 0
        Mut_Tricomplex(t) = 0

        # observable variables
        RAS_GTP_Eff_Total(t)
        RAS_GTP_Total(t)
        Tricomplex_Total(t)

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
        R15(t)
        R16(t)
        R17(t)
        
    end
    @equations begin

        # dynamically generate WT_kT and WT_kdEff
        WT_kGTP ~ WT_kGDP * WT_KMGTP * ((WT_kaGDP * WT_kdGTP) / (WT_kdGDP * WT_kaGTP)) / WT_KMGDP
        WT_kdEff ~ WT_KD * WT_kaEff

        # dynamically generate Mut_kaGTP Mut_kdEff
        Mut_kaGTP ~ Mut_kGDP * Mut_KMGTP * ((Mut_kaGDP * Mut_kdGTP) / (Mut_kdGDP * Mut_kGTP)) / Mut_KMGDP
        Mut_kdEff ~ Mut_KD * Mut_kaEff

        # dynamically generate k_offs
        k_off_1 ~ k_on_1*K_D_1
        WT_k_off_tricomplex ~ WT_k_on_tricomplex*WT_K_D_2
        Mut_k_off_tricomplex ~ Mut_k_on_tricomplex*Mut_K_D_2

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
        R15 ~ k_on_1*Drug*CYPA - k_off_1*Bicomplex
        R16 ~ WT_k_on_tricomplex*Bicomplex*WT_RAS_GTP - WT_k_off_tricomplex*WT_Tricomplex
        R17 ~ Mut_k_on_tricomplex*Bicomplex*Mut_RAS_GTP - Mut_k_off_tricomplex*Mut_Tricomplex

        # rate equations
        D(WT_RAS_GDP) ~ -R1+R2+R3-R4+R7
        D(WT_RAS_GTP) ~ R1-R2-R3-R5-R6-R16
        D(WT_RAS_0) ~ R4+R5
        D(Eff) ~ -R6+R7-R13+R14
        D(WT_RAS_GTP_Eff) ~ R6-R7
        D(Mut_RAS_GDP) ~ -R8+R9+R10-R11+R14
        D(Mut_RAS_GTP) ~ R8-R9-R10-R12-R13-R17
        D(Mut_RAS_0) ~ R11+R12
        D(Mut_RAS_GTP_Eff) ~ R13-R14
        D(Drug) ~ -R15
        D(CYPA) ~ -R15
        D(Bicomplex) ~ R15-R16-R17
        D(WT_Tricomplex) ~ R16
        D(Mut_Tricomplex) ~ R17

        # observable equations
        RAS_GTP_Eff_Total ~ WT_RAS_GTP_Eff + Mut_RAS_GTP_Eff
        RAS_GTP_Total ~ WT_RAS_GTP_Eff + Mut_RAS_GTP_Eff + WT_RAS_GTP + Mut_RAS_GTP
        Tricomplex_Total ~ WT_Tricomplex + Mut_Tricomplex
    end
end