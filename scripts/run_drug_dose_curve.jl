using RAS, Plots, Printf

prob = DoseResponseProblem(:G12V, 0.25)
ic50 = drug_dose_response_ic50!(prob)

@printf("IC50 = %.4g M\n", ic50)
plot(prob.doses, prob.ys, xscale=:log10, xlabel="Drug Dose (M)", ylabel="RAS_GTP_Eff_Total", title="Drug Dose Response Curve", legend=false)
vline!([ic50], label=@sprintf("IC50 = %.3g M", ic50), legend=true)
