using RAS

prob = SingleSimProblem(:G12V, 0.25)
u_ss = run_ss!(prob, 1e5) # steady state at drug dose = 1e5
