TODO
#

1) Make efficient IC50 finding function.
2) Find and lot IC50 of each mutant over the range of K_D_2s
3) Fit other mutants K_D_2s
4) Maximize RAS_GTP_Eff response and see changes.
Learnings
# 
- Need ROsenbrock23() as solver as it is a stiff system
- If using DynamicSS need higher tolerances:sol = solve(prob, DynamicSS(Rosenbrock23()); abstol=1-12, reltol=1e-10)