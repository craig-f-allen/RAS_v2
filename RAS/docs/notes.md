TODO
#

- IF CHANGING PARAMS THAT ARE USED FOR ICS: MUST USE u0 OBJECT AS WELL (can make function to construct it from params)

Learnings
# 
- Need ROsenbrock23() as solver as it is a stiff system
- If using DynamicSS need higher tolerances:sol = solve(prob, DynamicSS(Rosenbrock23()); abstol=1-12, reltol=1e-10)