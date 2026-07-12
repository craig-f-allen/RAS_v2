using RAS

file_path = normpath(joinpath(@__DIR__, "../RAS/data/kinetic_parameters.xlsx"))

mutants = get_mutant_params(file_path)

