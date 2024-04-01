using MAT
project_dir = "/Users/tyler/projects/waiting"
model_code_dir = joinpath(project_dir, "matlab_code/accumulation-model/")

include(joinpath(model_code_dir, "analytical_model.jl"))


overwrite   = true
res_dir     = joinpath(project_dir, "results")
data_dir    = joinpath(project_dir, "data")
ratlistfile = joinpath(data_dir, "included_rats.mat")
rats        = matread(ratlistfile)
rats        = rats["included_rats"]
println(rats)
hess = analytical_model.compute_hessian(rats,res_dir = res_dir, overwrite=overwrite)
