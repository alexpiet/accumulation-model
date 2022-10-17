model_code_dir = "/Users/oroville/projects/pbups_dyn/code/accumulation-model/"

include(joinpath(model_code_dir, "analytical_model.jl"))


overwrite   = true
project_dir = "/Users/oroville/projects/waiting"
res_dir     = joinpath(project_dir, "results")
data_dir    = joinpath(project_dir, "data")
ratlistfile = joinpath(data_dir, "included_rats.mat")
rats        = matread(ratlistfile)
rats        = rats["included_rats"]
println(rats)
hess = analytical_model.compute_hessian(rats,res_dir = res_dir, overwrite=overwrite)
