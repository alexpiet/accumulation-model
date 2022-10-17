"""
analytical_model is a package for fitting an analytical, unbounded version of Bing's model

The model has 8 parameters used to predict choices
"""
module analytical_model

using MAT
using ForwardDiff
using SpecialFunctions
using Glob

export compute_hessian

include("compute_LL.jl")
include("compute_hessian.jl")
end
