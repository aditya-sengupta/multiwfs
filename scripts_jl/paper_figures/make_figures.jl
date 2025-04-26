using Pkg
Pkg.activate(".") # this'll only be run from top-level

include("openloop_psds.jl")
include("single_ic_fulleval.jl")
