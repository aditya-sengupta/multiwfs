using Pkg
Pkg.activate(".") # this'll only be run from top-level

include("openloop_psds.jl")
include("single_ic_fulleval.jl")
include("double_ic_fulleval.jl")
include("double_ichpf_fulleval.jl")
include("param_heatmaps.jl")