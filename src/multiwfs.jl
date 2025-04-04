module multiwfs
    using ProgressMeter
    using Base.Threads
    using Plots
    const PROJECT_ROOT = pkgdir(multiwfs)

    include("von_karman.jl")
    include("zpkfilter.jl")
    include("controller.jl")
    include("simulation.jl")
    include("stability.jl")
    include("performance.jl")
    include("optimization.jl")
    include("plotting.jl")
    include("utils.jl")
end 
