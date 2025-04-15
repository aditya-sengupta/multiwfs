using Base.Threads: @threads
using Base: product
using FLoops
using ProgressMeter

"""
Grid search over simulation parameters

sim_generator - callable, k parameters
Generates a Simulation object for each set of parameters.

parameter_ranges - vector of k iterables
The parameters to try out; grid search will be over the Cartesian product of these.
"""
function grid_search(sim_generator, parameter_ranges)
    @floop ThreadedEx() for pars in product(parameter_ranges...)
        this_error = Inf
        try
            this_error = notched_error_X(sim_generator(pars...))
        catch e
        end
        @reduce() do (best_error; this_error), (optpars; pars)
            if isless(this_error, best_error)
                best_error = this_error
                optpars = pars
            end
        end
    end
    
    optpars
end

function grid_search_serial(sim_generator, parameter_ranges)
    best_error = Inf
    optpars = zeros(length(parameter_ranges))
    @showprogress for pars in product(parameter_ranges...)
        this_error = Inf
        try
            this_error = notched_error_X(sim_generator(pars...))
        catch e
        end
        if isless(this_error, best_error)
            best_error = this_error
            optpars = pars
        end
    end
    
    optpars
end

function grid_search_coarse_to_fine(sim_generator, parameter_limits; npoints_per_parameter=11, search="serial", niter=3)
    parameters_this_iteration = nothing
    parameter_ranges = nothing
    gs = nothing
    if search == "serial"
        gs = grid_search_serial
    else
        gs = grid_search
    end
    for _ in 1:niter
        parameter_steps = [(x[2] - x[1]) / (npoints_per_parameter - 1) for x in parameter_limits]
        parameter_ranges = [x[1]:st:x[2] for (x, st) in zip(parameter_limits, parameter_steps)]
        parameters_this_iteration = gs(sim_generator, parameter_ranges)
        parameter_limits = [[max(p-st,pl[1]), min(p+st,pl[2])] for (p, st, pl) in zip(parameters_this_iteration, parameter_steps, parameter_limits)]
    end
    return parameters_this_iteration, parameter_ranges
end

function parameter_resolution(sim_generator, parameter_limits; kwargs...)
    parameters_central, parameter_ranges = grid_search_coarse_to_fine(sim_generator, parameter_limits; kwargs...)
    
end

export grid_search, grid_search_serial, grid_search_coarse_to_fine