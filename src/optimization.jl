using Base.Threads: @threads
using Base: product
using FLoops

"""
Grid search over simulation parameters

sim_generator - callable, k parameters
Generates a Simulation object for each set of parameters.

parameter_ranges - vector of k iterables
The parameters to try out; grid search will be over the Cartesian product of these.
"""
function grid_search(sim_generator, parameter_ranges)
    @floop ThreadedEx() for pars in product(parameter_ranges...)
        this_error = notched_error_X(sim_generator(pars...))
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
    for pars in product(parameter_ranges...)
        this_error = notched_error_X(sim_generator(pars...))
        if isless(this_error, best_error)
            best_error = this_error
            optpars = pars
        end
    end
    
    optpars
end

export grid_search, grid_search_serial