using multiwfs

r0 = 0.1031
r0_ncp = 1.0
vk_atm = VonKarman(0.3 * 10.0 / 3.0, 0.25 * r0^(-5/3))
vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))

# fast_controller, slow_controller, R, vk_atm, vk_ncp, f_noise_crossover
f_loop = 1000.0
R = 10
fast_controller = FilteredIntegrator(0.4, 0.999, ar1_filter(15.0, f_loop, "high"), 1/f_loop)
slow_controller = FilteredIntegrator(1.4, 0.999, ZPKFilter(0, 0, 1), R/f_loop)
sim = Simulation(f_loop, fast_controller, slow_controller, 10, vk_atm, vk_ncp, 500.0)

notched_error_X(sim)
notched_error_Y(sim)

gains_slow = 0.0:0.1:2.0
gains_fast = 0.0:0.1:1.0
cutoff_freqs = 0.0:2.0:40.0

opt_xerrs, opt_yerrs, opt_controlpars = [], [], []
r0_ncp_vals = [0.6, 0.8, 1.0, 1.4, 1.6, 2.0, 4.0, 6.0]
for r0_ncp in r0_ncp_vals
    vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))
    X_errors = zeros(length(gains_slow), length(gains_fast), length(cutoff_freqs))
    Y_errors = zeros(length(gains_slow), length(gains_fast), length(cutoff_freqs))

    @showprogress for i in eachindex(gains_slow)
        for j in eachindex(gains_fast)
            @threads for k in eachindex(cutoff_freqs)
                X_errors[i,j,k], Y_errors[i,j,k] = notched_errors(gains_slow[i], gains_fast[j], cutoff_freqs[k], vk_ncp)
            end
        end
    end

    idx_xmin = argmin(X_errors)
    res = optimize(x -> notched_errors(x..., vk_ncp)[1], [gains_slow[idx_xmin[1]], gains_fast[idx_xmin[2]], cutoff_freqs[idx_xmin[3]]], NelderMead(), Optim.Options(); autodiff=:forward)

    xerr, yerr = notched_errors(res.minimizer..., vk_ncp)

    push!(opt_xerrs, xerr)
    push!(opt_yerrs, yerr)
    push!(opt_controlpars, res.minimizer)
end

opt_controlpars = hcat(opt_controlpars...)

r0_ncp_labels = [0.6, 1.0, 2.0, 3.0, 6.0]
p1 = scatter(r0_ncp_vals, opt_xerrs, xscale=:log10, xticks=(r0_ncp_labels, r0_ncp_labels), label=nothing, xlims=(0.5, 7.0), ylims=(0.5, 1.0), xlabel="NCP r₀ (m)", ylabel="CL error at X (rad)")
plot!(r0_ncp_vals, opt_xerrs, label=nothing)

begin
    pl = [p1]
    for (i, l) in enumerate(["Slow gain", "Fast gain", "HPF cutoff freq (Hz)"])
        p = plot(r0_ncp_vals, opt_controlpars[i,:], xlabel="NCP r₀ (m)", ylabel=l, xlims=(0.5, 7.0), xticks=(r0_ncp_labels, r0_ncp_labels), label=nothing, xscale=:log10)
        push!(pl, p)
    end
    plot(pl...)
end

"""begin
    pl = []
    for i in eachindex(res.minimizer)
        minres = copy(res.minimizer)
        par_sweep = (0.95:0.005:1.05) * minres[i]
        errs = []
        for p in par_sweep
            minres[i] = p
            push!(errs, notched_errors(minres...)[1])
        end
        push!(pl, plot(par_sweep, errs))
    end
    plot(pl...)
end"""

