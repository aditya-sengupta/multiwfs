using multiwfs
using multiwfs: Hcont, Hfilter, is_stable
using Plots
using Base: product
using Base.Meta: parse
using Base.Threads: @threads
using ProgressMeter: @showprogress
using Optim
using ForwardDiff: gradient, jacobian

f_loop = 1000.0
fr = 10 .^ (-2:0.01:log10(f_loop/2))
sr = 2π .* im .* fr
sT = sr / f_loop
no_filter = ZPKFilter(0, 0, 1)
r0 = 0.1031
vk_atm = VonKarman(0.3 * 10.0 / 3.0, 0.25 * r0^(-5/3))
f_noise_crossover = 500.0
noise_normalization = psd_von_karman(f_noise_crossover, vk_atm)

function systems(gain_slow, gain_fast, f_cutoff)
    ar1_high = ar1_filter(f_cutoff, f_loop, "high")
    # f_loop, frame_delay, gain, leak, fpf
    sys_fast = AOSystem(f_loop, 1.0, gain_fast, 0.999, 1, ar1_high)
    sys_slow = AOSystem(f_loop / 10, 0.1, gain_slow, 0.999, 1, no_filter)
    return sys_fast, sys_slow
end

function notched_errors(gain_slow, gain_fast, f_cutoff, vk_ncp)
    if gain_slow > 2.0
        return Inf, Inf
    end
    sys_fast, sys_slow = systems(gain_slow, gain_fast, f_cutoff)
    Cfast = sT -> Hcont(sys_fast, sT * f_loop) * Hfilter(sys_fast, sT * f_loop)
    Cslow = sT -> Hcont(sys_slow, sT * f_loop)
    margin_values = margins(Vector{AOSystem}([sys_fast, sys_slow]))
    if margin_values.gm < 2.5 || margin_values.pm < 45
        return Inf, Inf
    end
    err_X, err_Y = notched_error_X(Cfast, Cslow, 10, vk_atm, vk_ncp, f_noise_crossover, f_min=1.0), notched_error_Y(Cfast, Cslow, 10, vk_atm, vk_ncp, f_noise_crossover, f_min=1.0)
    return err_X, err_Y
end

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

