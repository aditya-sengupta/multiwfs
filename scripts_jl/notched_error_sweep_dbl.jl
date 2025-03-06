using multiwfs
using Plots
using Base: product
using Base.Threads: @threads
using ProgressMeter: @showprogress
using Optim
using ForwardDiff

# min |X / phi|^2 |PSD atm| + (|X / Lfast|^2 + |X / Lslow|^2) |PSD NCP| 
# subject to |Y / phi|^2 |PSD atm| + (|Y / Lfast|^2 + |Y / Lslow|^2) |PSD NCP| < threshold

f_loop = 1000.0
fr = 10 .^ (-2:0.01:log10(f_loop/2))
sr = 2π .* im .* fr
sT = sr / f_loop
no_filter = ZPKFilter(0, 0, 1)
vk_atm = VonKarman(0.3 * 10.0 / 3.0, 0.25 * (0.1031)^(-5/3))
vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25)
f_noise_crossover = 500.0
noise_normalization = psd_von_karman(f_noise_crossover, vk_atm)

function systems(α, β, gain_fast, f_cutoff)
    ar1_high = ar1_filter(f_cutoff, f_loop, "high")
    # f_loop, frame_delay, gain, leak, fpf
    sys_fast = AOSystem(f_loop, 1.0, gain_fast, 0.999, 1, ar1_high)
    sys_slow = make_double_integrator_system(α, β)
    return sys_fast, sys_slow
end

function notched_errors(α, β, gain_fast, f_cutoff)
    if α < 0 || gain_fast < 0 || f_cutoff < 0
        return Inf, Inf
    end
    sys_fast, sys_slow = systems(α, β, gain_fast, f_cutoff)
    Cfast = sT -> Hcont(sys_fast, sT * f_loop) * Hfilter(sys_fast, sT * f_loop)
    Cslow = sT -> Hcont(sys_slow, sT * f_loop) * Hfilter(sys_slow, sT * f_loop)
    if !is_stable((sys_fast, sys_slow))
       return Inf, Inf
    end
    err_X, err_Y = notched_error_X(Cfast, Cslow, 10, vk_atm, vk_ncp, f_noise_crossover, f_min=1.0), notched_error_Y(Cfast, Cslow, 10, vk_atm, vk_ncp, f_noise_crossover, f_min=1.0)
    return err_X, err_Y
end

notched_errors(0.995, 0.0, 0.4, 15.0)

nyquist_plot(systems(0.995, -0.1, 0.4, 15.0))

s = systems(0.995, -0.1, 0.4, 15.0)
ten_etf_plots(controller(s[1]), controller(s[2]), 10, vk_atm, vk_ncp, noise_normalization)

alphas = 0.99:0.001:0.999
betas = -0.1:0.01:0.01

X_errors = zeros(length(alphas), length(betas))
Y_errors = zeros(length(alphas), length(betas))

@showprogress for i in eachindex(alphas)
    @threads for m in eachindex(betas)
        X_errors[i,m], Y_errors[i,m] = notched_errors(alphas[i], betas[m], 0.52, 32.155)
    end
end

minimum(X_errors)
idx_xmin = argmin(X_errors)
alphas[idx_xmin[1]], betas[idx_xmin[2]]
Y_errors[argmin(X_errors)]

X_errors[Y_errors .>= 1.0] .= Inf
argmin(X_errors)
minimum(X_errors)

notched_errors(0.9999, -0.1, 0.52, 32.155)

pars = [0.995, -0.08, 0.4, 10.0]
sys_fast, sys_slow = systems(pars...)
Cfast = sT -> Hcont(sys_fast, sT * f_loop)(sys_fast, sT * f_loop)
Cslow = sT -> Hcont(sys_slow, sT * f_loop)
plot_integrands(Cfast, Cslow, 10, vk_atm, vk_ncp, noise_normalization)

begin
    plots_contrib = []
    pars = [0.995, -0.08, 0.4, 10.0]
    sys_fast, sys_slow = systems(pars...)
    Cfast = sT -> Hcont(sys_fast, sT * f_loop) * Hfilter(sys_fast, sT * f_loop)
    Cslow = sT -> Hcont(sys_slow, sT * f_loop) * Hfilter(sys_slow, sT * f_loop)
    errX, errY = notched_errors(pars...)
    p_v = plot(legend=:bottomleft, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="Closed-loop residual (rad/Hz)", title="X error = $(round(errX, digits=3)) rad, Y error = $(round(errY, digits=3)) rad", ylims=(1e-10, 1e2))
    for errsource in ["atm_error_at_f_X", "ncp_error_at_f_X", "ncp_error_at_f_Y", "noise_error_at_f_X"]
        err_source_fn = eval(parse(errsource))
        plot!(fr, err_source_fn.(fr, Cfast, Cslow, 10, Ref(vk_atm), Ref(vk_ncp), noise_normalization, f_loop), label=errsource)
    end
    p_v
end

begin
    plots_etf = []
    for v in ["X", "Y"]
        ne = eval(parse("notched_error_$v"))
        p_v = plot(legend=:bottomright, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="ETF", ylims=(1e-10, 1e2), title="$v error = $(round(ne(Cfast, Cslow, 10, vk_atm, vk_ncp, f_noise_crossover), digits=3)) rad")
        for (fname, c, s) in zip(["phi_to_$v", "Lfast_to_$v", "Lslow_to_$v", "Nfast_to_$v", "Nslow_to_$v"], [1, 2, 2, 3, 3], [:solid, :solid, :dash, :solid, :dash])
            f = eval(parse(fname))
            plot!(fr, abs2.(f.(sT, Cfast, Cslow, 10)), label="|$fname|²", c=c, ls=s)
        end
        push!(plots_etf, p_v)
    end
    plot(plots_etf...)
end