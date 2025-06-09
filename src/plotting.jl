using Plots
using DSP
using Base.Meta: parse

symbol_lookup = Dict(
    "phi" => "ϕ",
    "Lfast" => "Lfast",
    "Lslow" => "Lslow",
    "Nfast" => "Nfast",
    "Nslow" => "Nslow",
    "atm" => "Atmosphere",
    "ncp" => "NCP",
    "noise" => "Noise"
)

function fn_to_equation(fname)
    m = match(r"(.+)_to_(.)", fname)
    return m[2] * "/" * symbol_lookup[m[1]]
end

function errsource_to_equation(fname)
    m = match(r"(.+)_error_at_f_(.)", fname)
    return symbol_lookup[m[1]] * " at " * m[2]
end

function nyquist_plot(sim; mark_gm_pm=true, label="Nyquist plot", kwargs...)
    success = :green
    nyquist_contour, gm, gm_point, pm, pm_point = nyquist_and_margins(sim)
    plot(real(nyquist_contour), imag(nyquist_contour), xlim=(-1.1,1.1), ylim=(-1.1,1.1), aspect_ratio=:equal, label=label; legend=:topright, kwargs...)
    phasegrid = range(-π, π, length=500)
    xunit, yunit = cos.(phasegrid), sin.(phasegrid)
    if !is_stable(gm, pm)
        success = :red
    end
    if mark_gm_pm
        vline!([-1/2.5], ls=:dash, color=:grey, label=nothing)
        scatter!([real(gm_point)], [imag(gm_point)], label="GM = $(round(gm, digits=2))", color=:grey)
        plot!([-2,0,-2], [-2,0,2], ls=:dash, color=4, label=nothing)
        if !isnothing(pm_point)
            scatter!([real(pm_point)], [imag(pm_point)], label="PM = $(round(pm, digits=2))", color=4)
        end
    end
    plot!(xunit, yunit, ls=:dash, label=nothing, color=success)
end

function ten_etf_plots(sim, cs=[1,2,2,3,3]; kwargs...)
    plots_etf = []
    for v in ["X", "Y"]
        ne = eval(parse("notched_error_$v"))
        p_v = plot(;legend=:bottomright, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="|Error transfer function|²", ylims=(1e-10, 1e2), title="$v error = $(round(ne(sim), digits=3)) rad", xticks=exp10.(-3:2), kwargs...)
        for (fname, c, s) in zip(["phi_to_$v", "Lfast_to_$v", "Lslow_to_$v", "Nfast_to_$v", "Nslow_to_$v"], cs, [:solid, :dash, :solid, :dash, :solid])
            f = eval(parse(fname))
            label = fn_to_equation(fname)
            plot!(sim.fr, abs2.(f.(sim.sT, Ref(sim))) .+ 1e-30, label="|$label|²", c=c, ls=s; kwargs...)
        end
        push!(plots_etf, p_v)
    end
    plots_etf
end

function five_psd_plots(sim; cs=[1,2,3,4,5], kwargs...)
    ol_atm = psd_von_karman.(sim.fr, Ref(sim.vk_atm))
    plot(sim.fr, ol_atm, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="Power (rad²/Hz)", label="Open-loop atm", alpha=0.5, c=cs[1];xticks=exp10.(-3:2), legend=:topright, kwargs...)
    plot!(sim.fr, psd_von_karman.(sim.fr, Ref(sim.vk_ncp)), label="Open-loop NCP", alpha=0.5, c=cs[2]; kwargs...)
    hline!([psd_von_karman(sim.f_noise_crossover, sim.vk_atm)], label="Open-loop noise", alpha=0.5, c=cs[3]; kwargs...)
    plot!(sim.fr, error_at_f_X.(sim.fr, Ref(sim)), label="Closed loop at X", c=cs[4]; kwargs...)
    plot!(sim.fr, error_at_f_Y.(sim.fr, Ref(sim)), label="Closed loop at Y", c=cs[5]; kwargs...)
end

function plot_integrand(errsource, sim; kwargs...)
    p = plot()
    plot_integrand!(errsource, sim; kwargs...)
    p
end


function plot_integrand!(errsource, sim; kwargs...)
    label = errsource_to_equation(errsource)
    err_source_fn = eval(parse(errsource))
    plot!(sim.fr, err_source_fn.(sim.fr, Ref(sim)) .+ 1e-20, label=label, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="Closed-loop residual (rad²/Hz)"; ylims=(1e-5, 1e3), legend=:topright, yticks=exp10.(-5:2:5), xticks=exp10.(-3:2), kwargs...)
    # plot!(sim.fr, err_source_fn.(sim.fr, Ref(sim)) .+ 1e-20, label=label, xlabel="Frequency (Hz)", ylabel="Closed-loop residual (rad²/Hz)"; legend=:topright, kwargs...)
end

function plot_integrands(v, sim, cs; kwargs...)
    if v == "XY"
        vt = "X"
    else
        vt = v
    end
    plot_integrand("atm_error_at_f_$vt", sim; legend=:topright, c=cs[1], kwargs...)
    plot_integrand!("ncp_error_at_f_$vt", sim; c=cs[2], kwargs...)
    if v == "XY"
        plot_integrand!("ncp_error_at_f_Y", sim; c=cs[3], kwargs...)
    end
    plot_integrand!("noise_error_at_f_$vt", sim; c=cs[4], kwargs...)
    vline!([sim.f_min_cost], ls=:dash, color=:black, label="cost cutoff freq.")
end

function plot_psd(f, p; normalize=true, kwargs...)
    pl = plot()
    plot_psd!(f, p; normalize=normalize, kwargs...)
    pl
end

function plot_psd!(f, p; normalize=true, kwargs...)
    if normalize
        p /= p[1]
    end
    plot!(f, p, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="Power"; kwargs...)
end

function plot_psd(psd::DSP.Periodograms.Periodogram; kwargs...)
    f, p = freq(psd)[2:end], power(psd)[2:end]
    plot_psd(f, p; kwargs...)
end

function plot_psd!(psd::DSP.Periodograms.Periodogram; kwargs...)
    f, p = freq(psd)[2:end], power(psd)[2:end]
    plot_psd!(f, p; kwargs...)
end

function plot_psd_p!(f, p; kwargs...)
    plot!(f, p, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)"; kwargs...)
end

function plot_timedomain_vs_analytic(sim, etffn, ol_timeseries, cl_timeseries)
    N = length(ol_timeseries)
	local fig = Figure(size=(1000,500))
	local axs = [
		Axis(fig[1,1], xlabel="Timesteps", ylabel="Time series (rad)"), Axis(fig[1,2], xscale=log10, yscale=log10, xlabel="Frequency (Hz)", ylabel="Power spectrum (rad²/Hz)"), 
		Axis(fig[1,3], xscale=log10, yscale=log10, xlabel="Frequency (Hz)", ylabel="$etffn |ETF|²")
	]
	lines!(axs[1], 1:N, ol_timeseries, label="Open loop")
	lines!(axs[1], 1:N, cl_timeseries, label="Closed loop at X")
	local freqs, psd_ol = genpsd(ol_timeseries, sim.f_loop)
	_, psd_cl = genpsd(cl_timeseries, sim.f_loop)
	lines!(axs[2], freqs, psd_ol, label="Open loop")
	lines!(axs[2], freqs, psd_cl, label="Closed loop at X")
	lines!(axs[3], freqs, psd_cl ./ psd_ol, label="Time-domain ETF", color=:teal)
	lines!(axs[3], sim.fr, abs2.(etffn.(sim.sT, Ref(sim))), label="Analytic ETF", color=:violetred)
	for ax in axs axislegend(ax, position=:rb) end
	fig
end

export nyquist_plot, ten_etf_plots, plot_integrand, plot_integrand!, plot_integrands, plot_psd_p!, plot_psd, plot_psd!, five_psd_plots, plot_timedomain_vs_analytic