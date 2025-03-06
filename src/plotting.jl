using Plots
using DSP
using Base.Meta: parse

function nyquist_plot!(sys; mark_gm_pm=true, label="Nyquist plot", kwargs...)
    success = :green
    nyquist_contour, gm, gm_point, pm, pm_point = nyquist_and_margins(sys)
    plot!(real(nyquist_contour), imag(nyquist_contour), xlim=(-1.1,1.1), ylim=(-1.1,1.1), aspect_ratio=:equal, legend=:outertopright, label=label; kwargs...)
    phasegrid = range(-π, π, length=500)
    xunit, yunit = cos.(phasegrid), sin.(phasegrid)
    if !is_stable(gm, pm)
        success = :red
    end
    if mark_gm_pm
        vline!([-1/2.5], ls=:dash, label="Gain margin cutoff", color=:grey)
        scatter!([real(gm_point)], [imag(gm_point)], label="Gain margin = $(round(gm, digits=2))", color=:grey)
        plot!([-2,0,-2], [-2,0,2], ls=:dash, label="Phase margin cutoff", color=4)
        if !isnothing(pm_point)
            scatter!([real(pm_point)], [imag(pm_point)], label="Phase margin = $(round(pm, digits=2))", color=4)
        end
    end
    plot!(xunit, yunit, ls=:dash, label=nothing, color=success)
end

function nyquist_plot(sys; kwargs...)
    p = plot()
    nyquist_plot!(sys; kwargs...)
    p
end

function ten_etf_plots(Cfast, Cslow, R, vk_atm, vk_ncp, f_noise_crossover; fr=exp10.(-2:0.01:log10(500.0)))
    plots_etf = []
    for v in ["X", "Y"]
        ne = eval(parse("notched_error_$v"))
        p_v = plot(legend=:bottomright, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="ETF", ylims=(1e-10, 1e2), title="$v error = $(round(ne(Cfast, Cslow, R, vk_atm, vk_ncp, f_noise_crossover), digits=3)) rad")
        for (fname, c, s) in zip(["phi_to_$v", "Lfast_to_$v", "Lslow_to_$v", "Nfast_to_$v", "Nslow_to_$v"], [1, 2, 2, 3, 3], [:solid, :solid, :dash, :solid, :dash])
            f = eval(parse(fname))
            plot!(fr, abs2.(f.(sT, Cfast, Cslow, 10)), label="|$fname|²", c=c, ls=s)
        end
        push!(plots_etf, p_v)
    end
    plot(plots_etf...)
end

function plot_integrand(errsource, Cfast, Cslow, R, vk_atm, vk_ncp, noise_normalization)
    p = plot()
    plot_integrand!(errsource, Cfast, Cslow, R, vk_atm, vk_ncp, noise_normalization)
    p
end

function plot_integrand!(errsource, Cfast, Cslow, R, vk_atm, vk_ncp, noise_normalization; fr=exp10.(-2:0.01:log10(500.0)), f_loop=1000.0, kwargs...)

    err_source_fn = eval(parse(errsource))
    plot!(fr, err_source_fn.(fr, Cfast, Cslow, R, Ref(vk_atm), Ref(vk_ncp), noise_normalization, f_loop), label=errsource, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="Closed-loop residual (rad/Hz)"; legend=:bottomleft, ylims=(1e-10, 1e2), kwargs...)
end

function plot_integrands(Cfast, Cslow, R, vk_atm, vk_ncp, noise_normalization)
    plot_integrand("atm_error_at_f_X", Cfast, Cslow, R, vk_atm, vk_ncp, noise_normalization)
    plot_integrand!("ncp_error_at_f_X", Cfast, Cslow, R, vk_atm, vk_ncp, noise_normalization)
    plot_integrand!("ncp_error_at_f_Y", Cfast, Cslow, R, vk_atm, vk_ncp, noise_normalization)
    plot_integrand!("atm_error_at_f_X", Cfast, Cslow, R, vk_atm, vk_ncp, noise_normalization)
    p
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

export nyquist_plot, nyquist_plot!, ten_etf_plots, plot_integrand, plot_integrand!, plot_integrands, plot_psd_p!, plot_psd, plot_psd!