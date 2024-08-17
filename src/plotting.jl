using Plots

function nyquist_plot(sys)
    success = :green
    nyquist_contour, gm, gm_point, pm, pm_point = nyquist_and_margins(sys)
    p = plot(real(nyquist_contour), imag(nyquist_contour), xlim=(-1.1,1.1), ylim=(-1.1,1.1), aspect_ratio=:equal, legend=:outertopright, label="Nyquist plot")
    phasegrid = range(-π, π, length=500)
    xunit, yunit = cos.(phasegrid), sin.(phasegrid)
    vline!([-1/2.5], ls=:dash, label="Gain margin cutoff", color=:grey)
    if !is_stable(gm, pm)
        success = :red
    end
    scatter!([real(gm_point)], [imag(gm_point)], label="Gain margin = $(round(gm, digits=2))", color=:grey)
    plot!([-2,0,-2], [-2,0,2], ls=:dash, label="Phase margin cutoff", color=4)
    plot!(xunit, yunit, ls=:dash, label=nothing, color=success)
    if !isnothing(pm_point)
        scatter!([real(pm_point)], [imag(pm_point)], label="Phase margin = $(round(pm, digits=2))", color=4)
    end
    p
end

function plot_psd_p!(f, p; kwargs...)
    plot!(f, p, xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="Power"; kwargs...)
end

export nyquist_plot, plot_psd_p!