using multiwfs
using JLD2
using Plots
pgfplotsx()

include("color_scheme.jl")
vk_ncp = VonKarman(0.001, 0.25)
f_noise_crossover = 50.0

gain_fast, gain_slow = load("data/ic_opt_for_paper.jld2")["(1.0, 50.0)"][2]
sim = simgen_ichpf(gain_fast, gain_slow, 0.0, VonKarman(0.001, 0.25), 50.0)
n1 = nyquist_plot(sim; size=(350,350), tickfontsize=20, legendfontsize=20, label=nothing, color=general_color, lw=2)
ex, ey = ten_etf_plots(sim, [atm_color, ncp_color, ncp_color, noise_color, noise_color], titlefontsize=24, labelfontsize=20, tickfontsize=20, legendfontsize=20, lw=2)

parameter_names = ["Fast-WFS gain", "Slow-WFS gain"]
parameter_centers = [gain_fast, gain_slow]
fractional_change = 0.9:0.001:1.1
parameter_grids = [(0.9:0.001:1.1) .* x for x in [gain_fast, gain_slow]]

# resolved parameter minimum grid
begin
    minplot = plot(legend=:topright, xlabel="Fractional change in parameter", labelfontsize=20, tickfontsize=20, legendfontsize=20)
    for (i, (parname, pargrid, c)) in enumerate(zip(parameter_names, parameter_grids, param_colors))
        Xerrs_parsweep = []
        for par in pargrid
            params_this = [i == j ? par : parameter_centers[j] for j in 1:2]
            sim_subopt = simgen_ichpf(params_this..., 0.0, vk_ncp, f_noise_crossover)
            push!(Xerrs_parsweep, notched_error_X(sim_subopt))
        end
        plot!(fractional_change, Xerrs_parsweep, ylabel="Error at X (rad)", label=parname, c=c, lw=2)
        vline!([1.0], color=:black, ls=:dash, label=nothing)
    end
    minplot
end

integrands = plot_integrands("XY", sim, [atm_color, xsignal_color, ysignal_color, noise_color], labelfontsize=20, tickfontsize=20, legendfontsize=20, lw=2, xlims=(1.0, 500.0))
psds = five_psd_plots(sim, labelfontsize=20, tickfontsize=20, legendfontsize=20, cs=[atm_color, ncp_color, noise_color, xsignal_color, ysignal_color], lw=2, xlims=(1.0, 500.0))
doublewfs = plot(n1, ex, ey, integrands, psds, minplot, size=(1200,800))

Plots.savefig(n1, "externalization/figures_tex/nyquist_doublewfs.tex")
Plots.savefig(ex, "externalization/figures_tex/errx_doublewfs.tex")
Plots.savefig(ey, "externalization/figures_tex/erry_doublewfs.tex")
Plots.savefig(integrands, "externalization/figures_tex/integrands_doublewfs.tex")
Plots.savefig(psds, "externalization/figures_tex/psds_doublewfs.tex")
Plots.savefig(minplot, "externalization/figures_tex/minplot_doublewfs.tex")