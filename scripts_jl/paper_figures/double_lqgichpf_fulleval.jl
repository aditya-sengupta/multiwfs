using multiwfs
using JLD2
using Plots

include("color_scheme.jl")
vk_ncp = VonKarman(0.001, 0.25)
f_noise_crossover = 50.0

optpars = [1.107421875, -10.0, 6.6, 1.0625, 15.84375]
sim = simgen_lqgicbothhpf(optpars..., vk_ncp, f_noise_crossover)
n1 = nyquist_plot(sim; size=(350,350), tickfontsize=20, legendfontsize=20, labelfontsize=20, label=nothing, color=general_color, lw=2)
ex, ey = ten_etf_plots(sim, [atm_color, ncp_color, ncp_color, noise_color, noise_color], titlefontsize=24, labelfontsize=20, tickfontsize=20, legendfontsize=20, lw=2)

parameter_names = ["α slow num. nines", "log₁₀(slow noise)", "α fast num. nines", "log₁₀(fast noise)", "Cutoff frequency"]
parameter_centers = [x for x in optpars]
fractional_change = 0.9:0.01:1.1
parameter_grids = [fractional_change .* x for x in parameter_centers]

# resolved parameter minimum grid
begin
    minplot = plot(legend=:topright, xlabel="Fractional change in parameter", labelfontsize=20, tickfontsize=20, legendfontsize=20)
    for (i, (parname, pargrid, c)) in enumerate(zip(parameter_names, parameter_grids, param_colors))
        Xerrs_parsweep = []
        for par in pargrid
            params_this = [i == j ? par : parameter_centers[j] for j in 1:5]
            sim_subopt = simgen_lqgicbothhpf(params_this..., vk_ncp, f_noise_crossover)
            push!(Xerrs_parsweep, notched_error_X(sim_subopt))
        end
        plot!(fractional_change, Xerrs_parsweep, ylabel="Error at X (rad)", c=c, label=parname)
        vline!([1.0], color=:black, ls=:dash, label=nothing)
    end
    minplot
end


integrands = plot_integrands("XY", sim, [atm_color, xsignal_color, ysignal_color, noise_color], labelfontsize=20, tickfontsize=20, legendfontsize=20, lw=2, xlims=(1.0, 500.0))
psds = five_psd_plots(sim, labelfontsize=20, tickfontsize=20, legendfontsize=20, cs=[atm_color, ncp_color, noise_color, xsignal_color, ysignal_color], lw=2, xlims=(1.0, 500.0))

lqgic_doublewfs_hpf = plot(n1, ex, ey, integrands, psds, minplot, size=(1200,800))
Plots.savefig(n1, "externalization/figures_tex/nyquist_lqgic_doublewfs_hpf.tex")
Plots.savefig(ex, "externalization/figures_tex/errx_lqgic_doublewfs_hpf.tex")
Plots.savefig(ey, "externalization/figures_tex/erry_lqgic_doublewfs_hpf.tex")
Plots.savefig(integrands, "externalization/figures_tex/integrands_lqgic_doublewfs_hpf.tex")
Plots.savefig(psds, "externalization/figures_tex/psds_lqgic_doublewfs_hpf.tex")
Plots.savefig(minplot, "externalization/figures_tex/minplot_lqgic_doublewfs_hpf.tex")