using multiwfs
using JLD2
using Plots
pgfplotsx()

vk_ncp = VonKarman(0.001, 0.25)
f_noise_crossover = 50.0

gain_fast = load("data/ic_opt_for_paper.jld2")["(1.0, 50.0)"][1][1]
sim = simgen_ichpf(gain_fast, 0.0, 0.0, VonKarman(0.001, 0.25), 50.0)
n1 = nyquist_plot(sim; size=(350,350))
ex, ey = ten_etf_plots(sim)

parameter_names = ["Fast-WFS gain"]
parameter_centers = [gain_fast]
parameter_grids = [(gain_fast-0.1):0.001:(gain_fast+0.1)]
# resolved parameter minimum grid
min_plots = []
for (i, (parname, pargrid)) in enumerate(zip(parameter_names, parameter_grids))
    Xerrs_parsweep = []
    for par in pargrid
        params_this = [par]
        sim_subopt = simgen_ichpf(params_this..., 0.0, 0.0, vk_ncp, f_noise_crossover)
        push!(Xerrs_parsweep, notched_error_X(sim_subopt))
    end
    local p = plot(pargrid, Xerrs_parsweep, xlabel=parname, ylabel="Error at X (rad)", legend=nothing)
    vline!([parameter_centers[i]], color=:black, ls=:dash)
    push!(min_plots, p)
end

integrands = plot_integrands("XY", sim)
psds = five_psd_plots(sim)

min_plots[1]

singlewfs = plot(n1, ex, ey, integrands, psds, min_plots[1], size=(1200,800))
Plots.savefig(n1, "externalization/figures_tex/nyquist_singlewfs.tex")
Plots.savefig(ex, "externalization/figures_tex/errx_singlewfs.tex")
Plots.savefig(ey, "externalization/figures_tex/erry_singlewfs.tex")
Plots.savefig(integrands, "externalization/figures_tex/integrands_singlewfs.tex")
Plots.savefig(psds, "externalization/figures_tex/psds_singlewfs.tex")
Plots.savefig(min_plots[1], "externalization/figures_tex/minplot_singlewfs.tex")