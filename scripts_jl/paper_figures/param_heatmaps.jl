using multiwfs
using JLD2
using Base: product
pgfplotsx()

all_icopt = load("data/ic_opt_for_paper.jld2")
all_lqgicboth_opt = load("data/lqgicboth_opt.jld2")

r0_ncp_vals = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 4.0, 5.0, 6.0]
f_crossover_vals = [50.0, 100.0, 200.0, 500.0]

function heatmap_from_pars(k1,k2;kwargs...)
    pars = zeros(length(f_crossover_vals), length(r0_ncp_vals))
    for (i, r0_ncp) in enumerate(r0_ncp_vals) for (j, f_crossover) in enumerate(f_crossover_vals)
        key = "($r0_ncp, $f_crossover)"
        if k1 < 4
            pars[j,i] = all_icopt[key][k1][k2]
        else
            pars[j,i] = all_lqgicboth_opt[key][1][k2]
        end
    end end

    heatmap(
        pars,
        xticks=(eachindex(r0_ncp_vals),r0_ncp_vals), 
        yticks=(eachindex(f_crossover_vals), f_crossover_vals),
        xlabel="NCP r₀ (m)",
        ylabel="Crossover frequency (Hz)",
        grid=false; labelfontsize=24, titlefontsize=30, kwargs...
    )
end

heatmap_from_pars(1,1, title="Single-IC optimal fast gain", clims=(0,1))
Plots.savefig("externalization/figures_tex/singlewfs_fastgain_heatmap.tex")

heatmap_from_pars(2,1, title="Double-IC optimal fast gain", clims=(0,1))
Plots.savefig("externalization/figures_tex/doublewfs_fastgain_heatmap.tex")

heatmap_from_pars(2,2, title="Double-IC optimal slow gain",clims=(0,2))
Plots.savefig("externalization/figures_tex/doublewfs_slowgain_heatmap.tex")

heatmap_from_pars(3,1, title="Double-IC-HPF optimal fast gain",clims=(0,1))
Plots.savefig("externalization/figures_tex/doublewfshpf_fastgain_heatmap.tex")

heatmap_from_pars(3,2, title="Double-IC-HPF optimal slow gain", clims=(0,2))
Plots.savefig("externalization/figures_tex/doublewfshpf_slowgain_heatmap.tex")

heatmap_from_pars(3,3, title="Double-IC-HPF optimal cutoff frequency",clims=(0,50))
Plots.savefig("externalization/figures_tex/doublewfshpf_fcutoff_heatmap.tex")

heatmap_from_pars(4,1,clims=(1,5), title="LQG-IC-HPF optimal slow -log₁₀(1-α)")
Plots.savefig("externalization/figures_tex/lqgichpf_slownumnines_heatmap.tex")

heatmap_from_pars(4,2,clims=(-10,10), title="LQG-IC-HPF optimal slow log₁₀(noise)")
Plots.savefig("externalization/figures_tex/lqgichpf_slownoise_heatmap.tex")

heatmap_from_pars(4,3,clims=(1,5), title="LQG-IC-HPF optimal fast -log₁₀(1-α)")
Plots.savefig("externalization/figures_tex/lqgichpf_fastnumnines_heatmap.tex")

heatmap_from_pars(4,4,clims=(-10,10),title="LQG-IC-HPF optimal fast log₁₀(noise)")
Plots.savefig("externalization/figures_tex/lqgichpf_fastnoise_heatmap.tex")

heatmap_from_pars(4,5,clims=(0,50),title="LQG-IC-HPF optimal cutoff frequency")
Plots.savefig("externalization/figures_tex/lqgichpf_fcutoff_heatmap.tex")
