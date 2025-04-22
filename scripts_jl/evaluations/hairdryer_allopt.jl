using JLD2
using multiwfs
using Plots

gridres = load("data/all_opt.jld2")
gridres_y = load("data/all_opt_y.jld2")
gridres_lqgicslow = load("data/lqgicslow_opt.jld2")
gridres_lqgicboth = load("data/lqgicboth_opt.jld2")
r0_ncp_vals = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 4.0, 5.0, 6.0]
f_crossover_vals = [50.0, 100.0, 200.0, 500.0]

begin
    hairdryers = []
    for f_crossover in f_crossover_vals
        x_errs = [[], [], [], [], [], [], []]
        for r0_ncp in r0_ncp_vals
            vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))
            gf_oneint = gridres["($r0_ncp, $f_crossover)"][1][1]
            push!(x_errs[1], notched_error_X(simgen_ichpf(gridres["($r0_ncp, $f_crossover)"][1][1], 0.0, 0.0, vk_ncp, f_crossover)))
            push!(x_errs[2], notched_error_X(simgen_ichpf(gridres["($r0_ncp, $f_crossover)"][2]..., 0.0, vk_ncp, f_crossover)))
            push!(x_errs[3], notched_error_X(simgen_ichpf(gridres["($r0_ncp, $f_crossover)"][3]..., vk_ncp, f_crossover)))
            push!(x_errs[4], notched_error_X(simgen_ichpf(0.4, 1.4, 15.0, vk_ncp, f_crossover)))
            push!(x_errs[5], notched_error_X(simgen_lqgicfasthpf(gridres["($r0_ncp, $f_crossover)"][4]..., vk_ncp, f_crossover)))
            push!(x_errs[6], notched_error_X(simgen_lqgicslowhpf(gridres_lqgicslow["($r0_ncp, $f_crossover)"][1]..., vk_ncp, f_crossover)))
            push!(x_errs[7], notched_error_X(simgen_lqgicbothhpf(gridres_lqgicboth["($r0_ncp, $f_crossover)"][1]..., vk_ncp, f_crossover)))
        end
        push!(hairdryers, plot(r0_ncp_vals, x_errs, xscale=:log10, xticks=([0.6, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0], [0.6, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0]), label=["fast IC" "slow IC + fast IC" "slow IC + fast IC-HPF, opt" "slow IC + fast IC-HPF, Lisa" "slow IC + fast LQG-IC-HPF" "slow LQG-IC + fast IC-HPF" "slow LQG-IC + fast LQG-IC-HPF"], xlabel="r₀ NCP (m)", ylabel="X error (rad)", title="Noise crossover frequency = $(Int(f_crossover)) Hz", ylims=(0.5, 1.5), xtickfontsize=5, ytickfontsize=5, legendfontsize=5, titlefontsize=8, xlabelfontsize=8, ylabelfontsize=8, dpi=300, legend=(f_crossover == 500.0 ? :topright : nothing)))
    end
    plot(hairdryers...)
    # Plots.savefig("figures/hairdryer/hairdryer_ichpf.png")
end

begin
    f_cutoff_plots = []
    for f_crossover in f_crossover_vals
        f_cutoffs = [[], [], [], []]
        for r0_ncp in r0_ncp_vals
            push!(f_cutoffs[1], gridres["($r0_ncp, $f_crossover)"][3][3])
            push!(f_cutoffs[2], gridres["($r0_ncp, $f_crossover)"][4][3])
            push!(f_cutoffs[3], gridres_lqgicslow["($r0_ncp, $f_crossover)"][1][4])
            push!(f_cutoffs[4], gridres_lqgicboth["($r0_ncp, $f_crossover)"][1][5])
        end
        push!(f_cutoff_plots, plot(r0_ncp_vals, f_cutoffs, xscale=:log10, xticks=([0.6, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0], [0.6, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0]), label=["slow IC + fast IC-HPF" "slow IC + fast LQG-IC-HPF" "slow LQG-IC + fast IC-HPF" "slow LQG-IC + fast LQG-IC-HPF"], xlabel="r₀ NCP (m)", ylabel="Cutoff frequency (Hz)", title="Noise crossover frequency = $(Int(f_crossover)) Hz", xtickfontsize=5, ytickfontsize=5, legendfontsize=5, titlefontsize=8, xlabelfontsize=8, ylabelfontsize=8, dpi=300, legend=(f_crossover == 500.0 ? :topright : nothing), ylims=(0,50), color=[3 5 6 7]))
    end
    plot(f_cutoff_plots...)
end

begin
    hairdryers = []
    for f_crossover in f_crossover_vals
        x_errs = [[], []]
        for r0_ncp in r0_ncp_vals
            vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))
            push!(x_errs[1], notched_error_X(simgen_ichpf(0.4, 1.4, 15.0, vk_ncp, f_crossover)))
            push!(x_errs[2], notched_error_X(simgen_ichpf(gridres["($r0_ncp, $f_crossover)"][3]..., vk_ncp, f_crossover)))
        end
        push!(hairdryers, plot(r0_ncp_vals, x_errs, xscale=:log10, xticks=([0.6, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0], [0.6, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0]), label=["slow IC + fast IC-HPF, Lisa" "slow IC + fast IC-HPF, opt"], xlabel="r₀ NCP (m)", ylabel="X error (rad)", title="Noise crossover frequency = $(Int(f_crossover)) Hz", ylims=(0.5, 1.5), xtickfontsize=5, ytickfontsize=5, legendfontsize=5, titlefontsize=8, xlabelfontsize=8, ylabelfontsize=8, dpi=300))
    end
    plot(hairdryers...)
    # Plots.savefig("figures/hairdryer/hairdryer_ichpf.png")
end

begin
    hairdryers = []
    for f_crossover in f_crossover_vals
        x_errs = [[], []]
        for r0_ncp in r0_ncp_vals
            vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))
            push!(x_errs[1], notched_error_X(simgen_ichpf(gridres_y["($r0_ncp, $f_crossover)"][1][1], 0.0, 0.0, vk_ncp, f_crossover)))
            # push!(x_errs[2], notched_error_X(simgen_ichpf(gridres["($r0_ncp, $f_crossover)"][1][1], 0.0, 0.0, vk_ncp, f_crossover)))
            push!(
                x_errs[2], 
                min(
                    notched_error_X(simgen_ichpf(gridres["($r0_ncp, $f_crossover)"][1][1], 0.0, 0.0, vk_ncp, f_crossover)),
                    notched_error_X(simgen_ichpf(gridres["($r0_ncp, $f_crossover)"][2]..., 0.0, vk_ncp, f_crossover)),
                    notched_error_X(simgen_ichpf(gridres["($r0_ncp, $f_crossover)"][3]..., vk_ncp, f_crossover)),
                    notched_error_X(simgen_lqgicfasthpf(gridres["($r0_ncp, $f_crossover)"][4]..., vk_ncp, f_crossover))
                )
            )
        end
        push!(hairdryers, plot(r0_ncp_vals, x_errs, xscale=:log10, xticks=([0.6, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0], [0.6, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0]), label=["one-WFS optimal" "two-WFS optimal"], xlabel="r₀ NCP (m)", ylabel="X error (rad)", title="Noise crossover frequency = $(Int(f_crossover)) Hz", ylims=(0.5, 1.5), xtickfontsize=5, ytickfontsize=5, legendfontsize=5, titlefontsize=8, xlabelfontsize=8, ylabelfontsize=8, dpi=300))
    end
    plot(hairdryers...)
    Plots.savefig("figures/hairdryer/main_hairdryer.png")
end


begin
    hairdryers = []
    for f_crossover in f_crossover_vals
        bw = [[], [], [], [], []]
        for r0_ncp in r0_ncp_vals
            vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))
            gf_oneint = gridres["($r0_ncp, $f_crossover)"][1][1]
            push!(bw[1], zero_db_bandwidth(simgen_ichpf(gridres["($r0_ncp, $f_crossover)"][1][1], 0.0, 0.0, vk_ncp, f_crossover)))
            push!(bw[2], zero_db_bandwidth(simgen_ichpf(gridres["($r0_ncp, $f_crossover)"][2]..., 0.0, vk_ncp, f_crossover)))
            push!(bw[3], zero_db_bandwidth(simgen_ichpf(gridres["($r0_ncp, $f_crossover)"][3]..., vk_ncp, f_crossover)))
            push!(bw[4], zero_db_bandwidth(simgen_lqgicfasthpf(gridres["($r0_ncp, $f_crossover)"][4]..., vk_ncp, f_crossover)))
            push!(bw[5], zero_db_bandwidth(simgen_lqgicslowhpf(gridres_lqgicslow["($r0_ncp, $f_crossover)"][1]..., vk_ncp, f_crossover)))
        end
        push!(hairdryers, plot(r0_ncp_vals, bw, xscale=:log10, xticks=([0.6, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0], [0.6, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0]), label=["fast IC" "slow IC + fast IC" "slow IC + fast IC-HPF" "slow IC + fast LQG-IC-HPF" "slow LQG-IC + fast IC-HPF"], xlabel="r₀ NCP (m)", ylabel="Bandwidth (Hz)", title="Noise crossover frequency = $(Int(f_crossover)) Hz", xtickfontsize=5, ytickfontsize=5, legend=nothing, legendfontsize=5, titlefontsize=8, xlabelfontsize=8, ylabelfontsize=8, dpi=300, ylims=(30, 90)))
    end
    plot(hairdryers...)
end

"""[grid_search_coarse_to_fine(gf -> simgen_ichpf(gf, 0.0, 0.0, vk_ncp, f_crossover), [[0.0, 1.0]])[1],
grid_search_coarse_to_fine((gf, gs) -> simgen_ichpf(gf, gs, 0.0, vk_ncp, f_crossover), [[0.0, 1.0], [0.0, 2.0]])[1],
grid_search_coarse_to_fine((gf, gs, fc) -> simgen_ichpf(gf, gs, fc, vk_ncp, f_crossover), [[0.0, 1.0], [0.0, 2.0], [0.0, 200.0]])[1],
grid_search_coarse_to_fine((af, lnf, gs, fc) -> simgen_lqgicfasthpf(af, lnf, gs, fc, vk_ncp, f_crossover), [[0.0, 1.0], [-10.0, 10.0], [0.0, 2.0], [0.0, 200.0]])[1],
grid_search_coarse_to_fine((as, lns, gf, fc) -> simgen_lqgicfasthpf(as, lns, gf, fc, vk_ncp, f_crossover), [[0.0, 1.0], [-10.0, 10.0], [0.0, 2.0], [0.0, 200.0]])[1],
grid_search_coarse_to_fine((as, lns, af, lnf, fc) -> simgen_lqgicfasthpf(as, lns, af, lnf, fc, vk_ncp, f_crossover), [[0.0, 1.0], [-10.0, 10.0], [0.0, 1.0], [-10.0, 10.0], [0.0, 200.0]])[1]
]"""