using JLD2
using multiwfs
using Plots

gridres = load("all_opt.jld2")
r0_ncp_vals = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 4.0, 5.0, 6.0]
f_crossover_vals = [50.0, 100.0, 200.0, 500.0]

begin
    f_crossover = 200.0
    x_errs = [[], [], [], []]
    for r0_ncp in r0_ncp_vals
        vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))
        gf_oneint = gridres["($r0_ncp, $f_crossover)"][1][1]
        push!(x_errs[1], notched_error_X(simgen_ichpf(gridres["($r0_ncp, $f_crossover)"][1][1], 0.0, 0.0, vk_ncp, f_crossover)))
        push!(x_errs[2], notched_error_X(simgen_ichpf(gridres["($r0_ncp, $f_crossover)"][2]..., 0.0, vk_ncp, f_crossover)))
        push!(x_errs[3], notched_error_X(simgen_ichpf(gridres["($r0_ncp, $f_crossover)"][3]..., vk_ncp, f_crossover)))
        push!(x_errs[4], notched_error_X(simgen_lqgicfasthpf(gridres["($r0_ncp, $f_crossover)"][4]..., vk_ncp, f_crossover)))
    end
    plot(r0_ncp_vals, x_errs, xscale=:log10, xticks=(r0_ncp_vals, r0_ncp_vals), label=["fast IC" "slow IC + fast IC" "slow IC + fast IC-HPF" "slow IC + fast LQG-IC-HPF"], xlabel="r₀ NCP (m)", ylabel="X error (rad)", title="Hairdryer plot, crossover frequency = $f_crossover Hz")
end



[grid_search_coarse_to_fine(gf -> simgen_ichpf(gf, 0.0, 0.0, vk_ncp, f_crossover), [[0.0, 1.0]])[1],
grid_search_coarse_to_fine((gf, gs) -> simgen_ichpf(gf, gs, 0.0, vk_ncp, f_crossover), [[0.0, 1.0], [0.0, 2.0]])[1],
grid_search_coarse_to_fine((gf, gs, fc) -> simgen_ichpf(gf, gs, fc, vk_ncp, f_crossover), [[0.0, 1.0], [0.0, 2.0], [0.0, 200.0]])[1],
grid_search_coarse_to_fine((af, lnf, gs, fc) -> simgen_lqgicfasthpf(af, lnf, gs, fc, vk_ncp, f_crossover), [[0.0, 1.0], [-10.0, 10.0], [0.0, 2.0], [0.0, 200.0]])[1],
grid_search_coarse_to_fine((as, lns, gf, fc) -> simgen_lqgicfasthpf(as, lns, gf, fc, vk_ncp, f_crossover), [[0.0, 1.0], [-10.0, 10.0], [0.0, 2.0], [0.0, 200.0]])[1],
grid_search_coarse_to_fine((as, lns, af, lnf, fc) -> simgen_lqgicfasthpf(as, lns, af, lnf, fc, vk_ncp, f_crossover), [[0.0, 1.0], [-10.0, 10.0], [0.0, 1.0], [-10.0, 10.0], [0.0, 200.0]])[1]
]