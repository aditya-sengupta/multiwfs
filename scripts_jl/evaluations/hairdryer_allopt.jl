using JLD2
using multiwfs

gridres = load("all_opt.jld2")
r0_ncp_vals = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 4.0, 5.0, 6.0]
f_crossover_vals = [50.0, 100.0, 200.0, 500.0]

begin
    f_crossover = 50.0
    x_errs = []
    for r0_ncp in r0_ncp_vals
        vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))
        push!(x_errs, notched_error_X(simgen_ichpf(gridres["($r0_ncp, 50.0)"][1], 0.0, 0.0, vk_ncp, f_crossover)))
    end
end

[grid_search_coarse_to_fine(gf -> simgen_ichpf(gf, 0.0, 0.0, vk_ncp, f_crossover), [[0.0, 1.0]])[1],
grid_search_coarse_to_fine((gf, gs) -> simgen_ichpf(gf, gs, 0.0, vk_ncp, f_crossover), [[0.0, 1.0], [0.0, 2.0]])[1],
grid_search_coarse_to_fine((gf, gs, fc) -> simgen_ichpf(gf, gs, fc, vk_ncp, f_crossover), [[0.0, 1.0], [0.0, 2.0], [0.0, 200.0]])[1],
grid_search_coarse_to_fine((af, lnf, gs, fc) -> simgen_lqgicfasthpf(af, lnf, gs, fc, vk_ncp, f_crossover), [[0.0, 1.0], [-10.0, 10.0], [0.0, 2.0], [0.0, 200.0]])[1],
grid_search_coarse_to_fine((as, lns, gf, fc) -> simgen_lqgicfasthpf(as, lns, gf, fc, vk_ncp, f_crossover), [[0.0, 1.0], [-10.0, 10.0], [0.0, 2.0], [0.0, 200.0]])[1],
grid_search_coarse_to_fine((as, lns, af, lnf, fc) -> simgen_lqgicfasthpf(as, lns, af, lnf, fc, vk_ncp, f_crossover), [[0.0, 1.0], [-10.0, 10.0], [0.0, 1.0], [-10.0, 10.0], [0.0, 200.0]])[1]
]