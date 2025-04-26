using Pkg
Pkg.activate(".")
using multiwfs
using JLD2

r0 = 0.1031
vk_atm = VonKarman(0.3 * 10.0 / 3.0, 0.25 * r0^(-5/3))
r0_ncp_vals = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 4.0, 5.0, 6.0]
f_crossover_vals = [50.0, 100.0, 200.0, 500.0]
f_loop = 1000.0
leak = 0.999
R = 10

all_opt = Dict()
for r0_ncp in r0_ncp_vals
    vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))
    for f_crossover in f_crossover_vals
        println(string((r0_ncp,f_crossover)))
        all_opt[string((r0_ncp,f_crossover))] = [grid_search_coarse_to_fine(gf -> simgen_ichpf(gf, 0.0, 0.0, vk_ncp, f_crossover), [[0.0, 1.0]], search="parallel")[1],
        grid_search_coarse_to_fine((gf, gs) -> simgen_ichpf(gf, gs, 0.0, vk_ncp, f_crossover), [[0.0, 1.0], [0.0, 2.0]], search="parallel")[1],
        grid_search_coarse_to_fine((gf, gs, fc) -> simgen_ichpf(gf, gs, fc, vk_ncp, f_crossover), [[0.0, 1.0], [0.0, 2.0], [0.0, 200.0]], search="parallel")[1],
        ]
    end
end

save("ic_opt_for_paper.jld2", all_opt)