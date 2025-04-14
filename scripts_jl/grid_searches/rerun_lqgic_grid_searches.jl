using Pkg
Pkg.activate(".")
using multiwfs
using JLD2

r0_ncp_vals = [0.6]#, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 4.0, 5.0, 6.0]
f_crossover_vals = [50.0]#, 100.0, 200.0, 500.0]

all_opt = Dict()
for r0_ncp in r0_ncp_vals
    vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))
    for f_crossover in f_crossover_vals
        all_opt[string((r0_ncp,f_crossover))] = [
        grid_search_coarse_to_fine((las, lns, gf, fc) -> simgen_lqgicslowhpf(exp10(las), lns, gf, 15.0, vk_ncp, f_crossover), [[-1.0, 0.0], [-10.0, 10.0], [0.0, 2.0]])[1],
        ]
    end
end

save("lqgicslow_opt.jld2", all_opt)