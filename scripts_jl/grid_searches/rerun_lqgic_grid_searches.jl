using Pkg
Pkg.activate(".")
using multiwfs
using JLD2

r0 = 0.1031
r0_ncp_vals = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 4.0, 5.0, 6.0]
f_crossover_vals = [50.0, 100.0, 200.0, 500.0]

all_opt = Dict()
for r0_ncp in r0_ncp_vals
    vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))
    for f_crossover in f_crossover_vals
        all_opt[string((r0_ncp,f_crossover))] = [
        grid_search_coarse_to_fine((as, lns, gf, fc) -> simgen_lqgicslowhpf(as, lns, gf, fc, vk_ncp, f_crossover), [[0.0, 1.0], [-10.0, 10.0], [0.0, 2.0], [0.0, 200.0]])[1],
        grid_search_coarse_to_fine((as, lns, af, lnf, fc) -> simgen_lqgicbothhpf(as, lns, af, lnf, fc, vk_ncp, f_crossover), [[0.0, 1.0], [-10.0, 10.0], [0.0, 1.0], [-10.0, 10.0], [0.0, 200.0]])[1]
        ]
    end
end

save("lqgic_opt.jld2", all_opt)