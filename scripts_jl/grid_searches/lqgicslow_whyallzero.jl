using multiwfs
using LinearAlgebra: BLAS
BLAS.set_num_threads(1)
r0_ncp = 0.6
f_crossover = 50.0

vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))

grid_search_coarse_to_fine((nna, lns, gf, fc) -> simgen_lqgicslowhpf(nna, lns, gf, fc, vk_ncp, f_crossover), [[1.0, 5.0], [-10.0, 10.0], [0.0, 1.0], [0.0, 40.0]])[1]