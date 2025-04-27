using LinearAlgebra: BLAS
BLAS.set_num_threads(1)
using multiwfs
using Base: product
using Plots
using JLD2

r0 = 0.1031
r0_ncp = 1.0
vk_atm = VonKarman(0.3 * 10.0 / 3.0, 0.25 * r0^(-5/3))
vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))
f_noise_crossover = 50.0

load("data/lqgicboth_opt.jld2")["(1.0, 50.0)"]

grid_search_coarse_to_fine(
    (p1, p2, p3, p4) -> simgen_lqgicbothhpf(p1, -10.0, p2, p3, p4, vk_ncp, f_noise_crossover),
    [[1.0, 2.0], [4.0, 6.0], [0.9, 1.1], [10.0, 30.0]]; npoints_per_parameter=5
)

grid_search_coarse_to_fine(
    (p1, p2, p3, p4) -> simgen_lqgicbothhpf(p1, -10.0, p2, p3, p4, vk_ncp, f_noise_crossover),
    [[0.984375, 1.203125], [5.4, 6.6], [0.95625, 1.16875], [14.625, 17.875]]; npoints_per_parameter=5
)

