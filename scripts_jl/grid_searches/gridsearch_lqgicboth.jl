using LinearAlgebra: BLAS
BLAS.set_num_threads(1)
using multiwfs
using Base: product
using Plots

r0 = 0.1031
r0_ncp = 1.0
vk_atm = VonKarman(0.3 * 10.0 / 3.0, 0.25 * r0^(-5/3))
vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))
f_loop = 1000.0
f_noise_crossover = 50.0
R = 10
leak = 0.9999

function sim_generator_lqgicboth(log_noise_slow, log_noise_fast, f_cutoff)
    A = [[0. 0. 0. 0. 0. ]
    [1. 0. 0. 0. 0. ]
    [0. 0. leak 0. 0. ]
    [0. 0. leak 0. 0. ]
    [0. 0. 0. 1. 0. ]]
    B = [0.; 0; 1; 1; 0]
    C = [-1+1/R -1/R 0. 1-1/R 1/R]
    D = [1. 0 0 0 0]'
    Pw = hcat(1.0...)
    W = B * Pw * B'
    V = hcat(exp10(log_noise_slow)...)
    K = kalman_gain(A, C, W, V)
    Vv = [0 1. 0 0 -1]
    Q = Vv' * Vv
    Rlqg = zeros(1,1)
    L = lqr_gain(A, D, Q, Rlqg)
    slow_controller_lqg = FilteredLQG(LQG(A, D, C, K, L, R/f_loop), no_filter)
    A_ar1 = [leak 0; 1 0]
    L = A_DM(2)
    Ã = block_diag(L, A_ar1)
    C̃ = [0 -1 0 1]
    D̃ = [1 0 0 0]' 
    B = [0; 0; 1; 0]
    Pw = hcat(1...)
    W = B * Pw * B'
    V = hcat(exp10(log_noise_fast)...)
    K̃ = kalman_gain(Ã, C̃, W, V)
    Vv = [0 -1 0 1]
    Q = Vv' * Vv
    Rlqg = zeros(1,1)
    L = lqr_gain(Ã, D̃, Q, Rlqg)
    fast_controller_lqg = FilteredLQG(LQG(Ã, D̃, C̃, K̃, L, 1/f_loop), ar1_filter(f_cutoff, f_loop, "high"))
    return Simulation(f_loop, fast_controller_lqg, slow_controller_lqg, R, vk_atm, vk_ncp, f_noise_crossover)
end

@time grid_search(
    sim_generator_lqgicboth,
    [-8.0:1.0:8.0, -8.0:1.0:8.0, 5.0:5.0:100.0]
)

@time grid_search(
    sim_generator_ic,
    [0.0:0.125:1.0, 0.0:0.0625:1.0, 5.0:5.0:100.0]
)

@time grid_search_coarse_to_fine(
    sim_generator_lqgicboth,
    [[-20.0, 5.0], [-8.0, 8.0], [0.0, 100.0]];
    search="parallel"
)

opt_sim_icboth = sim_generator_ic(1.42, 0.23, 20.6);
opt_sim_lqgicboth = sim_generator_lqgicboth(-7, 1, 16.0);

plot(ten_etf_plots(opt_sim_icboth)..., suptitle="IC slow and fast, optimized")
plot(ten_etf_plots(opt_sim_lqgicboth)..., suptitle="LQG-IC slow and fast, optimized")

plot(
    plot_integrands("XY", opt_sim_icboth, title="IC both, optimized, X error = $(round(notched_error_X(opt_sim_icboth), digits=3)) rad"),
    plot_integrands("XY", opt_sim_lqgicboth, title="LQG-IC both, optimized, X error = $(round(notched_error_X(opt_sim_lqgicboth), digits=3)) rad"),
    size=(1000,500)
)

