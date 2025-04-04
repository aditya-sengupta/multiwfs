using multiwfs
using Base: product

r0 = 0.1031
r0_ncp = 1.0
vk_atm = VonKarman(0.3 * 10.0 / 3.0, 0.25 * r0^(-5/3))
vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))
f_loop = 1000.0
f_noise_crossover = 500.0
R = 10
leak = 0.999

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

s1 = sim_generator_lqgicboth(-8, -8, 1.0)
s2 = sim_generator_lqgicboth(-8, -7, 1.0)

@time notched_error_X(s1);
@time notched_error_X(s2);

@time grid_search(
    sim_generator_lqgicboth,
    [-8.0:1:-1, -8.0:1:1, 5:5:40.]
)
