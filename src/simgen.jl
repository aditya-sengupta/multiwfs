function simgen_ichpf(gain_fast, gain_slow, f_cutoff, vk_ncp, f_noise_crossover; leak=0.99, R=10, f_loop=1000.0, vk_atm=VonKarman(0.3 * 10.0 / 3.0, 0.25 * (0.1031)^(-5/3)))
    fast_controller = FilteredIntegrator(gain_fast, leak, ar1_filter(f_cutoff, f_loop, "high"), 1/f_loop)
	slow_controller = FilteredIntegrator(gain_slow, leak, no_filter, R/f_loop)
	return Simulation(f_loop, fast_controller, slow_controller, R, vk_atm, vk_ncp, f_noise_crossover)
end

function simgen_lqgicfasthpf(alpha_fast, log_noise_fast, gain_slow, f_cutoff, vk_ncp, f_noise_crossover; leak=0.999, R=10, f_loop=1000.0, vk_atm=VonKarman(0.3 * 10.0 / 3.0, 0.25 * (0.1031)^(-5/3)))
    A_ar1 = [alpha_fast 0; 1 0]
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
    fast_controller = FilteredLQG(LQG(Ã, D̃, C̃, K̃, L, 1/f_loop), ar1_filter(f_cutoff, f_loop, "high"))
    slow_controller = FilteredIntegrator(gain_slow, leak, no_filter, R/f_loop)
	return Simulation(f_loop, fast_controller, slow_controller, R, vk_atm, vk_ncp, f_noise_crossover)
end

function simgen_lqgicslowhpf(num_nines_alpha_slow, log_noise_slow, gain_fast, f_cutoff, vk_ncp, f_noise_crossover; leak=0.999, R=10, f_loop=1000.0, vk_atm=VonKarman(0.3 * 10.0 / 3.0, 0.25 * (0.1031)^(-5/3)))
    alpha_slow = 1 - exp10(-num_nines_alpha_slow)
    A = [[0. 0. 0. 0. 0. ]
    [1. 0. 0. 0. 0. ]
    [0. 0. alpha_slow 0. 0. ]
    [0. 0. alpha_slow 0. 0. ]
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
    fast_controller = FilteredIntegrator(gain_fast, leak, ar1_filter(f_cutoff, f_loop, "high"), 1/f_loop)
    slow_controller = FilteredLQG(LQG(A, D, C, K, L, R/f_loop), no_filter)
    return Simulation(f_loop, fast_controller, slow_controller, R, vk_atm, vk_ncp, f_noise_crossover)
end

function simgen_lqgicbothhpf(num_nines_alpha_slow, log_noise_slow, num_nines_alpha_fast, log_noise_fast, f_cutoff, vk_ncp, f_noise_crossover; R=10, f_loop=1000.0, vk_atm=VonKarman(0.3 * 10.0 / 3.0, 0.25 * (0.1031)^(-5/3)))
    alpha_slow = 1 - exp10(-num_nines_alpha_slow)
    alpha_fast = 1 - exp10(-num_nines_alpha_fast)
    A_ar1 = [alpha_fast 0; 1 0]
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
    fast_controller = FilteredLQG(LQG(Ã, D̃, C̃, K̃, L, 1/f_loop), ar1_filter(f_cutoff, f_loop, "high"))
    A = [[0. 0. 0. 0. 0. ]
    [1. 0. 0. 0. 0. ]
    [0. 0. alpha_slow 0. 0. ]
    [0. 0. alpha_slow 0. 0. ]
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
    slow_controller = FilteredLQG(LQG(A, D, C, K, L, R/f_loop), no_filter)
    return Simulation(f_loop, fast_controller, slow_controller, R, vk_atm, vk_ncp, f_noise_crossover)
end

export simgen_ichpf, simgen_lqgicbothhpf, simgen_lqgicfasthpf, simgen_lqgicslowhpf