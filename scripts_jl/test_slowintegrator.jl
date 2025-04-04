using multiwfs
using Plots

leak = 0.9
f_loop = 1000.0
R = 10
vk_atm = VonKarman(0.3 * 10.0 / 3.0, 0.25 * 0.1031^(-5/3))
vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25)
f_noise_crossover = 500.0
no_filter = ZPKFilter(0, 0, 1)
fast_controller = FilteredIntegrator(0.0, leak, no_filter, 1/f_loop)

begin
    slow_controller = FilteredIntegrator(1.7451, leak, no_filter, R/f_loop)
    sim = Simulation(f_loop, fast_controller, slow_controller, 10, vk_atm, vk_ncp, 500.0)

    A = [[0. 0. 0. 0. 0. ]
    [1. 0. 0. 0. 0. ]
    [0. 0. leak 0. 0. ]
    [0. 0. leak 0. 0. ]
    [0. 0. 0. 1. 0. ]]
    B = [0.; 0; 1; 1; 0]
    C = [-0.9 -0.1 0. 0.9 0.1]
    D = [1. 0 0 0 0]'
    Pw = hcat(1.0...)
    W = B * Pw * B'
    V = hcat(1e-6...)
    K = kalman_gain(A, C, W, V)
    Vv = [0 1. 0 0 -1]
    Q = Vv' * Vv
    Rlqg = zeros(1,1)
    L = lqr_gain(A, D, Q, Rlqg)
    lqg = LQG(A, D, C, K, L, R/f_loop)
    slow_controller_lqg = FilteredLQG(lqg, no_filter)

    sim_lqgic = Simulation(f_loop, fast_controller, slow_controller_lqg, R, vk_atm, vk_ncp, f_noise_crossover)
end

notched_error_X(sim)
notched_error_X(sim_lqgic)

# Then, a more useful way to show the better performance of LQG IC vs. IC 
# instead of matching stability margins is to add a noise component to the input PSD 
# and then find the optimal gain and leak for IC 
# And noise and alpha for LQG IC. 
# Adjust the noise level until the curve of gain vs. WFE is U-shaped 
# (i.e., the optimal gain is not due to stability limits). that case should clearly show a better LQG IC vs. IC. 

begin
    plot(sim.fr, abs2.(phi_to_X.(sim.sT, Ref(sim))), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="ETF", label="IC, bandwidth = $(round(zero_db_bandwidth(sim), digits=2)) Hz", legend=:topleft, size=(400,400))
    plot!(sim.fr, abs2.(phi_to_X.(sim.sT, Ref(sim_lqgic))), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="ETF", label="LQG-IC, bandwidth = $(round(zero_db_bandwidth(sim_lqgic), digits=2)) Hz", legend=:bottomright)
end