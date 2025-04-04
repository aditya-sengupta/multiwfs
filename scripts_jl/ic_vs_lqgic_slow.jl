using multiwfs
using Plots

r0 = 0.1031
r0_ncp = 1.0
vk_atm = VonKarman(0.3 * 10.0 / 3.0, 0.25 * r0^(-5/3))
vk_ncp = VonKarman(0.3 * 0.01 / 3.0, 0.25 * r0_ncp^(-5/3))

# fast_controller, slow_controller, R, vk_atm, vk_ncp, f_noise_crossover
f_loop = 1000.0
R = 10
leak = 0.999
fast_controller = FilteredIntegrator(0.0, leak, ZPKFilter(0, 0, 1), 1/f_loop)
slow_controller = FilteredIntegrator(0.43, leak, ZPKFilter(0, 0, 1), R/f_loop)
sim = Simulation(f_loop, fast_controller, slow_controller, R, vk_atm, vk_ncp, 500.0)

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
V = hcat(0.05...)
K = kalman_gain(A, C, W, V)
Vv = [0 1. 0 0 -1]
Q = Vv' * Vv
Rlqg = zeros(1,1)
L = lqr_gain(A, D, Q, Rlqg)
lqg = LQG(A, D, C, K, L, R/f_loop)
slow_controller_lqg = FilteredLQG(lqg, no_filter)

sim_lqgic = Simulation(f_loop, fast_controller, slow_controller_lqg, R, vk_atm, vk_ncp, 500.0)

begin
    plot(sim.fr, abs2.(phi_to_X.(sim.sT, Ref(sim))), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="ETF", label="IC, bandwidth = $(round(zero_db_bandwidth(sim), digits=2)) Hz", legend=:bottomright, size=(400,400))
    plot!(sim.fr, abs2.(phi_to_X.(sim.sT, Ref(sim_lqgic))), xscale=:log10, yscale=:log10, xlabel="Frequency (Hz)", ylabel="ETF", label="LQG-IC, bandwidth = $(round(zero_db_bandwidth(sim_lqgic), digits=2)) Hz")
end

plot_integrands("X", sim)
plot_integrands("X", sim_lqgic)